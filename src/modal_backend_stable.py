"""
This script defines a Modal application that serves a FastAPI backend for running
Nardini, a tool for analyzing Intrinsically Disordered Regions (IDRs) in proteins.

The backend uses Modal's job queue pattern with:
- A lightweight FastAPI web server for handling uploads and spawning jobs
- A heavy Nardini worker for processing sequences

The backend exposes a `/upload_fasta` endpoint that accepts a FASTA file,
spawns a Nardini analysis job, and returns a run_id for tracking.

Deployment:
-----------
modal deploy modal_backend.py

Example usage with curl:
-----------------------
# Replace with your deployed Modal app URL.
curl -X POST "YOUR_MODAL_APP_URL/upload_fasta" \
-F "file=@data/fasta_inputs/Halophile-pHtolerant-yeast-first16.fasta"
"""

import json
import logging
import os
import tempfile
import time
import zipfile
from pathlib import Path
from uuid import uuid4

import modal

# ---------------------- Environment configuration ---------------------- #
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()
logging.basicConfig(
    level=LOG_LEVEL,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("nardini_backend")

# Configurable parameters with env-variable overrides
VOLUME_DIR = Path(os.getenv("VOLUME_DIR", "/data"))
VOLUME_NAME = os.getenv("VOLUME_NAME", "nardini_volume_stable")
TIMEOUT_SECONDS = int(os.getenv("TIMEOUT_SECONDS", "43200"))  # default 12h
MAX_UPLOAD_MB = int(os.getenv("MAX_UPLOAD_MB", "10"))
MAX_FILE_SIZE = MAX_UPLOAD_MB * 1024 * 1024  # bytes
WORKER_COUNT = int(os.getenv("WORKER_COUNT", "16"))
# ---------------------------------------------------------------------- #

# Lightweight web image for FastAPI
web_image = modal.Image.debian_slim().pip_install(
    "fastapi[standard]==0.115.13",
    "python-multipart==0.0.20",
)

# Heavy processing image for Nardini
nardini_image = (
    modal.Image.debian_slim(python_version="3.9")
    .uv_pip_install(
        "pydantic==2.10.6",
        "numpy==1.19.5",
    )
    .uv_pip_install("nardini==1.1.1")
)

# Create a Modal application with shared volume
app = modal.App("nardini_backend_stable")
vol = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

# Web image imports (lightweight)
with web_image.imports():
    from fastapi import FastAPI, File, HTTPException, UploadFile
    from fastapi.responses import Response

# Nardini image imports (heavy)
with nardini_image.imports():
    import functools
    from concurrent.futures import ProcessPoolExecutor

    from nardini.constants import (
        DEFAULT_PREFIX_NAME,
        DEFAULT_RANDOM_SEED,
        NUM_SCRAMBLED_SEQUENCES,
        TYPEALL,
    )
    from nardini.score_and_plot import calculate_zscore_and_plot
    from nardini.utils import read_sequences_from_filename, set_random_seed


### NOTE Functions for stitching .zip results from each invocation together into final .zip
# helper fxn for mergeZips
def addtoSeqTSV(masterSeqTSV, tsvContent):
    import re

    lines = tsvContent.splitlines()
    for line in lines:
        line = line.strip()
        if not line:  # Skip empty lines
            continue
        # Split on 2+ consecutive spaces (the actual format used by Nardini)
        lineList = re.split(r" {2,}", line)
        # Skip header lines (check if first column is 'ID')
        if lineList and lineList[0] != "ID":
            with open(masterSeqTSV, "a") as f:
                # Join with tabs for proper TSV format in output
                f.write("\t".join(lineList) + "\n")
    return


def _get_unique_filename(original_filename, used_filenames, zip_index):
    """Generate a unique filename by adding a suffix if the original is already used"""
    if original_filename not in used_filenames:
        return original_filename

    # Extract filename and extension
    if "." in original_filename:
        name, ext = original_filename.rsplit(".", 1)
        ext = "." + ext
    else:
        name = original_filename
        ext = ""

    # Try adding zip index first
    candidate = f"{name}_zip{zip_index}{ext}"
    if candidate not in used_filenames:
        return candidate

    # If still duplicate, add a counter
    counter = 1
    while True:
        candidate = f"{name}_zip{zip_index}_{counter}{ext}"
        if candidate not in used_filenames:
            return candidate
        counter += 1


def mergeZips(zipList, destination_filepath):
    """Merge multiple zip files into a single zip, combining sequences.tsv files into master_sequences.tsv"""
    if not zipList:
        raise ValueError("No zip files provided to merge")

    destination_filepath = str(destination_filepath)
    if not destination_filepath.endswith(".zip"):
        destination_filepath = destination_filepath + ".zip"

    # Use a temporary file for the master TSV to avoid working directory issues
    with tempfile.NamedTemporaryFile(
        mode="w+", suffix=".tsv", delete=False
    ) as temp_tsv:
        masterTSVPath = temp_tsv.name
        # Write header
        header = [
            "ID",
            "original_seq",
            "most_similar_seq",
            "sum_abs_zscore_original_seq",
            "sum_abs_zscore_scrambled_seq",
        ]
        temp_tsv.write("\t".join(header) + "\n")

    try:
        # Create a new zip file with the collected files
        with zipfile.ZipFile(destination_filepath, "w") as merged_zip:
            used_filenames = (
                set()
            )  # Track filenames already added to prevent duplicates

            for zip_index, zipPath in enumerate(zipList):
                with zipfile.ZipFile(zipPath, "r") as zip_ref:
                    for file_info in zip_ref.infolist():
                        fileName = file_info.filename
                        if fileName != "sequences.tsv":
                            # Handle potential duplicate filenames
                            unique_filename = _get_unique_filename(
                                fileName, used_filenames, zip_index
                            )
                            used_filenames.add(unique_filename)

                            # Read the file content from the original zip and write to the new zip
                            with zip_ref.open(file_info) as source_file:
                                merged_zip.writestr(unique_filename, source_file.read())
                        else:
                            # Read the content of sequences.tsv from within the zip
                            with zip_ref.open(file_info) as seq_tsv_file:
                                seqTSVContent = seq_tsv_file.read().decode(
                                    "utf-8"
                                )  # Decode bytes to string
                                logger.info(f"Processing sequences.tsv from {zipPath}")
                                addtoSeqTSV(
                                    masterTSVPath, seqTSVContent
                                )  # Pass content to helper function

            # Write the content of the master_sequences.tsv into the zip while it's open
            with open(masterTSVPath, "r") as master_tsv_file:
                master_tsv_content = master_tsv_file.read()
                merged_zip.writestr("master_sequences.tsv", master_tsv_content)

    finally:
        # Clean up the temporary file
        if os.path.exists(masterTSVPath):
            os.unlink(masterTSVPath)

    if not Path(destination_filepath).exists():
        raise ValueError("Zip creation failed!")
    return destination_filepath


###
def process_single_sequence(
    seq_record,
    amino_acid_groupings,
    num_scrambled_sequences,
    random_seed,
    output_dir,
):
    """Process a single sequence with Nardini analysis and return the created zip path.

    Each call is executed in its own temporary directory so that concurrently running
    processes cannot clash while writing output files.  The produced zip archive is
    moved to *output_dir* and the final absolute path is returned to the caller.
    """

    import shutil  # Local import to ensure availability inside worker process

    # Create an isolated working directory for this worker process
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Temporarily switch the current working directory so that Nardini writes
        # inside the isolated directory.
        old_cwd = os.getcwd()
        os.chdir(tmp_path)

        try:
            # Run the calculation
            calculate_zscore_and_plot(
                [seq_record],  # Preserve original API (expects list)
                amino_acid_groupings,
                num_scrambled_sequences,
                random_seed,
            )

            # Locate the zip produced by Nardini (there should be exactly one)
            produced_zips = list(tmp_path.glob("nardini-data-*.zip"))
            if len(produced_zips) != 1:
                raise RuntimeError(
                    f"Expected one zip for sequence {seq_record.id}, found {len(produced_zips)}"
                )

            produced_zip = produced_zips[0]

            # Build a deterministic destination name to avoid collisions
            # Use a hash of the sequence content to create a unique filename
            import hashlib

            sequence_hash = hashlib.md5(str(seq_record.seq).encode()).hexdigest()[:8]
            dest_zip_name = f"nardini-{seq_record.id}-{sequence_hash}.zip"
            dest_path = Path(output_dir) / dest_zip_name

            # Move the zip file to the shared output directory
            shutil.move(str(produced_zip), dest_path)
            logger.info(
                f"Moved {produced_zip.name} to {dest_path.name} for sequence {seq_record.id}"
            )

        finally:
            # Restore the original cwd no matter what
            os.chdir(old_cwd)

    # Return the absolute path so the parent can build the mapping
    return str(dest_path)


# 24 hour timeout
@app.function(
    image=nardini_image,
    volumes={str(VOLUME_DIR): vol},
    # cpu=4.0,
    # memory=4096,
    timeout=TIMEOUT_SECONDS,
)
def process_nardini_job(sequences_data: dict, run_id: str) -> dict:
    """Heavy Nardini processing job that runs in parallel workers."""
    start_time = time.time()

    try:
        # Parse FASTA content (now we have access to nardini.utils in this image)
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_file:
            temp_file.write(sequences_data["file_content"])
            temp_file_path = temp_file.name

        try:
            parsed_sequences = read_sequences_from_filename(
                temp_file_path, DEFAULT_PREFIX_NAME, verbose=False
            )
        finally:
            if os.path.exists(temp_file_path):
                os.unlink(temp_file_path)

        if not parsed_sequences:
            raise ValueError("No valid sequences found in FASTA file")

        # Log parsed sequences info
        logger.info(f"Parsed {len(parsed_sequences)} sequences from FASTA file")

        # Read existing analyzed sequences from the modal volume
        existing_sequence_map = {}
        seq_zip_maps_dir = VOLUME_DIR / "seq_zip_maps"
        if not seq_zip_maps_dir.exists():
            seq_zip_maps_dir.mkdir(parents=True, exist_ok=True)

        cache_files = [f for f in os.listdir(seq_zip_maps_dir) if f.endswith(".json")]
        logger.info(f"Found {len(cache_files)} cache files")

        for file in cache_files:
            with open(seq_zip_maps_dir / file, "r") as f:
                # get the keys
                file_dict = json.load(f)
                for key, value in file_dict.items():
                    existing_sequence_map[key] = value["zip_filename"]
        existing_sequences = set(existing_sequence_map.keys())
        logger.info(f"Total cached sequences loaded: {len(existing_sequences)}")

        # Create an in-progress file in the modal volume
        in_progress_filename = f"{run_id}_in_progress.json"
        # Track progress for each requested sequence plus global run status.
        progress_dict = {
            "_status": "in_progress",  # overall run status
            "merged_zip_filename": None,  # will be filled in once complete
        }

        # Populate per-sequence status (None means still pending).
        cached_count = 0
        novel_count = 0

        for seq in parsed_sequences:
            sequence_string = str(seq.seq)  # seq is a Bio.SeqRecord object

            if sequence_string in existing_sequences:
                progress_dict[sequence_string] = existing_sequence_map[sequence_string]
                cached_count += 1
            else:
                progress_dict[sequence_string] = None
                novel_count += 1

        logger.info(
            f"Progress dict populated: {cached_count} cached, {novel_count} novel sequences"
        )
        logger.info(
            f"Total progress dict entries (excluding special keys): {len([k for k in progress_dict.keys() if not k.startswith('_')])}"
        )

        # Helper that atomically writes the progress file and persists the volume.
        def persist_progress():
            with open(VOLUME_DIR / in_progress_filename, "w") as fp:
                json.dump(progress_dict, fp)
            # Flush changes to the Modal volume so other processes/requests can see them.
            vol.commit()

        # Write the initial state so the client can poll immediately.
        persist_progress()

        # Nardini parameters.
        # NOTE: Should seed be randomized for each sequence calculation?
        random_seed = set_random_seed(DEFAULT_RANDOM_SEED)
        amino_acid_groupings = TYPEALL
        num_scrambled_sequences = NUM_SCRAMBLED_SEQUENCES

        # Run the analysis in parallel.
        calculation_start_time = time.time()

        # Directory where workers should place their individual zip files
        worker_output_dir = VOLUME_DIR / "zipfiles" / "by_idr"
        if not worker_output_dir.exists():
            worker_output_dir.mkdir(parents=True, exist_ok=True)

        # Parallelize the sequence processing
        with ProcessPoolExecutor(max_workers=WORKER_COUNT) as executor:
            # Create partial function with fixed parameters (except the sequence)
            process_func = functools.partial(
                process_single_sequence,
                amino_acid_groupings=amino_acid_groupings,
                num_scrambled_sequences=num_scrambled_sequences,
                random_seed=random_seed,
                output_dir=str(worker_output_dir),
            )
            novel_sequences = [
                seq
                for seq in parsed_sequences
                if str(seq.seq) not in existing_sequences
            ]
            logger.info(f"Processing {len(novel_sequences)} novel sequences.")

            # Submit all sequences and keep a mapping to retrieve results
            future_to_seq = {
                executor.submit(process_func, seq): seq for seq in novel_sequences
            }

            # Collect results as they complete
            from concurrent.futures import as_completed

            for future in as_completed(future_to_seq):
                seq = future_to_seq[future]
                sequence_string = str(seq.seq)
                try:
                    zip_path = future.result()
                    progress_dict[sequence_string] = zip_path
                    logger.info(
                        f"Processed sequence {seq.id} ({sequence_string[:50]}...) -> {zip_path}"
                    )
                except Exception as e:
                    logger.error(
                        f"Error processing sequence {seq.id} ({sequence_string[:50]}...): {e}"
                    )
                    progress_dict[sequence_string] = f"error: {e}"

                # Persist progress after each sequence completes so clients get near-real-time updates.
                persist_progress()
        calculation_end_time = time.time()
        logger.info(
            f"Nardini analysis took {calculation_end_time - calculation_start_time} seconds."
        )
        logger.info(f"Processed {len(novel_sequences)} sequences in parallel.")

        # Collect all zip files corresponding to the requested sequences
        # Filter out special keys (_status, merged_zip_filename) and only get actual sequence zip paths

        # Check each entry in progress_dict
        valid_entries = []
        invalid_entries = []
        for key, p in progress_dict.items():
            if key.startswith("_") or key == "merged_zip_filename":
                continue
            elif not p:
                logger.warning(f"Empty value for sequence: {key[:50]}...")
                invalid_entries.append((key, p))
                continue
            elif str(p).startswith("error"):
                logger.error(f"Error entry for sequence: {key[:50]}... -> {p}")
                invalid_entries.append((key, p))
                continue
            else:
                valid_entries.append((key, p))

        logger.info(
            f"Collecting zip files: {len(valid_entries)} valid, {len(invalid_entries)} invalid entries"
        )
        zip_files = [Path(p) for key, p in valid_entries]

        if not zip_files:
            raise RuntimeError("Nardini completed but no .zip archive was found.")

        final_output_dir = VOLUME_DIR / "zipfiles" / "by_fasta"
        if not final_output_dir.exists():
            final_output_dir.mkdir(parents=True, exist_ok=True)

        merged_zip_path = final_output_dir / f"{run_id}.zip"
        zip_file = mergeZips(zip_files, merged_zip_path)
        # logger.info(zip_files)
        logger.info(f"Created merged zip file: {zip_file}")

        # Update mappings for newly processed sequences before finalizing
        seq_mappings_to_update = {}
        for seq in parsed_sequences:
            sequence_string = str(seq.seq)
            if sequence_string not in existing_sequences and progress_dict.get(
                sequence_string
            ):
                # This is a newly processed sequence with a valid zip path
                seq_mappings_to_update[sequence_string] = {
                    "zip_filename": progress_dict[sequence_string],
                    "processed_at": time.time(),
                    "run_id": run_id,
                }

        # Write the updated mappings to the seq_zip_maps directory
        if seq_mappings_to_update:
            mapping_filename = f"{run_id}_mappings.json"
            mapping_file_path = seq_zip_maps_dir / mapping_filename
            with open(mapping_file_path, "w") as f:
                json.dump(seq_mappings_to_update, f, indent=2)
            logger.info(
                f"Updated mappings for {len(seq_mappings_to_update)} sequences in {mapping_filename}"
            )

        # Mark the run as completed and store final merged zip location
        progress_dict["_status"] = "completed"
        progress_dict["merged_zip_filename"] = str(zip_file)

        # Persist the final state (this includes a commit).
        persist_progress()

        end_time = time.time()
        return {
            "status": "completed",
            "zip_filename": zip_file,
            "calculation_time": calculation_end_time - calculation_start_time,
            "total_time": end_time - start_time,
        }
    except Exception as e:
        logger.error(f"Error in process_nardini_job: {e}")
        # Update progress to show error status
        progress_dict = {
            "_status": "failed",
            "error": str(e),
            "merged_zip_filename": None,
        }
        with open(VOLUME_DIR / f"{run_id}_in_progress.json", "w") as fp:
            json.dump(progress_dict, fp)
        vol.commit()
        raise


@app.function(
    image=web_image,
    volumes={str(VOLUME_DIR): vol.read_only()},
    # cpu=0.25,
    # memory=512,
)
@modal.asgi_app()
def fastapi_app():
    """Lightweight FastAPI application for handling uploads and job management."""
    api = FastAPI(title="Nardini Backend", version="1.0.0")

    @api.get("/health")
    async def health():
        """Health check endpoint."""
        return {"status": "healthy", "service": "nardini-backend"}

    @api.post("/upload_fasta", summary="Upload a FASTA file and spawn processing job")
    async def upload_fasta(file: UploadFile = File(...)) -> dict:
        """Validates the FASTA file, spawns processing job, and returns the run_id."""
        # Validate file type
        if not file.filename or not file.filename.lower().endswith(
            (".fasta", ".fa", ".fas")
        ):
            raise HTTPException(
                status_code=400,
                detail="Invalid file type. Please upload a FASTA file (.fasta, .fa, .fas).",
            )

        try:
            # Read and validate file content
            content = await file.read()
            if not content:
                raise HTTPException(status_code=400, detail="File is empty.")

            # MAX_FILE_SIZE is defined globally above
            if len(content) > MAX_FILE_SIZE:
                raise HTTPException(
                    status_code=413,
                    detail=f"File is too large. Limit is {MAX_FILE_SIZE / (1024 * 1024)}MB.",
                )

            # Pass raw content to worker for parsing (since nardini.utils is only in nardini_image)
            sequences_data = {
                "file_content": content.decode("utf-8"),  # Pass raw FASTA content
                "filename": file.filename,
            }

        except Exception as e:
            raise HTTPException(
                status_code=400, detail=f"Failed to parse FASTA file: {e}"
            )

        finally:
            pass  # No temp file cleanup needed since we pass content directly

        # Generate run ID and spawn processing job
        run_id = str(uuid4())

        # Spawn the heavy processing job
        job_call = process_nardini_job.spawn(sequences_data, run_id)

        return {
            "run_id": run_id,
            "job_id": job_call.object_id,
            "status": "submitted",
            "message": "Job submitted for processing",
        }

    @api.get("/status/{run_id}", summary="Check job status")
    async def check_status(run_id: str) -> dict:
        """Check the status of a processing job."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")

        try:
            # Read progress file from volume
            progress_file_path = f"{run_id}_in_progress.json"
            file_contents = vol.read_file(progress_file_path)

            if not file_contents:
                raise HTTPException(
                    status_code=404, detail="Data for this run not found"
                )

            # Read the file contents
            file_contents_bytes = b""
            for chunk in file_contents:
                file_contents_bytes += chunk

            progress_data = json.loads(file_contents_bytes)

            return {
                "run_id": run_id,
                "status": progress_data.get("_status", "unknown"),
                "progress": progress_data,
            }

        except Exception as e:
            logger.error(f"Error checking status for {run_id}: {e}")
            raise HTTPException(
                status_code=500, detail=f"Error checking status: {str(e)}"
            )

    @api.get("/download/{run_id}", summary="Download results")
    async def download_zip_file(run_id: str):
        """Download the results zip file."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")

        try:
            # Construct the path to the zip file in the volume
            volume_zip_path = f"zipfiles/by_fasta/{run_id}.zip"

            # Read the zip file from the volume
            file_contents = vol.read_file(volume_zip_path)

            if not file_contents:
                raise HTTPException(
                    status_code=404,
                    detail="Zip file not found. Job may still be in progress.",
                )

            # Collect all chunks into bytes
            zip_data = b""
            for chunk in file_contents:
                zip_data += chunk

            # Return the zip file as a downloadable response
            return Response(
                content=zip_data,
                media_type="application/zip",
                headers={"Content-Disposition": f"attachment; filename={run_id}.zip"},
            )

        except Exception as e:
            logger.error(f"Error downloading zip for {run_id}: {e}")
            raise HTTPException(
                status_code=500, detail=f"Error downloading file: {str(e)}"
            )

    return api
