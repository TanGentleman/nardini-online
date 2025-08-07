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
from pathlib import Path
import os
import time
import tempfile
import modal
import zipfile 
from uuid import uuid4
import logging
from typing import Any
import shutil
from typing_extensions import TypedDict

class SequenceInput(TypedDict):
    sequence: Any # Bio.SeqRecord object
    seq_uuid: str

class ProgressData(TypedDict):
    run_id: str
    status: str
    pending_sequences: list[str]

# ---------------------- Environment configuration ---------------------- #
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()
logging.basicConfig(
    level=LOG_LEVEL,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("nardini_backend")

# Configurable parameters with env-variable overrides
VOLUME_DIR = Path(os.getenv("VOLUME_DIR", "/data"))
VOLUME_NAME = os.getenv("VOLUME_NAME", "nardini_halophile_test_dev")
TIMEOUT_SECONDS = int(os.getenv("TIMEOUT_SECONDS", "10800"))  # default 3h
MAX_UPLOAD_MB = int(os.getenv("MAX_UPLOAD_MB", "10"))
MAX_FILE_SIZE = MAX_UPLOAD_MB * 1024 * 1024  # bytes
WORKER_COUNT = int(os.getenv("WORKER_COUNT", "16"))
# ---------------------------------------------------------------------- #

# Lightweight web image for FastAPI
web_image = (
    modal.Image.debian_slim()
    .pip_install(
        "fastapi[standard]==0.115.13",
        "python-multipart==0.0.20",
        "biopython==1.84",
    )
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
app = modal.App("nardini_halophile_test")
vol = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

# Web image imports (lightweight)
with web_image.imports():
    from fastapi import FastAPI, HTTPException, File, UploadFile
    from fastapi.responses import Response
    from Bio import SeqIO
    from io import StringIO

    # The following functions are copied from nardini.utils.py
    # NOTE: See https://github.com/mshinn23/nardini/blob/main/nardini/utils.py
    def read_sequences_from_filename(sequence_filename, default_name, verbose=False):
        """This is a helper function to read in sequences from a sequence FASTA file.
        This is a companion function to `read_sequences_from_string_list`.

        @param sequence_filename (str): The filepath of the file containing sequences.
                                        This file can contain FASTA records, or
                                        sequences that are saved on each newline.

        @param default_name (str):      This is the prefix name that's used when the
                                        file is found to only contain raw sequences
                                        on each new line.

        @param verbose (bool):          The extent to which information should be
                                        displayed throughout the analysis. Default:
                                        False.

        @returns seqio_sequences:       A list of sequence strings that were
                                        extracted from the sequence file.
        """
        seqio_sequences = list()
        if sequence_filename is None:
            raise RuntimeError('This parameter cannot be `None`. Exiting.')

        if type(sequence_filename) is not str:
            raise RuntimeError('This parameter can only be a string. Exiting.')

        if os.path.exists(sequence_filename):
            with open(sequence_filename, 'r') as seqfile:
                content = seqfile.read()
                if content.count('>') > 0:
                    seqio_sequences = list(SeqIO.parse(sequence_filename, 'fasta'))
                else:
                    # This means that there is a list of sequences that are separated
                    # by newlines. We split by any whitespace so as to capture
                    # carriage returns and line feeds.
                    raw_sequences = [s.strip() for s in content.split() if len(s.strip()) > 0]
                    seqio_sequences = read_sequences_from_string_list(raw_sequences, default_name)
        else:
            raise RuntimeError(f'Sequence filename: "{sequence_filename}" not found.')
        return seqio_sequences


    def read_sequences_from_string_list(list_of_sequences, default_name, verbose=False):
        """This is a helper function to read in sequences from a list of raw sequences.
        This is a companion function to `read_sequences_from_filename`.

        @param sequence_filename (str): The filepath of the file containing sequences.
                                        This file can contain FASTA records, or
                                        sequences that are saved on each newline.

        @param default_name (str):      This is the prefix name that's used when the
                                        file is found to only contain raw sequences
                                        on each new line.

        @param verbose (bool):          The extent to which information should be
                                        displayed throughout the analysis. Default:
                                        False.

        @returns seqio_sequences:       A list of sequence strings that were
                                        extracted from the sequence file.
        """
        sequences = list()
        # This means that we have to create a fake record using the sequence content.
        for index, sequence in enumerate(list_of_sequences, start=1):
            fasta = f'>{default_name}-{index}\n{sequence}'
            fake_record = SeqIO.read(StringIO(fasta), 'fasta')
            sequences.append(fake_record)

        if verbose:
            print('Number of sequences read: {num}'.format(num=len(sequences)), end='\n\n')
        return sequences

# Nardini image imports (heavy)
with nardini_image.imports():
    from nardini.constants import (
        NUM_SCRAMBLED_SEQUENCES,
        DEFAULT_RANDOM_SEED,
        TYPEALL,
    )
    from nardini.score_and_plot import calculate_zscore_and_plot
    from nardini.utils import set_random_seed

### NOTE Functions for stitching .zip results from each invocation together into final .zip
#helper fxn for mergeZips
def addtoSeqTSV(masterSeqTSV, tsvContent):
  import re
  lines = tsvContent.splitlines()
  for line in lines:
    line = line.strip()
    if not line:  # Skip empty lines
        continue
    # Split on 2+ consecutive spaces (the actual format used by Nardini)
    lineList = re.split(r' {2,}', line)
    # Skip header lines (check if first column is 'ID')
    if lineList and lineList[0] != 'ID':
      with open(masterSeqTSV, 'a') as f:
        # Join with tabs for proper TSV format in output
        f.write('\t'.join(lineList) + '\n')
  return


def _get_unique_filename(original_filename, used_filenames, zip_index):
  """Generate a unique filename by adding a suffix if the original is already used"""
  if original_filename not in used_filenames:
    return original_filename
  
  # Extract filename and extension
  if '.' in original_filename:
    name, ext = original_filename.rsplit('.', 1)
    ext = '.' + ext
  else:
    name = original_filename
    ext = ''
  
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
    destination_filepath = destination_filepath + '.zip'
  
  # Use a temporary file for the master TSV to avoid working directory issues
  with tempfile.NamedTemporaryFile(mode='w+', suffix='.tsv', delete=False) as temp_tsv:
    masterTSVPath = temp_tsv.name
    # Write header
    header = ["ID", "original_seq", "most_similar_seq", "sum_abs_zscore_original_seq", "sum_abs_zscore_scrambled_seq"]
    temp_tsv.write('\t'.join(header) + '\n')
  
  try:
    # Create a new zip file with the collected files
    with zipfile.ZipFile(destination_filepath, 'w') as merged_zip:
        used_filenames = set()  # Track filenames already added to prevent duplicates
        
        for zip_index, zipPath in enumerate(zipList):
            with zipfile.ZipFile(zipPath, 'r') as zip_ref:
                for file_info in zip_ref.infolist():
                    fileName = file_info.filename
                    if fileName != 'sequences.tsv':
                      # Handle potential duplicate filenames
                      unique_filename = _get_unique_filename(fileName, used_filenames, zip_index)
                      used_filenames.add(unique_filename)
                      
                      # Read the file content from the original zip and write to the new zip
                      with zip_ref.open(file_info) as source_file:
                        merged_zip.writestr(unique_filename, source_file.read())
                    else:
                      # Read the content of sequences.tsv from within the zip
                      with zip_ref.open(file_info) as seq_tsv_file:
                          seqTSVContent = seq_tsv_file.read().decode('utf-8') # Decode bytes to string
                          logger.info(f"Processing sequences.tsv from {zipPath}")
                          addtoSeqTSV(masterTSVPath, seqTSVContent) # Pass content to helper function

        # Write the content of the master_sequences.tsv into the zip while it's open
        with open(masterTSVPath, 'r') as master_tsv_file:
          master_tsv_content = master_tsv_file.read()
          merged_zip.writestr('master_sequences.tsv', master_tsv_content)
  
  finally:
    # Clean up the temporary file
    if os.path.exists(masterTSVPath):
      os.unlink(masterTSVPath)
  
  if not Path(destination_filepath).exists():
    raise ValueError("Zip creation failed!")
  logger.info(f"Merged zip file created at {destination_filepath}")
  return destination_filepath

def ensure_volume_directories():
  """Ensure the volume directories exist."""
  for dir in ["zipfiles/by_fasta", "zipfiles/by_idr", "runs"]:
    if not (VOLUME_DIR / dir).exists():
      logger.warning(f"Volume directory {dir} does not exist, creating it")
      (VOLUME_DIR / dir).mkdir(parents=True, exist_ok=True)


###
@app.function(
    image=nardini_image,
    volumes={str(VOLUME_DIR): vol},
    timeout=TIMEOUT_SECONDS,
    retries=2,
    scaledown_window=10,
    # enable_memory_snapshot=True
)
def process_single_sequence(
    sequence_input: SequenceInput,
) -> None:
    """Process a single sequence with Nardini analysis and return the created zip path.

    Each call is executed in its own temporary directory so that concurrently running
    processes cannot clash while writing output files. The produced zip archive is
    moved to the modal volume in zipfiles/by_idr/{seq_uuid}.zip.
    """
    sequence = sequence_input["sequence"]
    seq_uuid = sequence_input["seq_uuid"]
    amino_acid_groupings = TYPEALL
    num_scrambled_sequences = NUM_SCRAMBLED_SEQUENCES
    random_seed = set_random_seed(DEFAULT_RANDOM_SEED)
    output_dir = VOLUME_DIR / "zipfiles" / "by_idr"

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
                [sequence],  # Preserve original API (expects list)
                amino_acid_groupings,
                num_scrambled_sequences,
                random_seed,
            )

            # Locate the zip produced by Nardini (there should be exactly one)
            zip_files = [f for f in os.listdir(tmp_path) if f.startswith("nardini-data-") and f.endswith(".zip")]
            if len(zip_files) != 1:
                raise RuntimeError(
                    f"Expected one zip for sequence {sequence.id}, found {len(zip_files)}"
                )

            produced_zip_path = tmp_path / zip_files[0]

            # Build a deterministic destination name to avoid collisions
            # dest_zip_name = f"{seq_record.id}.zip"
            dest_path = Path(output_dir) / f"{seq_uuid}.zip"

            # Move the zip file to the shared output directory
            shutil.move(str(produced_zip_path), str(dest_path))
        
        except Exception as e:
            logger.error(f"Error processing sequence {sequence.id}: {e}")
            return None

        finally:
            # Restore the original cwd no matter what
            os.chdir(old_cwd)

    output_to_zip = dest_path
    if not output_to_zip.exists():
        raise RuntimeError(f"Output zip file {output_to_zip} does not exist")
    # Return the absolute path so the parent can build the mapping
    logging.info(f"Sequence {sequence.id} processed and added to volume at: {output_to_zip}")
    return None


@app.function(
    image=web_image,
    volumes={str(VOLUME_DIR): vol},
    # cpu=0.125,
    # memory=100,
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
        
        ensure_volume_directories()
        
        try:
            # Read and validate file content
            content = await file.read()
            if not content:
                raise HTTPException(status_code=400, detail="File is empty.")
            
            # MAX_FILE_SIZE is defined globally above
            if len(content) > MAX_FILE_SIZE:
                raise HTTPException(
                    status_code=413, 
                    detail=f"File is too large. Limit is {MAX_FILE_SIZE/(1024*1024)}MB."
                )

            # Parse sequences using the reference code with temporary file
            file_content = content.decode('utf-8')
            
            # Create temporary file to use with read_sequences_from_filename
            with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp_file:
                temp_file.write(file_content)
                temp_file_path = temp_file.name
            
            try:
                parsed_sequences = read_sequences_from_filename(
                    temp_file_path,
                    "sequence",  # default_name prefix  
                    verbose=False
                )
                
                if not parsed_sequences:
                    raise HTTPException(status_code=400, detail="No valid sequences found in FASTA file")
                    
            finally:
                # Clean up temporary file
                if os.path.exists(temp_file_path):
                    os.unlink(temp_file_path)
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"Failed to parse FASTA file: {e}")
        
        # Create set of current sequence strings
        current_sequence_strings = {str(seq.seq) for seq in parsed_sequences}
        logger.info(f"Processing {len(current_sequence_strings)} unique sequences")
        
        # Build sequence mapping and consolidate entries from existing .jsonl files
        runs_dir = VOLUME_DIR / "runs"
        
        # Consolidated entries for the new .jsonl file
        consolidated_entries = {}
        
        # Scan all existing .jsonl files for sequences in our current set
        for file in os.listdir(runs_dir):
            if file.endswith(".jsonl"):
                jsonl_path = runs_dir / file
                try:
                    with open(jsonl_path, "r") as f:
                        for line in f:
                            line = line.strip()
                            if line:  # Skip empty lines
                                entry = json.loads(line)
                                sequence_string = entry.get("sequence_string")
                                
                                # Only include sequences that are in our current upload
                                if sequence_string in current_sequence_strings:
                                    if sequence_string in consolidated_entries:
                                        if entry.get("status") == "cached":
                                            # Update newer entry for cached sequence
                                            consolidated_entries[sequence_string] = entry
                                    else:
                                        consolidated_entries[sequence_string] = entry

                except (json.JSONDecodeError, FileNotFoundError) as e:
                    logger.warning(f"Failed to read {jsonl_path}: {e}")
        
        # Generate run_id and prepare sequence inputs
        run_id = str(uuid4())
        sequence_inputs = []
        cached_sequences = 0
        novel_sequences = 0
        
        # Process sequences and create entries for any missing ones
        for seq in parsed_sequences:
            sequence_string = str(seq.seq)
            
            if sequence_string in consolidated_entries:
                entry = consolidated_entries[sequence_string]
                
                if entry["status"] == "cached":
                    cached_sequences += 1
                    logger.info(f"Found cached sequence: {seq.id}")
                else:
                    # Reuse existing UUID for pending sequence
                    # Consider it cached for now
                    seq_uuid = entry["seq_uuid"]
                    cached_sequences += 1
                    logger.info(f"Reusing existing UUID for pending sequence {seq.id}: {seq_uuid}")
            else:
                # Truly novel sequence, generate new UUID and entry
                seq_uuid = str(uuid4())
                novel_sequences += 1
                consolidated_entries[sequence_string] = {
                    "sequence_string": sequence_string,
                    "seq_uuid": seq_uuid,
                    "sequence_id": seq.id,
                    "status": "pending",
                    "zip_path": None
                }
                
                sequence_inputs.append(SequenceInput(
                    sequence=seq,
                    seq_uuid=seq_uuid
                ))
                logger.info(f"Novel sequence {seq.id} assigned new UUID: {seq_uuid}")
        
        # Create {run_id}.jsonl file with consolidated entries
        jsonl_path = runs_dir / f"{run_id}.jsonl"
        with open(jsonl_path, "w") as f:
           for entry in consolidated_entries.values():
               f.write(json.dumps(entry) + "\n")
            # for seq in parsed_sequences:
            #     sequence_string = str(seq.seq)
            #     entry = consolidated_entries[sequence_string]
            #     f.write(json.dumps(entry) + "\n")
        
        # Create clean .json metadata file (no full sequence data)
        run_metadata = {
            "status": "in_progress",
            "total_sequences": len(parsed_sequences),
            "cached_sequences": cached_sequences,
            "novel_sequences": novel_sequences,
            "merged_zip_filename": None,
            "started_at": time.time(),
            "completed_at": None
        }
        
        json_path = runs_dir / f"{run_id}.json"
        with open(json_path, "w") as f:
            json.dump(run_metadata, f, indent=2)
        
        # Spawn processing jobs for novel sequences using spawn_map
        if sequence_inputs:
            logger.info(f"Spawning {len(sequence_inputs)} processing jobs")
            await process_single_sequence.spawn_map.aio(sequence_inputs)
        else:
            logger.info("No novel sequences to process - all sequences are cached")
        
        logger.info(f"Created run {run_id} with {run_metadata['novel_sequences']} novel and {run_metadata['cached_sequences']} cached sequences")
        
        vol.commit()
        return {
            "run_id": run_id,
            "status": "submitted", 
            "message": "Job submitted for processing",
            "total_sequences": run_metadata["total_sequences"],
            "novel_sequences": run_metadata["novel_sequences"],
            "cached_sequences": run_metadata["cached_sequences"]
        }

    @api.get("/status/{run_id}", summary="Check job status")
    async def check_status(run_id: str) -> ProgressData:
        """Check the status of a processing job."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")
        
        try:
            # Read progress from volume
            # First read the json file
            runs_dir = VOLUME_DIR / "runs"
            json_path = runs_dir / f"{run_id}.json"
            with open(json_path, "r") as f:
                progress_data = json.load(f)
            
            pending_entries = []
            pending_sequences = []
            if progress_data.get("status") == "completed":
                return ProgressData(
                    run_id=run_id,
                    status="completed",
                    pending_sequences=[],
                )
            else:
                # Read the jsonl file
                jsonl_path = runs_dir / f"{run_id}.jsonl"
                with open(jsonl_path, "r") as f:
                    for line in f: 
                        entry = json.loads(line)
                        if entry.get("status") == "pending":
                            pending_entries.append(entry)
            
            if pending_entries:
                # Check if all of the pending sequences are in the by_idr directory
                idr_dir = VOLUME_DIR / "zipfiles" / "by_idr"
                idr_filenames = os.listdir(str(idr_dir))
                for entry in pending_entries:
                    if f"{entry.get('seq_uuid')}.zip" not in idr_filenames:
                        pending_sequences.append(entry.get("sequence_string"))
            
            if pending_sequences:
                return ProgressData(
                    run_id=run_id,
                    status="in_progress",
                    pending_sequences=pending_sequences
                )
            else:
                # Update the .json file to the completed status
                progress_data["status"] = "completed"
                progress_data["completed_at"] = time.time()
                with open(json_path, "w") as f:
                    json.dump(progress_data, f, indent=2)
                vol.commit()
                return ProgressData(
                    run_id=run_id,
                    status="completed",
                    pending_sequences=[]
                )
            
        except Exception as e:
            logger.error(f"Error checking status for {run_id}: {e}")
            raise HTTPException(status_code=500, detail=f"Error checking status: {str(e)}")

    @api.get("/download/{run_id}", summary="Download results")
    async def download_zip_file(run_id: str):
        """Download the results zip file."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")
        
        try:
            # Construct the path to the zip file in the volume
            merged_zip_paths = VOLUME_DIR / "zipfiles" / "by_fasta"
            # Read the zip file from the volume
            filenames = os.listdir(str(merged_zip_paths))
            for filename in filenames:
                if filename == f"{run_id}.zip":
                    file_contents = vol.read_file(f"zipfiles/by_fasta/{filename}")
                    break
            else:
                # Get all the seq_uuids or zip_paths from the runs directory
                runs_dir = VOLUME_DIR / "runs"
                pending_seq_uuids = []
                completed_zip_paths = []
                for file in os.listdir(runs_dir):
                    if file == f"{run_id}.jsonl":
                        jsonl_path = runs_dir / file
                        with open(jsonl_path, "r") as f:
                            for line in f:
                                entry = json.loads(line)
                                if entry.get("zip_path"):
                                    completed_zip_paths.append(entry.get("zip_path"))
                                else:
                                    pending_seq_uuids.append(entry.get("seq_uuid"))
                if pending_seq_uuids:
                    # Check if all of the pending seq_uuids are in the by_idr directory
                    idr_dir = VOLUME_DIR / "zipfiles" / "by_idr"
                    idr_filenames = os.listdir(str(idr_dir))
                    for seq_uuid in pending_seq_uuids:
                        if f"{seq_uuid}.zip" not in idr_filenames:
                            logger.warning(f"Zip file for seq_uuid {seq_uuid} not found. Download request will fail.")
                            raise HTTPException(status_code=404, detail="Zip file not found. Job may still be in progress.")
                        else:
                            completed_zip_paths.append(str(idr_dir / f"{seq_uuid}.zip"))
                
                # Merge the zip files
                merged_zip_path = VOLUME_DIR / "zipfiles" / "by_fasta" / f"{run_id}.zip"
                zip_file = mergeZips(completed_zip_paths, str(merged_zip_path))
                assert zip_file == str(merged_zip_path)
                assert zip_file.endswith(f"{run_id}.zip")
                vol.commit()
                file_contents = vol.read_file(f"zipfiles/by_fasta/{run_id}.zip")
            
            # Collect all chunks into bytes
            zip_data = b""
            for chunk in file_contents:
                zip_data += chunk
            
            # Return the zip file as a downloadable response
            return Response(
                content=zip_data,
                media_type="application/zip",
                headers={"Content-Disposition": f"attachment; filename={run_id}.zip"}
            )
            
        except Exception as e:
            logger.error(f"Error downloading zip for {run_id}: {e}")
            raise HTTPException(status_code=500, detail=f"Error downloading file: {str(e)}")

    return api