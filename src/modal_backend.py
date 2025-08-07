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

Volume layout and files written:
--------------------------------
- runs/{run_id}.json         -> metadata: status, counts, timestamps
- runs/{run_id}.jsonl        -> entries: one per sequence with fields: 
                                sequence_string, seq_uuid, sequence_id, status, zip_path
- zipfiles/by_idr/{seq_uuid}.zip  -> per-sequence Nardini output
- zipfiles/by_fasta/{run_id}.zip  -> merged archive of all per-sequence zips

Key environment variables (with defaults):
-----------------------------------------
- VOLUME_DIR=/data
- VOLUME_NAME=nardini_halophile_test_dev
- TIMEOUT_SECONDS=10800
- MAX_UPLOAD_MB=10
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
from typing import Any, Optional, List, Dict
import shutil
from typing_extensions import TypedDict

class SequenceInput(TypedDict):
    sequence: Any # Bio.SeqRecord object
    seq_uuid: str

class ProgressData(TypedDict):
    run_id: str
    status: str
    pending_sequences: list[str]

class SequenceEntry(TypedDict):
    sequence_string: str
    seq_uuid: str
    sequence_id: str
    status: str  # "pending" | "cached"
    zip_path: Optional[str]

# ---------------------- Environment configuration ---------------------- #
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()
logging.basicConfig(
    level=LOG_LEVEL,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("nardini_backend")

# Configurable parameters with env-variable overrides
VOLUME_DIR = Path(os.getenv("VOLUME_DIR", "/data"))
VOLUME_NAME = os.getenv("VOLUME_NAME", "nardini_volume_dev")
TIMEOUT_SECONDS = int(os.getenv("TIMEOUT_SECONDS", "10800"))  # default 3h
MAX_UPLOAD_MB = int(os.getenv("MAX_UPLOAD_MB", "10"))
MAX_FILE_SIZE = MAX_UPLOAD_MB * 1024 * 1024  # bytes
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
app = modal.App("nardini_backend_dev")
vol = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

# ---------------------- Helper utilities (paths, IO, scans) ---------------------- #
def runs_dir() -> Path:
    return VOLUME_DIR / "runs"


def zip_by_idr_dir() -> Path:
    return VOLUME_DIR / "zipfiles" / "by_idr"


def zip_by_fasta_dir() -> Path:
    return VOLUME_DIR / "zipfiles" / "by_fasta"


def run_json_path(run_id: str) -> Path:
    return runs_dir() / f"{run_id}.json"


def run_jsonl_path(run_id: str) -> Path:
    return runs_dir() / f"{run_id}.jsonl"


def read_json(path: Path) -> Dict[str, Any]:
    with open(path, "r") as f:
        return json.load(f)


def write_json(path: Path, data: Dict[str, Any]) -> None:
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def read_jsonl_entries(path: Path) -> List[Dict[str, Any]]:
    entries: List[Dict[str, Any]] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                entries.append(json.loads(line))
    return entries


def write_jsonl_entries(path: Path, entries: List[Dict[str, Any]]) -> None:
    with open(path, "w") as f:
        for entry in entries:
            f.write(json.dumps(entry) + "\n")


def idr_zip_filenames() -> set[str]:
    """Return the set of zip filenames present in `zipfiles/by_idr`.

    Uses pathlib for clarity and to avoid OS-specific listdir nuances.
    """
    return {path.name for path in zip_by_idr_dir().iterdir() if path.is_file()}


def missing_pending_sequences(pending_entries: List[SequenceEntry]) -> List[str]:
    existing = idr_zip_filenames()
    missing: List[str] = []
    for entry in pending_entries:
        if f"{entry['seq_uuid']}.zip" not in existing:
            missing.append(entry["sequence_string"])
    return missing


def zip_paths_for_run(run_id: str) -> List[str]:
    completed_zip_paths: List[str] = []
    try:
        entries = read_jsonl_entries(run_jsonl_path(run_id))
    except FileNotFoundError:
        entries = []
    pending_seq_uuids: List[str] = []
    for e in entries:
        if e.get("zip_path"):
            completed_zip_paths.append(e["zip_path"])  # type: ignore[arg-type]
        else:
            pending_seq_uuids.append(e["seq_uuid"])  # type: ignore[index]
    existing = idr_zip_filenames()
    for seq_uuid in pending_seq_uuids:
        expected = f"{seq_uuid}.zip"
        if expected not in existing:
            raise FileNotFoundError(f"Zip for {seq_uuid} not found")
        completed_zip_paths.append(str(zip_by_idr_dir() / expected))
    return completed_zip_paths


def consolidate_entries_for_sequences(current_sequence_strings: set[str]) -> Dict[str, SequenceEntry]:
    """Scan existing run JSONL files and return latest entries for the provided sequences.

    Prefers entries with status "cached" when duplicates exist across multiple runs.
    """
    consolidated: Dict[str, SequenceEntry] = {}
    for path in runs_dir().iterdir():
        if path.suffix == ".jsonl" and path.is_file():
            try:
                for entry in read_jsonl_entries(path):
                    seq_str = entry.get("sequence_string")
                    if seq_str in current_sequence_strings:
                        if seq_str in consolidated:
                            if entry.get("status") == "cached":
                                consolidated[seq_str] = entry  # type: ignore[assignment]
                        else:
                            consolidated[seq_str] = entry  # type: ignore[assignment]
            except (json.JSONDecodeError, FileNotFoundError) as e:
                logger.warning(f"Failed to read {path}: {e}")
    return consolidated


def seqrecord_from_entry(entry: SequenceEntry):
    # Local import to avoid cross-image import issues
    from Bio.SeqRecord import SeqRecord  # type: ignore
    from Bio.Seq import Seq  # type: ignore
    return SeqRecord(Seq(entry["sequence_string"]), id=entry["sequence_id"], description="")

# ------------------------------------------------------------------------------- #

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
def add_to_master_sequences_tsv(master_tsv_path: str, tsv_content: str) -> None:
    """Append rows from a per-zip sequences.tsv into the master TSV, skipping headers."""
    import re
    lines = tsv_content.splitlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        columns = re.split(r" {2,}", line)
        if columns and columns[0] != "ID":
            with open(master_tsv_path, "a") as f:
                f.write("\t".join(columns) + "\n")


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


def merge_zip_archives(zip_list: List[str], destination_filepath: str) -> str:
    """Merge multiple zip files into a single zip.

    - Copies all files from each input zip into the merged archive, resolving
      filename collisions by appending a suffix.
    - Concatenates all per-zip `sequences.tsv` into a `master_sequences.tsv` in
      the merged archive.
    """
    if not zip_list:
        raise ValueError("No zip files provided to merge")

    destination_filepath = str(destination_filepath)
    if not destination_filepath.endswith(".zip"):
        destination_filepath = destination_filepath + ".zip"

    # Use a temporary file for the master TSV to avoid working directory issues
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".tsv", delete=False) as temp_tsv:
        master_tsv_path = temp_tsv.name
        header = [
            "ID",
            "original_seq",
            "most_similar_seq",
            "sum_abs_zscore_original_seq",
            "sum_abs_zscore_scrambled_seq",
        ]
        temp_tsv.write("\t".join(header) + "\n")

    try:
        with zipfile.ZipFile(destination_filepath, "w") as merged_zip:
            used_filenames = set()

            for zip_index, zip_path in enumerate(zip_list):
                with zipfile.ZipFile(zip_path, "r") as zip_ref:
                    for file_info in zip_ref.infolist():
                        file_name = file_info.filename
                        if file_name != "sequences.tsv":
                            unique_filename = _get_unique_filename(file_name, used_filenames, zip_index)
                            used_filenames.add(unique_filename)
                            with zip_ref.open(file_info) as source_file:
                                merged_zip.writestr(unique_filename, source_file.read())
                        else:
                            with zip_ref.open(file_info) as seq_tsv_file:
                                seq_tsv_content = seq_tsv_file.read().decode("utf-8")
                                logger.info(f"Processing sequences.tsv from {zip_path}")
                                add_to_master_sequences_tsv(master_tsv_path, seq_tsv_content)

            with open(master_tsv_path, "r") as master_tsv_file:
                master_tsv_content = master_tsv_file.read()
                merged_zip.writestr("master_sequences.tsv", master_tsv_content)

    finally:
        # Clean up the temporary file
        path_obj = Path(master_tsv_path)
        if path_obj.exists():
            path_obj.unlink()

    if not Path(destination_filepath).exists():
        raise ValueError("Zip creation failed!")
    logger.info(f"Merged zip file created at {destination_filepath}")
    return destination_filepath

def ensure_volume_directories() -> None:
    """Ensure the required directories exist in the mounted volume."""
    for directory in [zip_by_fasta_dir(), zip_by_idr_dir(), runs_dir()]:
        if not directory.exists():
            logger.warning(f"Volume directory {directory} does not exist, creating it")
            directory.mkdir(parents=True, exist_ok=True)

@app.function(
    image=nardini_image,
    volumes={str(VOLUME_DIR): vol},
    timeout=TIMEOUT_SECONDS,
    scaledown_window=5,
    # enable_memory_snapshot=True
)  
def retry_pending_sequences(run_id: str):
    """Retry processing for sequences that are still pending."""
    json_path = run_json_path(run_id)
    progress_data = read_json(json_path)

    if progress_data.get("status") == "completed":
        logger.warning(f"Run {run_id} is already completed, skipping retry")
        return

    pending_entries: List[SequenceEntry] = []
    for entry in read_jsonl_entries(run_jsonl_path(run_id)):
        if entry.get("status") == "pending":
            pending_entries.append(entry)  # type: ignore[arg-type]

    idr_filenames = idr_zip_filenames()
    pending_sequences: List[SequenceInput] = []
    for entry in pending_entries:
        if f"{entry['seq_uuid']}.zip" not in idr_filenames:
            seq_record = seqrecord_from_entry(entry)
            pending_sequences.append(SequenceInput(sequence=seq_record, seq_uuid=entry["seq_uuid"]))

    # Spawn processing jobs for novel sequences using spawn_map
    if pending_sequences:
        logger.info(f"Spawning {len(pending_sequences)} processing jobs")
        process_single_sequence.spawn_map(pending_sequences)
    else:
        logger.info("No novel sequences to process - all sequences are cached")

	
    

###
#processing_cpu_request_and_limit = (0.125, 0.15)
@app.function(
    image=nardini_image,
    volumes={str(VOLUME_DIR): vol},
    timeout=TIMEOUT_SECONDS,
    scaledown_window=10,
    #cpu=processing_cpu_request_and_limit,
    #max_containers=10
    # enable_memory_snapshot=True
)
def process_single_sequence(
    sequence_input: SequenceInput,
) -> None:
    """Process a single sequence with Nardini analysis and return the created zip path.

    Each call is executed in its own temporary directory so that concurrently running
    processes cannot clash while writing output files. The produced zip archive is
    moved to the modal volume in `zipfiles/by_idr/{seq_uuid}.zip`.
    """
    sequence = sequence_input["sequence"]
    seq_uuid = sequence_input["seq_uuid"]
    amino_acid_groupings = TYPEALL
    num_scrambled_sequences = NUM_SCRAMBLED_SEQUENCES
    random_seed = set_random_seed(DEFAULT_RANDOM_SEED)
    output_dir = zip_by_idr_dir()

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
            zip_files = list(Path(tmp_path).glob("nardini-data-*.zip"))
            if len(zip_files) != 1:
                raise RuntimeError(
                    f"Expected one zip for sequence {sequence.id}, found {len(zip_files)}"
                )

            produced_zip_path = zip_files[0]

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
        # Consolidate entries from existing runs that match current sequences
        consolidated_entries = consolidate_entries_for_sequences(current_sequence_strings)
        
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
        jsonl_path = run_jsonl_path(run_id)
        write_jsonl_entries(jsonl_path, list(consolidated_entries.values()))
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
        
        json_path = run_json_path(run_id)
        write_json(json_path, run_metadata)
        
        # Spawn processing jobs for novel sequences using spawn_map
        job_ids = []
        if sequence_inputs:
            logger.info(f"Spawning {len(sequence_inputs)} processing jobs")
            jobs = await process_single_sequence.spawn_map.aio(sequence_inputs)
            job_ids = [job.object_id for job in jobs]
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
            "cached_sequences": run_metadata["cached_sequences"],
            "job_ids": job_ids
        }

    @api.get("/status/{run_id}", summary="Check job status")
    async def check_status(run_id: str) -> ProgressData:
        """Check the status of a processing job."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")
        
        try:
            # Read progress from volume
            json_path = run_json_path(run_id)
            if not json_path.exists():
                raise HTTPException(status_code=404, detail="Run not found")
            progress_data = read_json(json_path)
            
            pending_entries: List[SequenceEntry] = []
            pending_sequences: List[str] = []
            if progress_data.get("status") == "completed":
                return ProgressData(
                    run_id=run_id,
                    status="completed",
                    pending_sequences=[],
                )
            else:
                # Read the jsonl file and gather pending entries
                for entry in read_jsonl_entries(run_jsonl_path(run_id)):
                    if entry.get("status") == "pending":
                        pending_entries.append(entry)  # type: ignore[arg-type]
            
            if pending_entries:
                pending_sequences = missing_pending_sequences(pending_entries)
            
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
                write_json(json_path, progress_data)
                vol.commit()
                return ProgressData(
                    run_id=run_id,
                    status="completed",
                    pending_sequences=[]
                )
            
        except HTTPException:
            raise
        except Exception as e:
            logger.error(f"Error checking status for {run_id}: {e}")
            raise HTTPException(status_code=500, detail=f"Error checking status: {str(e)}")

    @api.get("/download/{run_id}", summary="Download results")
    async def download_zip_file(run_id: str):
        """Download the results zip file."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")
        
        try:
            merged_zip_path = zip_by_fasta_dir() / f"{run_id}.zip"
            if not merged_zip_path.exists():
                try:
                    completed_zip_paths = zip_paths_for_run(run_id)
                except FileNotFoundError as e:
                    logger.warning(str(e))
                    raise HTTPException(status_code=404, detail="Zip file not found. Job may still be in progress.")

                zip_file = merge_zip_archives(completed_zip_paths, str(merged_zip_path))
                assert zip_file == str(merged_zip_path)
                assert zip_file.endswith(f"{run_id}.zip")
                vol.commit()

            # Read the zip file directly from the filesystem using pathlib
            zip_data = merged_zip_path.read_bytes()
            
            # Return the zip file as a downloadable response
            return Response(
                content=zip_data,
                media_type="application/zip",
                headers={"Content-Disposition": f"attachment; filename={run_id}.zip"}
            )
            
        except Exception as e:
            logger.error(f"Error downloading zip for {run_id}: {e}")
            raise HTTPException(status_code=500, detail=f"Error downloading file: {str(e)}")

    @api.get("/retry/{run_id}", summary="Retry processing for sequences that are still pending")
    async def retry_sequences(run_id: str):
        """Retry processing for sequences that are still pending."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")
        retry_pending_sequences.spawn(run_id)

    return api