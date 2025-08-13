"""
This script defines a Modal application that serves a FastAPI backend for running
Nardini, a tool for analyzing Intrinsically Disordered Regions (IDRs) in proteins.

The backend uses Modal's job queue pattern with:
- A lightweight FastAPI web server for handling uploads and spawning jobs
- A heavy Nardini worker for processing sequences

Endpoints:
- `/health` — service health check
- `/upload_fasta` — accepts a FASTA file, spawns jobs, and returns a run_id
- `/status/{run_id}` — returns pending sequence IDs or completion state
- `/download/{run_id}` — returns merged results zip when complete
- `/retry/{run_id}` — re-checks cache and re-spawns work

Deployment:
------------
modal deploy modal_backend.py

Example usage with curl:
-----------------------
# Replace with your deployed Modal app URL.
curl -X POST "YOUR_MODAL_APP_URL/upload_fasta" \
-F "file=@data/fasta_inputs/Halophile-pHtolerant-yeast-first16.fasta"

Volume layout and files written:
--------------------------------
- runs/{run_id}.json -> metadata for the run including:
    - status, timestamps (submitted_at, completed_at), other metadata
    - sequences: mapping of sequence_string -> {
        sequence_id, seq_uuid, status, start_time, end_time, zip_path, job_id
      }
- zipfiles/by_idr/{seq_uuid}.zip  -> per-sequence Nardini output
- zipfiles/by_fasta/{run_id}.zip  -> merged archive of all per-sequence zips

Key environment variables (with defaults):
-----------------------------------------
- APP_NAME=nardini_backend_dev
- VOLUME_DIR=/data
- VOLUME_NAME=nardini_volume_dev
- TIMEOUT_SECONDS=10800
- MAX_UPLOAD_MB=10
"""

import json
import logging
import os
import shutil
import tempfile
import time
import zipfile
from pathlib import Path
from typing import Any, Dict, List
from uuid import uuid4

import modal
from typing_extensions import TypedDict


class SequenceInput(TypedDict):
    sequence: Any  # Bio.SeqRecord object
    seq_uuid: str


# ---------------------- Environment configuration ---------------------- #
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()
logging.basicConfig(
    level=LOG_LEVEL,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

# Configurable parameters with env-variable overrides
APP_NAME = os.getenv("APP_NAME", "nardini_backend_dev")
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
    .add_local_python_source("schemas")
)

# Heavy processing image for Nardini
nardini_image = (
    modal.Image.debian_slim(python_version="3.9")
    .uv_pip_install(
        "pydantic==2.10.6",
        "numpy==1.19.5",
    )
    .uv_pip_install("nardini==1.1.1")
).add_local_python_source("schemas")

# Create a Modal application with shared volume
app = modal.App(APP_NAME)
vol = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

with web_image.imports():
    from schemas import RunData, SequenceData, SequencesMapping


# ---------------------- Helper utilities (paths, IO, scans) ---------------------- #
def get_runs_dir() -> Path:
    """Absolute path to the directory storing per-run JSON metadata files."""
    return VOLUME_DIR / "runs"


def get_zip_by_idr_dir() -> Path:
    """Absolute path to per-sequence zip outputs, keyed by `seq_uuid`."""
    return VOLUME_DIR / "zipfiles" / "by_idr"


def get_zip_by_fasta_dir() -> Path:
    """Absolute path to merged zip outputs, keyed by `run_id`."""
    return VOLUME_DIR / "zipfiles" / "by_fasta"


def get_run_json_path(run_id: str) -> Path:
    return get_runs_dir() / f"{run_id}.json"


def read_json(path: Path) -> dict:
    with open(path, "r") as f:
        return json.load(f)


def write_json(path: Path, data: dict) -> None:
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def get_run_metadata(run_id: str) -> RunData:
    path = get_run_json_path(run_id)
    return read_json(path)


def write_run_metadata_to_volume(run_id: str, data: RunData) -> None:
    path = get_run_json_path(run_id)
    write_json(path, data)


def create_sequences_data(sequence_strings: set[str]) -> SequencesMapping:
    """Build a mapping of sequence string -> metadata by scanning prior runs.

    If a sequence appears in multiple runs, a cached status is preferred over pending.
    Any sequences not found are initialized as pending with a fresh UUID.
    """
    # Copy to avoid mutating the caller's set
    sequence_strings = set(sequence_strings)
    sequences_data: SequencesMapping = {}
    for path in get_runs_dir().iterdir():
        if path.suffix == ".json" and path.is_file():
            try:
                data = get_run_metadata(path.stem)
                sequences = data["sequences"]
                if not sequences:
                    continue
                for seq_str, seq_data in sequences.items():
                    if seq_str not in sequence_strings:
                        continue
                    sequence_strings.remove(seq_str)

                    seq_status = seq_data["status"]
                    if seq_status == "pending" or seq_status == "pending_external":
                        seq_data["status"] = "pending_external"
                    elif seq_status == "cached" or seq_status == "complete":
                        seq_data["status"] = "cached"
                    else:
                        raise ValueError(f"Invalid sequence status: {seq_status}")

                    if seq_str in sequences_data:
                        # Prefer cached over pending
                        if sequences_data[seq_str]["status"] != "cached":
                            sequences_data[seq_str] = seq_data
                    else:
                        sequences_data[seq_str] = seq_data
            except (json.JSONDecodeError, FileNotFoundError) as e:
                logger.warning(f"Failed to read {path}: {e}")
            except Exception as e:
                logger.warning(f"Unexpected error reading {path}: {e}")

    for seq_str in sequence_strings:
        sequences_data[seq_str] = SequenceData(
            sequence_id="TEMP_ID",
            status="pending",
            start_time=None,
            end_time=None,
            seq_uuid=str(uuid4()),
            zip_path=None,
            job_id=None,
        )
    return sequences_data


def create_seqrecord(sequence_string: str, sequence_id: str):
    """Create a BioPython SeqRecord from a raw sequence string and identifier."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    return SeqRecord(Seq(sequence_string), id=sequence_id)


# ------------------ New helper functions for sequence status handling ------------------ #


def get_pending_uuids(sequence_dict: SequencesMapping) -> List[str]:
    """Return list of `seq_uuid` values for sequences still pending processing."""
    pending_seq_uuids = []
    for sequence_data in sequence_dict.values():
        if sequence_data["status"] in ("pending", "pending_external"):
            assert sequence_data["seq_uuid"], f"Error in Cache! {sequence_data}"
            pending_seq_uuids.append(sequence_data["seq_uuid"])
    return pending_seq_uuids


def find_completed_in_cache(pending_seq_uuids: List[str]) -> Dict[str, float]:
    """Locate expected zipfiles in the cache directory and return mapping of filename -> ctime."""
    expected_zipfiles = [f"{seq_uuid}.zip" for seq_uuid in pending_seq_uuids]
    idr_zip_dir = get_zip_by_idr_dir()
    completed: Dict[str, float] = {}
    for file in idr_zip_dir.iterdir():
        if file.name in expected_zipfiles:
            completed[file.name] = file.stat().st_ctime
    return completed


def update_sequences_with_completed(
    sequences_data: SequencesMapping,
    completed_mapping: Dict[str, float],
) -> List[str]:
    """Mutate `sequences_data` in place with completion times.

    Returns a list of `sequence_id` values that are still pending.
    """
    pending_sequence_ids: List[str] = []
    idr_zip_dir = get_zip_by_idr_dir()
    for seq_data in sequences_data.values():
        seq_uuid = seq_data["seq_uuid"]
        expected = f"{seq_uuid}.zip" if seq_uuid else None
        if expected and expected in completed_mapping:
            seq_data["zip_path"] = str(idr_zip_dir / expected)
            seq_data["seq_uuid"] = None
            seq_data["end_time"] = completed_mapping[expected]
            seq_data["status"] = "complete"
        elif seq_data["status"] in ("pending", "pending_external"):
            pending_sequence_ids.append(seq_data["sequence_id"])
    return pending_sequence_ids


def get_completed_zip_paths(
    sequences_data: SequencesMapping, require_all_complete: bool = False
) -> List[str]:
    """Return absolute paths for zipfiles that are available for the provided sequences."""
    completed_zip_paths: List[str] = []
    idr_zip_dir = get_zip_by_idr_dir()
    for seq_data in sequences_data.values():
        if seq_data["zip_path"]:
            completed_zip_paths.append(seq_data["zip_path"])
            continue
        seq_uuid = seq_data["seq_uuid"]
        if seq_uuid:
            candidate = idr_zip_dir / f"{seq_uuid}.zip"
            if candidate.exists():
                completed_zip_paths.append(str(candidate))
        elif require_all_complete:
            raise ValueError(f"Sequence {seq_data['sequence_id']} is not complete")
    return completed_zip_paths


# ------------------------------------------------------------------------------- #


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
    from Bio import SeqIO

    seqio_sequences = list()
    if sequence_filename is None:
        raise RuntimeError("This parameter cannot be `None`. Exiting.")

    if type(sequence_filename) is not str:
        raise RuntimeError("This parameter can only be a string. Exiting.")

    if os.path.exists(sequence_filename):
        with open(sequence_filename, "r") as seqfile:
            content = seqfile.read()
            if content.count(">") > 0:
                seqio_sequences = list(SeqIO.parse(sequence_filename, "fasta"))
            else:
                # This means that there is a list of sequences that are separated
                # by newlines. We split by any whitespace so as to capture
                # carriage returns and line feeds.
                raw_sequences = [
                    s.strip() for s in content.split() if len(s.strip()) > 0
                ]
                seqio_sequences = read_sequences_from_string_list(
                    raw_sequences, default_name
                )
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
    from io import StringIO

    from Bio import SeqIO

    sequences = list()
    # This means that we have to create a fake record using the sequence content.
    for index, sequence in enumerate(list_of_sequences, start=1):
        fasta = f">{default_name}-{index}\n{sequence}"
        fake_record = SeqIO.read(StringIO(fasta), "fasta")
        sequences.append(fake_record)

    if verbose:
        print("Number of sequences read: {num}".format(num=len(sequences)), end="\n\n")
    return sequences


### NOTE Functions for stitching .zip results from each invocation together into final .zip
# helper fxn for mergeZips
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
    with tempfile.NamedTemporaryFile(
        mode="w+", suffix=".tsv", delete=False
    ) as temp_tsv:
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
                            unique_filename = _get_unique_filename(
                                file_name, used_filenames, zip_index
                            )
                            used_filenames.add(unique_filename)
                            with zip_ref.open(file_info) as source_file:
                                merged_zip.writestr(unique_filename, source_file.read())
                        else:
                            with zip_ref.open(file_info) as seq_tsv_file:
                                seq_tsv_content = seq_tsv_file.read().decode("utf-8")
                                logger.info(f"Processing sequences.tsv from {zip_path}")
                                add_to_master_sequences_tsv(
                                    master_tsv_path, seq_tsv_content
                                )

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


def ensure_volume_directories() -> bool:
    """Ensure the required directories exist in the mounted volume."""
    was_created = False
    for directory in [get_zip_by_fasta_dir(), get_zip_by_idr_dir(), get_runs_dir()]:
        if not directory.exists():
            logger.warning(f"Volume directory {directory} does not exist, creating it")
            directory.mkdir(parents=True, exist_ok=True)
            was_created = True
    return was_created


def sanitize_output_filename(filename: str) -> str:
    if not filename:
        raise ValueError("Filename is required")
    if len(filename) > 250:
        filename = filename[:250]
    # Remove any potential .fasta extension
    suffix = filename.split(".")[-1]
    if suffix in ["fasta", "fa", "fas"]:
        suffix_len = len(suffix) + 1
        filename = filename[:-suffix_len]
    if not filename.endswith(".zip"):
        filename = filename + ".zip"
    return filename


@app.function(
    image=nardini_image,
    volumes={str(VOLUME_DIR): vol},
    timeout=TIMEOUT_SECONDS,
    scaledown_window=5,
)
def retry_pending_sequences(run_id: str):
    """Retry processing for sequences that are still pending."""
    progress_data = get_run_metadata(run_id)

    if progress_data.get("status") == "complete":
        logger.warning(f"Run {run_id} is already complete, skipping retry")
        return

    # NOTE: This function will be refactored if errors persist
    sequences_data = progress_data["sequences"]
    pending_uuids = get_pending_uuids(sequences_data)
    if not pending_uuids:
        logger.error(
            "No pending sequences...run status should have been marked complete!"
        )
        return
    completed_mapping = find_completed_in_cache(pending_uuids)
    pending_sequence_ids = update_sequences_with_completed(
        sequences_data, completed_mapping
    )
    logger.info(f"Retrying {len(pending_sequence_ids)} pending sequences")
    logger.info(f"Pending ids: {pending_sequence_ids}")
    job_ids = []
    for seq_str, data in sequences_data.items():
        if data["status"] == "pending":
            seq_uuid = data["seq_uuid"]
            if not seq_uuid:
                logger.error(f"Sequence {seq_str} has no seq_uuid")
                continue
            seq_record = create_seqrecord(seq_str, data.get("sequence_id"))
            call = process_single_sequence.spawn(
                SequenceInput(sequence=seq_record, seq_uuid=seq_uuid)
            )
            job_ids.append(call.object_id)
            # Update sequences map with job details for each sequence in the batch
            sequences_data[seq_str] = SequenceData(
                start_time=time.time(),
                end_time=None,
                seq_uuid=seq_uuid,
                sequence_id=seq_record.id,
                status="pending",
                zip_path=None,
                job_id=call.object_id,
            )
    logger.info(f"Spawned {len(job_ids)} jobs. Ids: {job_ids}")
    return


@app.function(
    image=nardini_image,
    volumes={str(VOLUME_DIR): vol},
)
def migrate_run_metadata() -> None:
    """Modify the run metadata all JSON files."""
    for json_file in get_runs_dir().glob("*.json"):
        run_metadata = get_run_metadata(json_file.stem)
        run_metadata["output_filename"] = sanitize_output_filename(
            run_metadata["output_filename"]
        )
        logger.info(
            f"Modifying run metadata for {json_file.stem} to {run_metadata['output_filename']}"
        )
        write_run_metadata_to_volume(json_file.stem, run_metadata)
    vol.commit()
    return


###
@app.function(
    image=nardini_image,
    volumes={str(VOLUME_DIR): vol},
    timeout=TIMEOUT_SECONDS,
    scaledown_window=10,
    # max_containers=50,
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
    try:
        output_path = _process_one_sequence_to_volume(sequence_input)
        sequence = sequence_input["sequence"]
        logging.info(
            f"Sequence {sequence.id} processed and added to volume at: {output_path}"
        )
    except Exception as e:
        sequence = sequence_input.get("sequence")
        seq_id = getattr(sequence, "id", "<unknown>")
        logger.error(f"Error processing sequence {seq_id}: {e}")
    return None


def _process_one_sequence_to_volume(sequence_input: SequenceInput) -> str:
    """Shared helper to run Nardini for one sequence and write the zip into the volume.

    Returns absolute path to the written zip in `zipfiles/by_idr`.
    Raises on error.
    """
    from nardini.constants import (
        DEFAULT_RANDOM_SEED,
        NUM_SCRAMBLED_SEQUENCES,
        TYPEALL,
    )
    from nardini.score_and_plot import calculate_zscore_and_plot
    from nardini.utils import set_random_seed

    sequence = sequence_input["sequence"]
    seq_uuid = sequence_input["seq_uuid"]

    amino_acid_groupings = TYPEALL
    num_scrambled_sequences = NUM_SCRAMBLED_SEQUENCES
    random_seed = set_random_seed(DEFAULT_RANDOM_SEED)
    output_dir = get_zip_by_idr_dir()

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        old_cwd = os.getcwd()
        os.chdir(tmp_path)
        try:
            calculate_zscore_and_plot(
                [sequence],
                amino_acid_groupings,
                num_scrambled_sequences,
                random_seed,
            )

            files = list(tmp_path.iterdir())
            zip_files = [f for f in files if f.name.endswith(".zip")]
            if len(zip_files) != 1:
                raise RuntimeError(
                    f"Expected one zip for sequence {sequence.id}, found {len(zip_files)}"
                )

            produced_zip_path = zip_files[0]
            dest_path = Path(output_dir) / f"{seq_uuid}.zip"
            shutil.move(str(produced_zip_path), str(dest_path))
        finally:
            os.chdir(old_cwd)

    if not Path(dest_path).exists():
        raise RuntimeError(f"Output zip file {dest_path} does not exist")
    return str(dest_path)


@app.function(
    image=nardini_image,
    volumes={str(VOLUME_DIR): vol},
    timeout=TIMEOUT_SECONDS,
    scaledown_window=10,
)
def process_16_sequences(sequence_inputs: List[SequenceInput]) -> None:
    """Process up to 16 sequences in parallel within a single worker.

    Each sequence is executed in its own temporary directory to avoid file clashes.
    The produced zip archives are moved to the modal volume in `zipfiles/by_idr/{seq_uuid}.zip`.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    if not sequence_inputs:
        return None

    # Use fixed 16 workers, capped by input size
    max_local_workers = min(16, len(sequence_inputs))

    with ProcessPoolExecutor(max_workers=max_local_workers) as executor:
        future_to_input = {
            executor.submit(_process_one_sequence_to_volume, si): si
            for si in sequence_inputs
        }
        for future in as_completed(future_to_input):
            si = future_to_input[future]
            sequence = si["sequence"]
            try:
                _ = future.result()
                logger.info(f"Processed sequence {sequence.id}")
            except Exception as e:
                logger.error(f"Error processing sequence {sequence.id}: {e}")
                continue
    return None


@app.function(
    image=web_image,
    volumes={str(VOLUME_DIR): vol},
    # cpu=0.125,
    # memory=100,
)
@modal.asgi_app(requires_proxy_auth=True)
def fastapi_app():
    from fastapi import FastAPI, File, HTTPException, UploadFile
    from fastapi.responses import Response

    from schemas import (
        HealthResponse,
        RetryResponse,
        StatusResponse,
        UploadFastaResponse,
    )

    """Lightweight FastAPI application for handling uploads and job management."""
    api = FastAPI(title="Nardini Backend", version="1.0.0")

    @api.get("/health")
    async def health():
        """Health check endpoint."""
        return HealthResponse(status="healthy")

    # TODO: Add more HTTPException handling and logging for errors
    @api.post("/upload_fasta", summary="Upload a FASTA file and spawn processing job")
    async def upload_fasta(
        file: UploadFile = File(...), output_filename: str = None
    ) -> dict:
        """Validates the FASTA file, spawns processing job, and returns the run_id."""
        # Validate file type
        logger.debug(f"Uploading file: {file.filename} at {time.time()}")
        original_filename = file.filename
        if not original_filename or not original_filename.lower().endswith((
            ".fasta",
            ".fa",
            ".fas",
        )):
            raise HTTPException(
                status_code=400,
                detail="Invalid file type. Please upload a FASTA file (.fasta, .fa, .fas).",
            )
        if output_filename:
            zip_output_filename = sanitize_output_filename(output_filename)
        else:
            zip_output_filename = sanitize_output_filename(original_filename)

        was_created = ensure_volume_directories()
        if was_created:
            vol.commit()

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

            # Parse sequences using the reference code with a temporary file
            file_content = content.decode("utf-8")

            # Create temporary file to use with read_sequences_from_filename
            with tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".fasta"
            ) as temp_file:
                temp_file.write(file_content)
                temp_file_path = temp_file.name

            logger.debug(f"Starting to parse sequences at {time.time()}")
            try:
                parsed_sequences = read_sequences_from_filename(
                    temp_file_path,
                    "sequence",  # default_name prefix
                    verbose=False,
                )

                if not parsed_sequences:
                    raise HTTPException(
                        status_code=400, detail="No valid sequences found in FASTA file"
                    )

            finally:
                # Clean up temporary file
                if os.path.exists(temp_file_path):
                    os.unlink(temp_file_path)
        except Exception as e:
            raise HTTPException(
                status_code=400, detail=f"Failed to parse FASTA file: {e}"
            )

        # Create set of current sequence strings
        sequence_strings = {str(seq.seq) for seq in parsed_sequences}
        logger.info(f"Processing {len(sequence_strings)} unique sequences")

        # Build sequence map from existing runs that match current sequences
        sequences_data = create_sequences_data(sequence_strings)

        # Generate run_id and prepare sequence inputs
        run_id = str(uuid4())
        novel_sequences = []

        for sequence in parsed_sequences:
            seq_str = str(sequence.seq)
            seq_data = sequences_data[seq_str]
            if seq_data.get("status") == "pending":
                novel_sequences.append(sequence)
            else:
                if seq_data.get("status") == "pending_external":
                    logger.warning(
                        f"We must ensure external sequence uuid {seq_data.get('seq_uuid')} is not stale/failed!"
                    )
                    continue
                assert seq_data.get("status") == "cached"

        novel_sequences.sort(key=lambda x: len(str(x.seq)), reverse=True)

        job_ids = []
        if novel_sequences:
            logger.debug(f"First sequence length: {len(str(novel_sequences[0].seq))}")
            logger.debug(f"Last sequence length: {len(str(novel_sequences[-1].seq))}")
            logger.info(f"Spawning {len(novel_sequences)} sequences in batches of 16")
            logger.debug(f"Starting to spawn at {time.time()}")

            USE_BATCH_PROCESSING = True
            BATCH_SIZE = 16
            if USE_BATCH_PROCESSING:

                def chunked(seq_list, size):
                    for i in range(0, len(seq_list), size):
                        yield seq_list[i : i + size]

                for batch in chunked(novel_sequences, BATCH_SIZE):
                    batch_inputs: List[SequenceInput] = []
                    for seq in batch:
                        seq_str = str(seq.seq)
                        seq_data = sequences_data[seq_str]
                        assert seq_data.get("status") == "pending", (
                            f"Sequence {seq.id} is not pending"
                        )
                        assert seq_data.get("seq_uuid"), (
                            f"Sequence {seq.id} has no seq_uuid"
                        )
                        seq_uuid = seq_data.get("seq_uuid")
                        batch_inputs.append(
                            SequenceInput(sequence=seq, seq_uuid=seq_uuid)
                        )

                    logger.info(f"Spawning batch of {len(batch_inputs)} sequences")
                    call = process_16_sequences.spawn(batch_inputs)
                    job_ids.append(call.object_id)

                    # Update sequences map with job details for each sequence in the batch
                    for seq in batch:
                        seq_str = str(seq.seq)
                        seq_uuid = sequences_data[seq_str]["seq_uuid"]
                        sequences_data[seq_str] = SequenceData(
                            start_time=time.time(),
                            end_time=None,
                            seq_uuid=seq_uuid,
                            sequence_id=seq.id,
                            status="pending",
                            zip_path=None,
                            job_id=call.object_id,
                        )
            else:
                # Spawn single sequence jobs
                for seq in novel_sequences:
                    seq_str = str(seq.seq)
                    seq_data = sequences_data[seq_str]
                    assert seq_data.get("status") == "pending", (
                        f"Sequence {seq.id} is not pending"
                    )
                    assert seq_data.get("seq_uuid"), (
                        f"Sequence {seq.id} has no seq_uuid"
                    )
                    seq_uuid = seq_data.get("seq_uuid")
                    call = process_single_sequence.spawn(
                        SequenceInput(sequence=seq, seq_uuid=seq_uuid)
                    )
                    job_ids.append(call.object_id)
                    sequences_data[seq_str] = SequenceData(
                        start_time=time.time(),
                        end_time=None,
                        seq_uuid=seq_uuid,
                        sequence_id=seq.id,
                        status="pending",
                        zip_path=None,
                        job_id=call.object_id,
                    )
        else:
            logger.info("No novel sequences to process - all sequences are cached")

        total_sequence_count = len(parsed_sequences)
        novel_sequence_count = len(novel_sequences)
        cached_sequence_count = total_sequence_count - novel_sequence_count
        # Create .json metadata file for the run
        run_metadata = RunData(
            status="pending",
            sequences=sequences_data,
            fasta_filename=original_filename,
            output_filename=zip_output_filename,
            total_sequences=total_sequence_count,
            cached_sequences=cached_sequence_count,
            merged_zip_filename=None,
            submitted_at=time.time(),
            completed_at=None,
        )

        def validate_run_metadata(run_metadata: RunData) -> None:
            # TODO: validate JSON before writing to volume
            assert run_metadata["status"] == "pending", "Run status must be pending"
            # ... add remaining validation here

        validate_run_metadata(run_metadata)

        write_run_metadata_to_volume(run_id, run_metadata)

        logger.info(
            f"Submitted run {run_id}. {run_metadata['cached_sequences']}/{run_metadata['total_sequences']} sequences are cached"
        )

        vol.commit()
        return UploadFastaResponse(
            run_id=run_id,
            status="submitted" if novel_sequences else "ready",
            message=f"Successfully uploaded {original_filename} with {len(novel_sequences)}/{total_sequence_count} novel sequences.",
            job_ids=job_ids,
        )

    # TODO: Add stale check for failed jobs (by excessive time or job_id data.)
    @api.get("/status/{run_id}", summary="Check job status")
    async def check_status(run_id: str) -> StatusResponse:
        """Check the status of a processing job."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")

        try:
            json_path = get_run_json_path(run_id)
            if not json_path.exists():
                raise HTTPException(status_code=404, detail="Run not found")

            progress_data = get_run_metadata(run_id)
            sequences_data = progress_data["sequences"]

            # If the run is already complete we can short-circuit the work below.
            if progress_data["status"] == "complete":
                return StatusResponse(
                    run_id=run_id, status="complete", pending_sequences=[]
                )

            # 1. Determine which sequences are still pending.
            pending_seq_uuids = get_pending_uuids(sequences_data)

            # 2. Look in the cache directory for those expected zip files.
            completed_mapping = find_completed_in_cache(pending_seq_uuids)

            # 3. Update the in-memory metadata for any sequences we just discovered.
            pending_sequence_ids = update_sequences_with_completed(
                sequences_data, completed_mapping
            )

            if pending_sequence_ids:
                # Persist intermediate progress and return an in-progress response.
                write_run_metadata_to_volume(run_id, progress_data)
                vol.commit()
                return StatusResponse(
                    run_id=run_id,
                    status="pending",
                    pending_sequences=pending_sequence_ids,
                )

            # 4. If no pending sequences remain, mark the run as complete.
            progress_data["status"] = "complete"
            progress_data["completed_at"] = (
                max(completed_mapping.values()) if completed_mapping else time.time()
            )
            write_run_metadata_to_volume(run_id, progress_data)
            vol.commit()
            return StatusResponse(
                run_id=run_id, status="complete", pending_sequences=[]
            )

        except HTTPException:
            raise
        except Exception as e:
            logger.error(f"Error checking status for {run_id}: {e}")
            raise HTTPException(
                status_code=500, detail=f"Error checking status: {str(e)}"
            )

    @api.get("/download/{run_id}", summary="Download results")
    async def download_zip_file(run_id: str, output_filename: str = None):
        """Download the results zip file."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")

        try:
            merged_zip_path = None
            run_metadata = get_run_metadata(run_id)
            if output_filename:
                zip_output_filename = sanitize_output_filename(output_filename)
            else:
                zip_output_filename = run_metadata["output_filename"]

            assert zip_output_filename.endswith(".zip"), (
                "Output filename must end with .zip"
            )
            merged_zip_path = get_zip_by_fasta_dir() / f"{run_id}.zip"
            # (1) Build the merged archive if it does not yet exist.
            if not merged_zip_path.exists():
                sequences_data = run_metadata["sequences"]
                pending_uuids = get_pending_uuids(sequences_data)
                completed_mapping = find_completed_in_cache(pending_uuids)
                pending_sequence_ids = update_sequences_with_completed(
                    sequences_data, completed_mapping
                )
                if pending_sequence_ids:
                    raise HTTPException(
                        status_code=404,
                        detail="Results are not ready yet. Some sequences are still pending.",
                    )
                completed_zip_paths = get_completed_zip_paths(
                    sequences_data, require_all_complete=True
                )
                merged_zip_path = Path(
                    merge_zip_archives(completed_zip_paths, str(merged_zip_path))
                )

                run_metadata["merged_zip_filename"] = str(merged_zip_path)
                run_metadata["completed_at"] = time.time()
                run_metadata["status"] = "complete"
                write_run_metadata_to_volume(run_id, run_metadata)
                vol.commit()

            # (2) Read the merged archive and send it back to the client.
            zip_data = merged_zip_path.read_bytes()
            return Response(
                content=zip_data,
                media_type="application/zip",
                headers={
                    "Content-Disposition": f"attachment; filename={zip_output_filename}"
                },
            )

        except HTTPException:
            raise
        except Exception as e:
            logger.error(f"Error downloading zip for {run_id}: {e}")
            raise HTTPException(
                status_code=500, detail=f"Error downloading file: {str(e)}"
            )

    @api.get(
        "/retry/{run_id}",
        summary="Retry processing for sequences that are still pending",
    )
    async def retry_sequences(run_id: str):
        """Retry processing for sequences that are still pending."""
        if not run_id:
            raise HTTPException(status_code=400, detail="Run ID is required")
        # retry_pending_sequences.spawn(run_id)
        logger.info(type(RetryResponse))
        return RetryResponse(run_id=run_id, status="retry_submitted")

    return api
