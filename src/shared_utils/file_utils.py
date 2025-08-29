import logging
import json
from pathlib import Path

from shared_utils.schemas import RunData

def clean_filename(filename: str) -> str:
    """Sanitize a filename to prevent path traversal."""
    return Path(filename).name

def get_volume_dir() -> Path:
    return Path("/data")

def get_runs_dir() -> Path:
    """Absolute path to the directory storing per-run JSON metadata files."""
    return get_volume_dir() / "runs"


def get_zip_by_idr_dir() -> Path:
    """Absolute path to per-sequence zip outputs, keyed by `seq_uuid`."""
    return get_volume_dir() / "zipfiles" / "by_idr"


def get_zip_by_fasta_dir() -> Path:
    """Absolute path to merged zip outputs, keyed by `run_id`."""
    return get_volume_dir() / "zipfiles" / "by_fasta"

def get_fasta_dir() -> Path:
    """Absolute path to per-sequence zip outputs, keyed by `seq_uuid`."""
    return get_volume_dir() / "fasta_inputs"


def get_run_json_path(run_id: str) -> Path:
    return get_runs_dir() / f"{clean_filename(run_id)}.json"


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

def ensure_volume_directories() -> bool:
    """Ensure the required directories exist in the mounted volume."""
    was_created = False
    for directory in [get_zip_by_fasta_dir(), get_zip_by_idr_dir(), get_runs_dir(), get_fasta_dir()]:
        if not directory.exists():
            logging.warning(f"Volume directory {directory} does not exist, creating it")
            directory.mkdir(parents=True, exist_ok=True)
            was_created = True
    return was_created


def save_fasta_to_volume(content: str, filename: str, run_id: str) -> None:
    path = get_fasta_dir() / clean_filename(run_id) / clean_filename(filename)
    if not path.parent.exists():
        logging.info(f"Saving fasta to folder {path.parent}")
        path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(content)
