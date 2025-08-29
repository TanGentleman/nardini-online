# ONLY include: get_pending_uuids, find_completed_in_cache, update_sequences_with_completed, sanitize_output_filename
from typing import Dict, List

from shared_utils.file_utils import get_zip_by_idr_dir
from shared_utils.schemas import SequencesMapping

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
