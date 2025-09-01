"""
TypedDict schemas for FastAPI endpoints and backend metadata.
These mirror the structures used in `modal_backend.py`.
"""

from typing import Any, Dict, List, Optional

from typing_extensions import Literal, TypedDict

# ---- Public API response schemas ----


class ErrorResponse(TypedDict):
    error: str


class HealthResponse(TypedDict):
    status: Literal["healthy"]


class UploadFastaResponse(TypedDict):
    run_id: str
    status: Literal["submitted", "ready"]
    message: str
    job_ids: List[str]


class StatusResponse(TypedDict):
    run_id: str
    status: Literal["pending", "complete"]
    pending_sequences: List[str]


class RetryResponse(TypedDict):
    run_id: str
    status: Literal["retry_submitted"]


class SimplifiedDownloadResponse(TypedDict):
    run_id: str
    destination_filepath: str


# ---- Backend metadata schemas (stored in volume as JSON) ----


class SequenceInput(TypedDict):
    sequence: Any  # Bio.SeqRecord object
    seq_uuid: str


class SequenceData(TypedDict):
    sequence_id: str
    status: Literal["pending", "cached", "pending_external", "complete"]
    start_time: Optional[float]
    end_time: Optional[float]
    seq_uuid: Optional[str]
    zip_path: Optional[str]
    job_id: Optional[str]


SequenceString = str
SequencesMapping = Dict[SequenceString, SequenceData]


class RunData(TypedDict):
    status: Literal["pending", "complete"]
    fasta_filename: str
    output_filename: str
    sequences: SequencesMapping
    total_sequences: int
    cached_sequences: int
    merged_zip_filename: Optional[str]
    submitted_at: float
    completed_at: Optional[float]
