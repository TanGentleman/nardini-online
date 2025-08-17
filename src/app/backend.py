import os
import tempfile
from pathlib import Path
from typing import List
from uuid import uuid4
import logging
import time
import modal

from app.common import app, nardini_image, web_image, vol, MAX_FILE_SIZE, TIMEOUT_SECONDS
from shared_utils.schemas import SequenceInput
from shared_utils.file_utils import get_volume_dir, write_run_metadata_to_volume, get_run_json_path, get_run_metadata, get_zip_by_fasta_dir, ensure_volume_directories, get_zip_by_idr_dir, get_runs_dir
from shared_utils.utils import get_pending_uuids, find_completed_in_cache, update_sequences_with_completed, sanitize_output_filename
logger = logging.getLogger(__name__)

VOLUME_DIR = get_volume_dir()

@app.function(
    image=nardini_image,
    volumes={str(VOLUME_DIR): vol},
    timeout=TIMEOUT_SECONDS,
    scaledown_window=5,
)
def retry_pending_sequences(run_id: str):
    """Retry processing for sequences that are still pending."""
    from backend_utils.utils import create_seqrecord
    from shared_utils.schemas import RunData, SequenceData
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
    from backend_utils.utils import _process_one_sequence_to_volume
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
    from backend_utils.utils import _process_one_sequence_to_volume
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


REQUIRE_AUTH = False
@app.function(
    image=web_image,
    volumes={str(VOLUME_DIR): vol},
    # cpu=0.125,
    # memory=100,
)
@modal.asgi_app(requires_proxy_auth=REQUIRE_AUTH)
def fastapi_app():
    from fastapi import FastAPI, File, HTTPException, UploadFile # type: ignore
    from fastapi.responses import Response # type: ignore

    from shared_utils.schemas import (
        HealthResponse,
        RetryResponse,
        StatusResponse,
        UploadFastaResponse,
        RunData,
        SequenceData
    )

    from fastapi_utils.utils import read_sequences_from_filename, create_sequences_data, get_completed_zip_paths, merge_zip_archives

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
