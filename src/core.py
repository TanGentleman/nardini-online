# ## Setup and Configuration
# Import required packages for NARDINI analysis

import datetime
import json
import logging
import os
from pathlib import Path

import requests
from dotenv import load_dotenv

from shared_utils.schemas import (
    ErrorResponse,
    HealthResponse,
    RetryResponse,
    SimplifiedDownloadResponse,
    StatusResponse,
    UploadFastaResponse,
)

# Set up logger
logger = logging.getLogger(__name__)

REQUIRE_AUTH = False


def sanitize_path(user_path: str, base_dir: Path) -> Path:
    """Sanitize user path to prevent traversal."""
    safe_path = Path(user_path).resolve()
    if not safe_path.is_relative_to(base_dir):
        raise ValueError("Path traversal detected")
    return safe_path


def _get_auth_headers():
    load_dotenv()
    token_id = os.getenv("MODAL_TOKEN_ID")
    token_secret = os.getenv("MODAL_TOKEN_SECRET")
    if not token_id or not token_secret:
        raise ValueError("MODAL_TOKEN_ID and MODAL_TOKEN_SECRET must be set")
    return {"Modal-Key": token_id, "Modal-Secret": token_secret}


# Test the health endpoint
def test_health(url: str):
    """Test if the NARDINI backend service is healthy."""
    headers = _get_auth_headers() if REQUIRE_AUTH else None
    try:
        health_response = requests.get(f"{url}/health", headers=headers)
        if health_response.ok:
            res = health_response.json()
            return HealthResponse(status=res["status"])
        else:
            return ErrorResponse(
                error=f"Error: {health_response.status_code} {health_response.text}"
            )
    except Exception as e:
        return ErrorResponse(error=f"Connection error: {e}")


# Main function to run Nardini
# TODO: Add param for output_filename
def upload_fasta(
    url: str, fasta_filepath: Path | str
) -> UploadFastaResponse | ErrorResponse:
    """Submit a FASTA file for NARDINI analysis."""
    safe_path = sanitize_path(str(fasta_filepath), Path.cwd())
    if not safe_path.exists():
        raise FileNotFoundError(f"File {safe_path} does not exist")

    headers = _get_auth_headers() if REQUIRE_AUTH else None
    with open(safe_path, "rb") as f:
        files = {"file": f}
        response = requests.post(f"{url}/upload_fasta", files=files, headers=headers)
    if response.ok:
        res = response.json()
        return UploadFastaResponse(
            run_id=res["run_id"],
            status=res["status"],
            message=res["message"],
            job_ids=res["job_ids"],
        )
    else:
        return ErrorResponse(error=f"Error: {response.status_code} {response.text}")


def get_run_status(url: str, run_id: str):
    """Check the status of a NARDINI analysis run."""
    headers = _get_auth_headers() if REQUIRE_AUTH else None
    try:
        status_response = requests.get(f"{url}/status/{run_id}", headers=headers)
        if status_response.ok:
            res = status_response.json()
            return StatusResponse(
                run_id=res["run_id"],
                status=res["status"],
                pending_sequences=res["pending_sequences"],
            )
        else:
            return f"Error: {status_response.status_code} {status_response.text}"
    except Exception as e:
        return f"Connection error: {e}"


# TODO: Add param for output_filename
def download_zip(
    url: str, run_id: str, destination_dir: Path | str
) -> SimplifiedDownloadResponse | ErrorResponse:
    """Download the results zip file for a completed analysis."""
    if not run_id:
        return ErrorResponse(error="Please provide a valid Run ID.")

    safe_dest = sanitize_path(str(destination_dir), Path.cwd())
    if not safe_dest.exists():
        raise FileNotFoundError(f"Destination directory {safe_dest} does not exist")

    headers = _get_auth_headers() if REQUIRE_AUTH else None
    try:
        response = requests.get(
            f"{url}/download/{run_id}", stream=True, headers=headers
        )
        if response.ok:
            # Extract filename from response headers
            content_disposition = response.headers.get("content-disposition", "")
            if "filename=" in content_disposition:
                filename = content_disposition.split("filename=")[1].strip('"')
                destination_filepath = safe_dest / filename
            else:
                destination_filepath = safe_dest / f"{run_id}.zip"

            with open(destination_filepath, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            logger.info(f"Downloaded results to: {destination_filepath}")
            logger.info(
                f"File size: {destination_filepath.stat().st_size / (1024 * 1024):.1f} MB"
            )
            return SimplifiedDownloadResponse(
                run_id=run_id, destination_filepath=str(destination_filepath)
            )
        else:
            logger.info("Analysis is likely still in progress!")
            return ErrorResponse(
                error=f"Error downloading file: {response.status_code} {response.text}"
            )
    except Exception as e:
        return ErrorResponse(error=f"Download error: {e}")


def retry_sequences(url: str, run_id: str):
    """Retry processing for sequences that are still pending."""
    if not run_id:
        return ErrorResponse(error="Please provide a valid Run ID.")

    headers = _get_auth_headers() if REQUIRE_AUTH else None
    try:
        response = requests.get(f"{url}/retry/{run_id}", headers=headers)
        if response.ok:
            res = response.json()
            return RetryResponse(run_id=res["run_id"], status=res["status"])
        else:
            return ErrorResponse(
                error=f"Error retrying sequences: {response.status_code} {response.text}"
            )
    except Exception as e:
        return ErrorResponse(error=f"Retry error: {e}")


# TODO: Move this to a JSON file on client-side, or associate runs with users in Modal Volume
def save_run_info(
    run_id: str, fasta_filename: str, output_filepath: Path | str | None = None
):
    """Save run information to a JSON file for reference."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    safe_output = sanitize_path(str(output_filepath), Path.cwd())
    if not safe_output.parent.exists():
        raise FileNotFoundError(
            f"Destination directory {safe_output.parent} does not exist"
        )

    run_info = {
        "title": "NARDINI Analysis Run Information",
        "timestamp": timestamp,
        "fasta_file": fasta_filename,
        "run_id": run_id,
    }

    with open(safe_output, "w") as f:
        json.dump(run_info, f, indent=2)

    return str(safe_output)


def get_available_runs(json_path: Path | str):
    """Get all available runs from the JSON file."""
    safe_json = sanitize_path(str(json_path), Path.cwd())
    if not safe_json.exists():
        return []
    with open(safe_json, "r") as f:
        return json.load(f)
