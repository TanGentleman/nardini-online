# ## Setup and Configuration
# Import required packages for NARDINI analysis

import datetime
import logging
import os
from pathlib import Path
import json
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
    headers = _get_auth_headers()
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
    if not Path(fasta_filepath).exists():
        raise FileNotFoundError(f"File {fasta_filepath} does not exist")

    headers = _get_auth_headers()
    with open(fasta_filepath, "rb") as f:
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
    headers = _get_auth_headers()
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

    destination_dir = Path(destination_dir)
    if not destination_dir.exists():
        raise FileNotFoundError(
            f"Destination directory {destination_dir} does not exist"
        )

    headers = _get_auth_headers()
    try:
        response = requests.get(
            f"{url}/download/{run_id}", stream=True, headers=headers
        )
        if response.ok:
            # Extract filename from response headers
            content_disposition = response.headers.get("content-disposition", "")
            if "filename=" in content_disposition:
                filename = content_disposition.split("filename=")[1].strip('"')
                destination_filepath = destination_dir / filename
            else:
                destination_filepath = destination_dir / f"{run_id}.zip"

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

    headers = _get_auth_headers()
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
def save_run_info(run_id: str, fasta_filename: str, output_filepath: Path | str | None = None):
    """Save run information to a JSON file for reference."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    output_filepath = Path(output_filepath)
    if not output_filepath.parent.exists():
        raise FileNotFoundError(
            f"Destination directory {output_filepath.parent} does not exist"
        )

    run_info = {
        "title": "NARDINI Analysis Run Information",
        "timestamp": timestamp,
        "fasta_file": fasta_filename,
        "run_id": run_id
    }

    with open(output_filepath, "w") as f:
        json.dump(run_info, f, indent=2)

    return str(output_filepath)

def get_available_runs(json_path: Path | str):
    """Get all available runs from the JSON file."""
    if not Path(json_path).exists():
        return []
    with open(json_path, "r") as f:
        return json.load(f)