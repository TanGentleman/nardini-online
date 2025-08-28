import os
import modal
# Configurable parameters with env-variable overrides
APP_NAME = os.getenv("APP_NAME", "nardini_online")
VOLUME_NAME = os.getenv("VOLUME_NAME", "run_fasta_volume")
TIMEOUT_SECONDS = int(os.getenv("TIMEOUT_SECONDS", "21600"))  # default 6h
MAX_UPLOAD_MB = int(os.getenv("MAX_UPLOAD_MB", "10"))
MAX_FILE_SIZE = MAX_UPLOAD_MB * 1024 * 1024  # bytes
# ---------------------------------------------------------------------- #
# Create a Modal application with shared volume
app = modal.App(APP_NAME)
vol = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)
# Lightweight web image for FastAPI
web_image = (
    modal.Image.debian_slim()
    .uv_pip_install(
        "fastapi[standard]==0.115.13",
        "python-multipart==0.0.20",
        "biopython==1.84",
    )
    .add_local_dir("shared_utils", remote_path="/root/shared_utils")
    .add_local_dir("fastapi_utils", remote_path="/root/fastapi_utils")
)

# Heavy processing image for Nardini
nardini_image = (
    modal.Image.debian_slim(python_version="3.9")
    .uv_pip_install(
        "pydantic==2.10.6",
        "numpy==1.19.5",
    )
    .uv_pip_install("nardini==1.1.1")
    .add_local_dir("shared_utils", remote_path="/root/shared_utils")
    .add_local_dir("backend_utils", remote_path="/root/backend_utils")
)