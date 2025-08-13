import os
from pathlib import Path

# Configuration
ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "data"
INPUTS_DIR = DATA_DIR / "fasta_inputs"
OUTPUTS_DIR = DATA_DIR / "zip_outputs"


def get_backend_url(dev=False):
    deployed_url = os.getenv(
        "DEPLOYED_URL",
        "https://tangentleman--nardini-backend-dev-fastapi-app.modal.run",
    )
    if dev:
        dev_url = deployed_url.replace(".modal.run", "-dev.modal.run")
        return dev_url
    return deployed_url


BACKEND_URL = get_backend_url(dev=False)


def ensure_dirs():
    if not INPUTS_DIR.exists():
        INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    if not OUTPUTS_DIR.exists():
        OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
