# nardini-online: implementation

Built for the Woolford Lab at CMU, whose users go through the
[Colab notebook](https://colab.research.google.com/drive/1jKDu1AfOHI-4BM1hHELIDFQCqbPz86M9).
[GUIDE.md](GUIDE.md) is the API and storage reference.

## Requirements

- Python 3.11+ and [uv](https://docs.astral.sh/uv/) for the app and its offline tests.
- The real-worker test needs Python 3.9: NARDINI 1.1.1 requires `scipy<1.7`, which keeps the worker
  image on Python 3.9 and `numpy==1.19.5`.
- Deploying needs a Modal account. `modal setup` stores the credentials; `MODAL_TOKEN_ID` and
  `MODAL_TOKEN_SECRET` (see `.env.example`) work instead.

## Usage

| Command | What it does |
| --- | --- |
| `uv run pytest -q` | Offline tests: API contract, zip merging, helpers, image pins. Skips the real-worker test. |
| `uv run ruff check .` | Lint. |
| `uv run ruff format --check src tests` | Format check. |
| `cd src && modal deploy app/backend.py` | Deploy the app. Run it from `src/`, or `app.common` does not import. |

Replace `$URL` with the deployed app URL:

| Method and path | Returns |
| --- | --- |
| `GET /health` | `{"status": "healthy"}` |
| `POST /upload_fasta` (form field `file`) | `run_id`, `status` (`submitted` or `ready`), `status_url`, `download_url` |
| `GET /status/{run_id}` | `status` (`pending` or `complete`), `pending_sequences`, `download_url` |
| `GET /download/{run_id}` | `<fasta name>_results.zip` |
| `GET /retry/{run_id}` | Respawns the run's pending sequences. |

```bash
curl -F "file=@data/fasta_inputs/Halophile-pHtolerant-yeast-first16.fasta" "$URL/upload_fasta"
curl "$URL/status/$RUN_ID"
curl -OJ "$URL/download/$RUN_ID"
```

Read by `src/app/common.py` when you deploy:

| Variable | Default | Purpose |
| --- | --- | --- |
| `APP_NAME` | `nardini_online` | Modal app name. |
| `VOLUME_NAME` | `run_fasta_volume` | Modal volume for runs and cached results. Created if missing. |
| `TIMEOUT_SECONDS` | `21600` | Worker function timeout. |
| `MAX_UPLOAD_MB` | `10` | Upload limit; larger files get 413. |
| `DEPLOYED_URL` | `https://tangentleman--nardini-online-fastapi-app.modal.run` | API URL used by `src/config.py` and `notebooks/nardini_jupyter.ipynb`. |

## How it works

`src/app/backend.py` is the Modal app: the `fastapi_app` web function and the NARDINI workers
(`process_16_sequences`, `retry_pending_sequences`). The web function spawns workers in batches of
16 sequences and caches each sequence's result on the Modal volume, so a sequence seen in an
earlier run is not recomputed. GUIDE.md has the step-by-step flow.

`src/app/common.py` defines the Modal images, volume and limits. `src/shared_utils`,
`src/fastapi_utils` and `src/backend_utils` are the helpers for each image.
`tests/requirements-worker.in` holds the worker image's pins; `.txt` is the
`uv pip compile --python-version 3.9` lock.

The tests build the app with `fastapi_app.local()` and replace the Modal workers and volume with
fakes, so they need no Modal account and no network.

## Tests

```bash
uv sync --locked
uv run ruff check .
uv run ruff format --check src tests
uv run pytest -q
```

CI (`.github/workflows/ci.yml`) runs these on Python 3.11, 3.12 and 3.15. Its `worker` job runs
real NARDINI on two sample IDRs on the worker image's Python 3.9 stack and merges the zips:

```bash
uv venv -p 3.9 .venv-worker
uv pip install -p .venv-worker -r tests/requirements-worker.txt
MPLBACKEND=Agg .venv-worker/bin/python -m pytest tests/test_worker.py -q
```

## Limits

- The tests and the worker test use no network and no credentials.
- Each upload of new sequences spends Modal compute, and so does `/retry`.
- The deployed app is public: `REQUIRE_AUTH = False` in `src/app/backend.py` sets
  `requires_proxy_auth=False`.
- Not verified live: no deploy has been run from this code, including the web image on Python 3.12.
  The URL in `src/config.py` returned 404 on 2026-11-08. The Colab notebook in Drive has not been
  run; `notebooks/nardini_colab.ipynb` is a copy.

## Credits

Made by Tanuj Vasudeva and Ethan Caine, 2025. Thanks to Dr. John Woolford at Carnegie Mellon
University for supporting the project, Modal for hosting, and Katherine Parry for advice.

Cohan, M. C., Shinn, M. K., Lalmansingh, J. M., & Pappu, R. V. (2021). Uncovering non-random binary
patterns within sequences of intrinsically disordered proteins. *Journal of Molecular Biology*,
434(2), 167373.
