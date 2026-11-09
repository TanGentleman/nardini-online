# Roadmap

Nardini Online runs NARDINI (IDR patterning statistics) on a FASTA file as a
Modal-hosted REST API. Woolford Lab biologists use it through a Colab
notebook; nobody installs anything.

## Status (checked 2026-11-08 against upstream `ddcf372`)

- Live API is down: `/health` on
  `https://tangentleman--nardini-online-fastapi-app.modal.run` returns 404.
- The API runs offline: `fastapi_app.local()` builds the FastAPI app without
  Modal credentials, and `uv run pytest` runs the API contract against fake
  workers and a temp-dir volume.
- Deploy works only from `src/` (`cd src && modal deploy app/backend.py`).
- `tests/test_worker.py` runs real NARDINI 1.1.1 on two sample IDRs on the
  Python 3.9 worker stack and merges the zips (about 2 minutes).
- Dependencies (`-7`): web image and dev group at fastapi 0.162.2,
  python-multipart 0.0.33, biopython 1.88, web image on Python 3.12; worker
  at pydantic 2.13.5 and pytest 8.4.2 (the last releases for 3.9). NARDINI
  1.1.1 requires `scipy<1.7`, which keeps the worker on Python 3.9 and
  `numpy==1.19.5`.

## Signal vs noise

Keep: `src/app/` (Modal app), `shared_utils`, `fastapi_utils`,
`backend_utils`, the Colab notebook copy, the sample FASTA.

Cut (`-2`): `src/single-container-app/` and its notebook (older app on a
separate volume, not what the notebook calls), `migrate_run_metadata`
(one-off), the `USE_BATCH_PROCESSING`/`REPLACE_IDR_RANGES` flags that are
always off, duplicate helpers in `fastapi_utils`.

## Stack

- `nardini-online-1` (done): import, uv + lockfile, offline tests (API
  contract with a fake worker), CI, plain README, GUIDE deploy fix, Jupyter
  status fix.
- `nardini-online-2` (done): cut the noise above (about 1,300 lines).
- `nardini-online-3` (done): real NARDINI worker test on Python 3.9, own CI job.
- `nardini-online-4` (done): 413 for oversize uploads, 404 for unknown runs
  on `/download`, uploaded ids in `master_sequences.tsv`.
- `nardini-online-5` (done): make `/retry` respawn pending sequences; 400/404
  for bad or unknown run ids instead of a blind `retry_submitted`.
- `nardini-online-6` (done): issue #11, `status_url` and `download_url` in
  the upload and status responses (new fields only).
- `nardini-online-7`: current deps and toolchain (Python 3.11+, uv 0.13.7,
  pinned actions), warnings as errors, a fix for the unclosed file Biopython
  1.88 leaves in `SeqIO.parse(path)`, a test that the offline pins match the
  Modal images, README and GUIDE rewritten.
- `nardini-online-8`: README cut to 27 words; details moved to IMPLEMENTATION.md.

## Needs Tan

- Redeploy (`cd src && modal deploy app/backend.py`) and check `/health`;
  the lab's notebook cannot reach the API right now.
- GitHub Actions only runs once applied upstream (it lives in
  `nardini-online/.github` here).
- Whether the deployed app should keep `requires_proxy_auth=False`
  (public, unauthenticated uploads).
- `-4` changes error codes (oversize 400 -> 413, unknown download 500 -> 404)
  and the ids in `master_sequences.tsv`. Check that no lab script depends on
  the old behavior.
- `-5` makes `/retry` spend Modal compute; a sequence whose worker is still
  running gets a second one. Bad or unknown run ids now get 400/404.
- `-6` returns relative paths (`/status/<id>`), not absolute URLs, because I
  could not check what scheme Modal's proxy forwards.
- `-2` removes `modal_backend_stable.py`; if `nardini-backend-stable` is
  still deployed on Modal, it can be stopped.
- The Colab notebook lives in Drive; only the repo copy was checked.
- `-7` moves the web image from fastapi 0.115 to 0.162 and Python 3.12.
  Tested offline only; check `/health` and one upload after redeploying.
- The `worker` CI job has not run on GitHub yet (about 2 minutes locally).
