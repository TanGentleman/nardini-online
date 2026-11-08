# Roadmap

Nardini Online runs NARDINI (IDR patterning statistics) on a FASTA file as a
Modal-hosted REST API. Woolford Lab biologists use it through a Colab
notebook; nobody installs anything.

## Status (checked 2026-11-08 against upstream `ddcf372`)

- Live API is down: `curl https://tangentleman--nardini-online-fastapi-app.modal.run/health`
  returns `404 modal-http: invalid function call` (same for the `-dev` and
  `nardini-backend-stable` URLs). The Colab notebook link itself loads (200).
- Code imports and the API runs offline: `fastapi_app.local()` builds the
  FastAPI app without Modal credentials. `uv run pytest` runs 25 tests on
  it as of `-6` (the real-worker test is skipped unless NARDINI is
  installed).
- `modal deploy src/app/backend.py` from the repo root (as GUIDE.md said)
  fails with `ModuleNotFoundError: No module named 'app'`; from `src/`,
  `modal deploy app/backend.py` gets as far as the auth check.
- NARDINI 1.1.1 still installs on Python 3.9 with `numpy==1.19.5` (the
  worker image pins) and produces the zip layout the merge code expects:
  `tests/test_worker.py` runs the worker on two sample IDRs and merges them
  (about 2.5 minutes; see its docstring for the commands).
- Dependencies: `modal>=1.1.1` resolves to 1.6.3 and works. Web image pins
  `fastapi==0.115.13` (latest 0.162), `biopython==1.84` (1.88),
  `python-multipart==0.0.20` (0.0.33). Left as is: the deployed image is the
  contract. `uv.lock` was gitignored; pyproject was named `runFasta`.
- Bugs found by the tests (fixed on `-4` unless noted):
  - an oversize upload returns 400, not the documented 413 (the parse handler
    rewraps every HTTPException);
  - `/download/<unknown run>` returns 500 with the volume path in the message;
  - `/retry` only returns `retry_submitted`: the spawn is commented out,
    `retry_pending_sequences` imports a `create_seqrecord` that does not
    exist and never saves what it spawns (fixed on `-5`);
  - `master_sequences.tsv` keeps the old sequence id for cached sequences
    while the per-file names get the new one.
- `notebooks/nardini_jupyter.ipynb` checked for `in_progress`/`completed`,
  which the API never returns (fixed on `-1`).

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

## Apply upstream

```bash
# in kranti-sandbox: patches after the import commit
base=$(git log --format=%H -1 --grep='^Import nardini-online upstream' origin/nardini-online-1)
git format-patch --relative=nardini-online "$base"..origin/nardini-online-1 -o /tmp/n1
# in a clean clone of TanGentleman/nardini-online at ddcf372
git am /tmp/n1/*.patch
uv sync --locked && uv run pytest
```

For `-N`, use `origin/nardini-online-(N-1)..origin/nardini-online-N`.

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
- The Colab notebook lives in Drive; only the repo copy was checked. The web
  image pins (`fastapi==0.115.13` etc.) were left as deployed.
- The `worker` CI job has not run on GitHub yet (about 3 minutes locally).
