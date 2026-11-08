# nardini-online: technical guide

[IMPLEMENTATION.md](IMPLEMENTATION.md) covers requirements, commands,
environment variables, tests and limits. This file is the API and storage reference.

## Processing flow

1. `POST /upload_fasta` parses the FASTA file with Biopython. A file without
   `>` headers is read as one sequence per line, named `sequence-<n>`.
2. Each unique sequence string is looked up in earlier runs on the volume. A
   completed one is reused (`cached`). One still pending in another run is
   not spawned again (`pending_external`).
3. Novel sequences are sorted longest first and spawned in batches of 16 to
   `process_16_sequences`, which runs NARDINI per sequence in a thread pool
   on the Python 3.9 worker image and writes one zip per sequence.
4. `GET /status` reports the run `complete` when no sequence is pending.
   `GET /download` merges the per-sequence zips into one run zip.
5. `GET /retry` respawns `retry_pending_sequences` for the run's pending
   sequences.

## Deploy

```bash
modal setup                  # or MODAL_TOKEN_ID / MODAL_TOKEN_SECRET
cd src && modal deploy app/backend.py
```

Modal prints the URL of the `fastapi_app` function. `src/app/common.py`
defines both images:

| Image | Python | Pins |
| --- | --- | --- |
| `web_image` | 3.12 | `fastapi[standard]==0.162.2`, `python-multipart==0.0.33`, `biopython==1.88` |
| `nardini_image` | 3.9 | `pydantic==2.13.5`, `numpy==1.19.5`, `nardini==1.1.1` |

`tests/test_pins.py` checks that the dev group in `pyproject.toml` and
`tests/requirements-worker.in` use the same pins.

## Endpoints

### `GET /health`

```json
{"status": "healthy"}
```

### `POST /upload_fasta`

Multipart form: `file` (required, `.fasta`, `.fa` or `.fas`) and
`output_filename` (optional query parameter, the name of the run zip).

```json
{
  "run_id": "uuid",
  "status": "submitted",
  "message": "...",
  "job_ids": ["modal-call-id"],
  "status_url": "/status/uuid",
  "download_url": null
}
```

`status` is `ready` when no sequence needed a new worker; `download_url` is
then `/download/uuid`. Both URLs are paths relative to the API URL.

### `GET /status/{run_id}`

```json
{
  "run_id": "uuid",
  "status": "pending",
  "pending_sequences": ["sequence id"],
  "download_url": null
}
```

`download_url` is set once `status` is `complete`.

### `GET /download/{run_id}`

Optional query parameter `output_filename`. Returns
`<fasta name>_results.zip` with each sequence's `.tsv` and `.png` files and
a `master_sequences.tsv` that uses the ids from the uploaded file.

### `GET /retry/{run_id}`

```json
{"run_id": "uuid", "status": "retry_submitted"}
```

## Status codes

Errors return `{"detail": "..."}`.

| Code | When |
| --- | --- |
| 400 | Wrong file extension, empty file, no sequences, FASTA parse error, malformed `run_id` |
| 404 | Unknown `run_id`; `/download` while sequences are still pending |
| 413 | Upload larger than `MAX_UPLOAD_MB` |
| 500 | Unexpected error in `/status` or `/download` |

## Volume layout

The Modal volume `VOLUME_NAME` is mounted at `/data`:

```
/data/
├── runs/{run_id}.json           # run metadata and per-sequence status
├── fasta_inputs/{run_id}/       # the uploaded FASTA file
└── zipfiles/
    ├── by_idr/{seq_uuid}.zip    # one NARDINI result per sequence
    └── by_fasta/{run_id}.zip    # merged run result
```
