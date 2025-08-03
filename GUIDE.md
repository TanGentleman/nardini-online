# RunFasta – Comprehensive Guide

## 1  What is RunFasta?

RunFasta is a cloud-backed, **blazing-fast interface for running the NARDINI statistical analysis** on Intrinsically Disordered Regions (IDRs) in proteins.  It offers two ways to drive the analysis:

1. **Interactive notebook** – Zero-install, guided workflow (ideal for researchers).
2. **REST API** – Upload a FASTA file, poll for completion, download a results ZIP (ideal for programmatic pipelines).

All heavy computation is executed remotely on [Modal](https://modal.com) using a scalable worker image; the local machine therefore only needs to send/receive files.

---

## 2  Quick Start (Notebook – Recommended)

1. **Open the notebook**  
   Either in local Jupyter _or_ directly on Google Colab:
   ```bash
   https://github.com/TanGentleman/RunFasta/blob/main/notebooks/demo.ipynb
   ```

2. **Install prerequisites**  (only needed locally – Colab already has them):
   ```bash
   pip install requests ipykernel
   ```

3. **Place your FASTA file** in one of:
   * `notebooks/` (same dir as notebook)
   * `data/fasta_inputs/`
   * a custom path you provide inside the notebook

4. **Run cells top-to-bottom**.  The notebook will:
   * sanity-check the backend health endpoint
   * let you pick a `.fasta` / `.fa` / `.fas` file
   * POST it to the `/upload_fasta` endpoint
   * expose a **Run ID** & persist run info
   * poll `/status/{run_id}` until finished
   * download `/download/{run_id}` into `data/zip_outputs/`

5. **Inspect results** – each ZIP contains:
   * `*.tsv` files – z-scores & statistics
   * `*.png` plots – visualisations
   * `master_sequences.tsv` – merged summary across all sequences in the upload

---

## 3  Deploying the Backend on Modal

> Skip this section if you only intend to **use** the public instance at  
> `https://tangentleman--nardini-backend-fastapi-app.modal.run`.

### 4.1  Prerequisites

* Modal account + billing set-up (free tier works fine)
* Modal CLI ≥ 0.54 ( `pip install modal` )
* Python 3.9 – 3.11

### 4.2  Steps

```bash
# 1. Authenticate (generates ~/.modal
modal token new

# 2. (Optional) set your own Modal volume name in src/modal_backend.py
#    by editing: vol = modal.Volume.from_name("nardini_volume", create_if_missing=True)

# 3. Deploy – this builds two container images:
#    • web_image     (FastAPI + upload logic)
#    • nardini_image (heavy worker with NARDINI + NumPy)
cd src
modal deploy modal_backend.py

# 4. Grab the URL printed by Modal → something like
#    https://<username>--nardini-backend-fastapi-app.modal.run
```

> **Tip:**  During development you can run `modal serve modal_backend.py` which hot-reloads code changes locally.

### 4.3  Updating
Simply re-run `modal deploy src/modal_backend.py` after editing the code.  The CLI uses caching so pushes are incremental.

---

## 4  REST API Reference

| Method & Path              | Purpose                                    | Body / Params                                                    |
|----------------------------|--------------------------------------------|------------------------------------------------------------------|
| `GET /health`              | Liveness probe                             | –                                                                |
| `POST /upload_fasta`       | Submit a FASTA file for analysis           | **form-data** `file=@your.fasta`  (≤ 10 MB by default)            |
| `GET /status/{run_id}`     | Poll progress                              | –  Returns JSON with per-sequence fields + `_status`             |
| `GET /download/{run_id}`   | Fetch merged results ZIP (once completed)  | –  Returns `application/zip` attachment                           |

### Example – cURL
```bash
# Upload a FASTA file
curl -X POST "${BACKEND_URL}/upload_fasta" \
     -F "file=@Halophile-pHtolerant-yeast-first16.fasta"

# → { "run_id": "<uuid>", ... }

# Check status
curl "${BACKEND_URL}/status/<uuid>"

# Download when ready
curl -L "${BACKEND_URL}/download/<uuid>" -o results.zip
```

### Example – Python Snippet
```python
import requests, time
url = "https://<your>.modal.run"

# 1. Upload
with open("sample.fasta", "rb") as f:
    r = requests.post(f"{url}/upload_fasta", files={"file": f})
    run_id = r.json()["run_id"]

# 2. Poll
while True:
    status = requests.get(f"{url}/status/{run_id}").json()["status"]
    if status == "completed":
        break
    time.sleep(10)

# 3. Download
zip_bytes = requests.get(f"{url}/download/{run_id}").content
open("results.zip", "wb").write(zip_bytes)
```

---

## 5  How It Works (Under the Hood)

```
┌──────────┐  POST /upload_fasta          ┌────────────┐   fork workers    ┌──────────┐
│  Client  │─────────────────────────────▶│  FastAPI    │──────────────────▶│ Process  │
└──────────┘                              │  (web_image)│17×parallel seq.   │  Pool    │
    ▲                                     └─────┬──────┘                   └────┬─────┘
    │   GET /status/{run_id}                    │   spawn process_nardini_job     │
    │   GET /download/{run_id}                  ▼                                │
┌──────────┐                              ┌────────────┐   NARDINI / NumPy   │
│  Browser │◀─────────────────────────────│ process_*  │─────────────────────┘
└──────────┘                              └────────────┘
```

1. **FastAPI layer** (`fastapi_app`) lives in a lightweight container (`web_image`) and simply enqueues work.
2. **Heavy worker** (`process_nardini_job`) lives in `nardini_image` (Python 3.9 + NARDINI + NumPy).
3. Worker parses the FASTA, skips sequences already cached in a **Modal Volume** (`/data`).
4. Novel sequences are executed in parallel via a `ProcessPoolExecutor` (default 16 cores).
5. Each sequence run produces its own `nardini-data-*.zip`; once all finish they are merged (`mergeZips`) into one archive per FASTA upload.
6. Progress is checkpointed to `<run_id>_in_progress.json` inside the volume for real-time status polling.

---

## 6  Configuration & Customisation

| Setting | Location | Default | Description |
|---------|----------|---------|-------------|
| **MAX_FILE_SIZE** | `fastapi_app.upload_fasta` | 10 MB | Hard cap on upload size |
| **NUM_SCRAMBLED_SEQUENCES** | NARDINI constant | 100k (see library) | Controls statistical rigour |
| **ProcessPool max_workers** | `process_nardini_job` | 16 | Increase for larger Modal CPUs |
| **Volume name** | `modal_backend.py` line 54 | `nardini_volume` | Rename if sharing projects |

> **Changing Nardini parameters** – edit the `TYPEALL`, `NUM_SCRAMBLED_SEQUENCES`, etc. import section and redeploy.

---

## 7  Caching Strategy

* **Sequence-level cache** – each unique amino-acid string gets its own ZIP.  Hash → ZIP mapping is stored as JSON files in `/data/seq_zip_maps/` on the Modal volume.
* **Fasta-level cache** – merged ZIPs live under `/data/zipfiles/by_fasta/` keyed by `run_id`.
* Re-uploading a FASTA with mixtures of cached/novel sequences only recomputes the new ones, dramatically saving time.

---

## 8  Troubleshooting

| Symptom | Probable Cause | Fix |
|---------|----------------|-----|
| `413 Payload Too Large` on upload | File > 10 MB | Split FASTA or raise `MAX_FILE_SIZE` constant & redeploy |
| `/status` returns `failed` | Nardini internal error (rare invalid sequence) | Check `error` field in JSON  → remove problematic sequence |

---

## 9  Optimizing the Engine (Parallelism & Performance)

RunFasta already parallelises per-sequence calculations using a Python
`ProcessPoolExecutor`.  Here are some levers you can tune to squeeze
the most performance out of your Modal deployment:

1. **Number of workers (`max_workers`)**  
   `process_nardini_job` currently runs with `max_workers=16`.  
   • This means a queue is constructed with a .fasta file with over 16 sequences.
   • This can be optimized by using other techniques for batching. See (https://modal.com/docs/guide/scale)

---



