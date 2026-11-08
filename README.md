# Nardini Online

Runs [NARDINI](https://github.com/mshinn23/nardini) (Cohan et al., 2021) on a
FASTA file of intrinsically disordered regions (IDRs) and returns one zip of
`.tsv` and `.png` results. Built for the Woolford Lab at CMU. Sequences run in
parallel on [Modal](https://modal.com), and results are cached per sequence,
so a sequence seen before is not recomputed.

## Use it

Open the [Colab notebook](https://colab.research.google.com/drive/1jKDu1AfOHI-4BM1hHELIDFQCqbPz86M9),
upload a FASTA file, and run the cells. You can close the notebook while the
analysis runs and come back for the results.

Input looks like [`data/fasta_inputs/Halophile-pHtolerant-yeast-first16.fasta`](data/fasta_inputs/Halophile-pHtolerant-yeast-first16.fasta):

```
>Drs1_IDR_2_[660-751]
MLKHKKEIQARPRRTWFQSESDKKNSKVLGALSRNKKVTNSKKRKREEAKADGNGARSYRKTKTDRIADQERTFKKQKSTNSNKKKGFKSR
```

## REST API

The notebook calls this API. Replace `$URL` with the deployed app URL.

```bash
curl -F "file=@sequences.fasta" "$URL/upload_fasta"   # -> {"run_id": ..., "status": "submitted" | "ready", ...}
curl "$URL/status/$RUN_ID"                             # -> {"status": "pending" | "complete", "pending_sequences": [...]}
curl -OJ "$URL/download/$RUN_ID"                       # -> <fasta name>_results.zip
```

See [GUIDE.md](GUIDE.md) for every endpoint and deployment.

## Develop

```bash
uv sync
uv run pytest            # offline: fakes Modal and the NARDINI worker
cd src && modal deploy app/backend.py
```

- `src/app/backend.py`: the Modal app (FastAPI web function and NARDINI workers)
- `src/app/common.py`: Modal images, volume and limits
- `src/shared_utils`, `src/fastapi_utils`, `src/backend_utils`: helpers for each image
- `notebooks/nardini_colab.ipynb`: copy of the Colab notebook

## Credits

Made by Tanuj Vasudeva and Ethan Caine, 2025. Thanks to Dr. John Woolford at
Carnegie Mellon University for supporting the project, Modal for hosting, and
Katherine Parry for advice.

Cohan, M. C., Shinn, M. K., Lalmansingh, J. M., & Pappu, R. V. (2021).
Uncovering non-random binary patterns within sequences of intrinsically
disordered proteins. *Journal of Molecular Biology*, 434(2), 167373.

MIT License.
