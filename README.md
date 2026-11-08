# nardini-online

Runs [NARDINI](https://github.com/mshinn23/nardini) 1.1.1 on a FASTA file of disordered regions
(IDRs) and returns one zip of results. A FastAPI app on Modal; deploying needs a Modal account.

```bash
git clone https://github.com/TanGentleman/nardini-online && cd nardini-online
uv sync --locked
uv run pytest -q    # offline API tests, no Modal account
```

Needs Python 3.11+ and [uv](https://docs.astral.sh/uv/). Details:
[IMPLEMENTATION.md](IMPLEMENTATION.md), API reference: [GUIDE.md](GUIDE.md), plan:
[ROADMAP.md](ROADMAP.md).

MIT licensed.
