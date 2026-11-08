"""Runs real NARDINI through the worker code. Needs the Python 3.9 worker env:

uv venv -p 3.9 .venv-worker
uv pip install -p .venv-worker -r tests/requirements-worker.txt
MPLBACKEND=Agg .venv-worker/bin/python -m pytest tests/test_worker.py -q
"""

import zipfile
from pathlib import Path

import pytest

pytest.importorskip("nardini")

from backend_utils.utils import _process_one_sequence_to_volume  # noqa: E402
from fastapi_utils.utils import (  # noqa: E402
    merge_zip_archives,
    read_sequences_from_filename,
)
from shared_utils.file_utils import get_zip_by_idr_dir  # noqa: E402

FASTA = (
    Path(__file__).parent.parent
    / "data/fasta_inputs/Halophile-pHtolerant-yeast-first16.fasta"
)


def test_worker_output_merges_into_one_zip(volume, tmp_path):
    records = read_sequences_from_filename(str(FASTA), "sequence")
    shortest = sorted(records, key=lambda r: len(r.seq))[:2]

    sequences = {}
    zips = []
    for i, record in enumerate(shortest):
        uuid = f"u{i}"
        path = _process_one_sequence_to_volume({"sequence": record, "seq_uuid": uuid})
        assert path == str(get_zip_by_idr_dir() / f"{uuid}.zip")
        with zipfile.ZipFile(path) as zf:
            assert set(zf.namelist()) == {
                f"regular-{record.id}.png",
                f"scrambled-{record.id}.png",
                f"zscore-original-sequence-{record.id}.tsv",
                f"zscore-scrambled-sequence-{record.id}.tsv",
                "sequences.tsv",
            }
        zips.append(path)
        sequences[str(record.seq)] = {
            "sequence_id": record.id,
            "status": "complete",
            "seq_uuid": uuid,
            "zip_path": None,
        }

    out = merge_zip_archives(zips, str(tmp_path / "merged.zip"), sequences)

    with zipfile.ZipFile(out) as zf:
        names = set(zf.namelist())
        master = zf.read("master_sequences.tsv").decode().splitlines()
    assert len(names) == 9
    for record in shortest:
        assert f"zscore-original-sequence-{record.id}.tsv" in names
    assert len([line for line in master if line.strip()]) == 3
