import zipfile
from pathlib import Path

import pytest

import shared_utils.file_utils as file_utils


@pytest.fixture
def volume(tmp_path, monkeypatch):
    """Point the Modal volume paths at a temp dir and create its folders."""
    monkeypatch.setattr(file_utils, "get_volume_dir", lambda: tmp_path)
    file_utils.ensure_volume_directories()
    return tmp_path


def write_nardini_zip(path: Path, seq_id: str, seq: str) -> Path:
    """Write a zip shaped like NARDINI 1.1.1 output for one sequence.

    File names and the space-aligned sequences.tsv match a real
    calculate_zscore_and_plot run.
    """
    header = (
        "ID  original_seq  most_similar_seq  "
        "sum_abs_zscore_original_seq  sum_abs_zscore_scrambled_seq"
    )
    row = f"{seq_id}  {seq}  {seq[::-1]}  0.734527  0.452121"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(f"regular-{seq_id}.png", b"png")
        zf.writestr(f"scrambled-{seq_id}.png", b"png")
        zf.writestr(f"zscore-original-sequence-{seq_id}.tsv", "z\n")
        zf.writestr(f"zscore-scrambled-sequence-{seq_id}.tsv", "z\n")
        zf.writestr("sequences.tsv", f"{header}\n{row}\n")
    return path
