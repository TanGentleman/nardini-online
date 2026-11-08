import zipfile

from conftest import write_nardini_zip

from fastapi_utils.utils import merge_zip_archives
from shared_utils.file_utils import get_zip_by_idr_dir


def test_merge_renames_files_to_current_ids_and_builds_master_tsv(volume, tmp_path):
    idr_dir = get_zip_by_idr_dir()
    # A cached zip produced under an older id, and a fresh one from this run.
    write_nardini_zip(idr_dir / "old.zip", "OldName_IDR_1", "MKVLLE")
    write_nardini_zip(idr_dir / "u2.zip", "Drs1_IDR_2_[660-751]", "AAAA")
    sequences = {
        "MKVLLE": {
            "sequence_id": "Drs1_IDR_1_[0-221]",
            "status": "cached",
            "seq_uuid": None,
            "zip_path": str(idr_dir / "old.zip"),
        },
        "AAAA": {
            "sequence_id": "Drs1_IDR_2_[660-751]",
            "status": "complete",
            "seq_uuid": "u2",
            "zip_path": None,
        },
    }

    out = merge_zip_archives(
        [str(idr_dir / "old.zip"), str(idr_dir / "u2.zip")],
        str(tmp_path / "merged"),
        sequences,
    )

    assert out.endswith("merged.zip")
    with zipfile.ZipFile(out) as zf:
        names = set(zf.namelist())
        master = zf.read("master_sequences.tsv").decode().splitlines()
    assert "regular-Drs1_IDR_1_[0-221].png" in names
    assert "zscore-original-sequence-Drs1_IDR_1_[0-221].tsv" in names
    assert not any("OldName" in n for n in names)
    assert "sequences.tsv" not in names
    assert len(names) == 9
    assert master[0].split("\t")[0] == "ID"
    assert [line.split("\t")[0] for line in master[1:]] == [
        "Drs1_IDR_1_[0-221]",
        "Drs1_IDR_2_[660-751]",
    ]
    assert all(len(line.split("\t")) == 5 for line in master)
