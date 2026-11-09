import pytest

from fastapi_utils.utils import create_sequences_data, read_sequences_from_filename
from shared_utils.file_utils import get_zip_by_idr_dir, write_run_metadata_to_volume
from shared_utils.utils import (
    find_completed_in_cache,
    get_pending_uuids,
    sanitize_output_filename,
    update_sequences_with_completed,
)


@pytest.mark.parametrize(
    "name, expected",
    [
        ("yeast.fasta", "yeast.zip"),
        ("yeast.fa", "yeast.zip"),
        ("yeast.fas", "yeast.zip"),
        ("results.zip", "results.zip"),
        ("my results", "my results.zip"),
        ("../../etc/passwd", "passwd.zip"),
        ("/abs/dir/run.fasta", "run.zip"),
    ],
)
def test_sanitize_output_filename(name, expected):
    assert sanitize_output_filename(name) == expected


def test_sanitize_output_filename_truncates_and_requires_name():
    assert len(sanitize_output_filename("a" * 400)) == 254
    with pytest.raises(ValueError):
        sanitize_output_filename("")


def test_read_fasta(tmp_path):
    path = tmp_path / "in.fasta"
    path.write_text(">IDR_1_[0-5]\nMKV\nLLE\n>IDR_2\nAAAA\n")
    records = read_sequences_from_filename(str(path), "sequence")
    assert [(r.id, str(r.seq)) for r in records] == [
        ("IDR_1_[0-5]", "MKVLLE"),
        ("IDR_2", "AAAA"),
    ]


def test_read_raw_sequences_get_default_names(tmp_path):
    path = tmp_path / "in.fasta"
    path.write_text("MKV\n\nAAAA\r\n")
    records = read_sequences_from_filename(str(path), "sequence")
    assert [(r.id, str(r.seq)) for r in records] == [
        ("sequence-1", "MKV"),
        ("sequence-2", "AAAA"),
    ]


def test_read_missing_file_raises(tmp_path):
    with pytest.raises(RuntimeError):
        read_sequences_from_filename(str(tmp_path / "nope.fasta"), "sequence")


def _seq(status, seq_uuid=None, zip_path=None, sequence_id="s"):
    return {
        "sequence_id": sequence_id,
        "status": status,
        "start_time": None,
        "end_time": None,
        "seq_uuid": seq_uuid,
        "zip_path": zip_path,
        "job_id": None,
    }


def test_create_sequences_data_reuses_prior_runs(volume):
    write_run_metadata_to_volume(
        "11111111-1111-1111-1111-111111111111",
        {
            "sequences": {
                "DONE": _seq("complete", zip_path="/data/x.zip"),
                "RUNNING": _seq("pending", seq_uuid="u-running"),
            }
        },
    )
    data = create_sequences_data({"DONE", "RUNNING", "NEW"})
    assert data["DONE"]["status"] == "cached"
    assert data["DONE"]["zip_path"] == "/data/x.zip"
    assert data["RUNNING"]["status"] == "pending_external"
    assert data["RUNNING"]["seq_uuid"] == "u-running"
    assert data["NEW"]["status"] == "pending"
    assert data["NEW"]["seq_uuid"]


def test_pending_sequences_complete_when_zip_appears(volume):
    sequences = {
        "AAA": _seq("pending", seq_uuid="u1", sequence_id="a"),
        "CCC": _seq("pending_external", seq_uuid="u2", sequence_id="c"),
        "GGG": _seq("cached", zip_path="/old.zip", sequence_id="g"),
    }
    assert get_pending_uuids(sequences) == ["u1", "u2"]

    (get_zip_by_idr_dir() / "u1.zip").write_bytes(b"")
    completed = find_completed_in_cache(get_pending_uuids(sequences))
    assert list(completed) == ["u1.zip"]

    still_pending = update_sequences_with_completed(sequences, completed)
    assert still_pending == ["c"]
    assert sequences["AAA"]["status"] == "complete"
    assert sequences["AAA"]["seq_uuid"] is None
    assert sequences["AAA"]["zip_path"] == str(get_zip_by_idr_dir() / "u1.zip")
