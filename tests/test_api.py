"""The REST API the Colab notebook calls, run offline.

`fastapi_app.local()` builds the FastAPI app without Modal credentials. The
volume is a temp dir, `vol.commit()` is a no-op and `process_16_sequences.spawn`
records batches instead of starting workers.
"""

import io
import uuid
import warnings
import zipfile

import pytest
from conftest import write_nardini_zip
from fastapi.testclient import TestClient

import app.backend as backend
from shared_utils.file_utils import get_run_metadata, get_zip_by_idr_dir

FASTA = (
    ">Drs1_IDR_1_[0-10]\nMVVGTKKYSN\n"
    ">Drs1_IDR_2_[660-680]\nMLKHKKEIQARPRRTWFQSE\n"
    ">Drs1_IDR_3_[0-5]\nMKVLL\n"
)


class FakeVolume:
    def commit(self):
        pass

    def reload(self):
        pass


class FakeCall:
    def __init__(self, n):
        self.object_id = f"fc-{n}"


class FakeWorker:
    def __init__(self):
        self.batches = []

    def spawn(self, batch):
        self.batches.append(batch)
        return FakeCall(len(self.batches))

    def finish_all(self):
        """Write a result zip for every spawned sequence, like the worker does."""
        for batch in self.batches:
            for si in batch:
                rec = si["sequence"]
                write_nardini_zip(
                    get_zip_by_idr_dir() / f"{si['seq_uuid']}.zip",
                    rec.id,
                    str(rec.seq),
                )


@pytest.fixture
def worker(volume, monkeypatch):
    worker = FakeWorker()
    monkeypatch.setattr(backend, "vol", FakeVolume())
    monkeypatch.setattr(backend, "process_16_sequences", worker)
    return worker


@pytest.fixture
def client(worker):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        api = backend.fastapi_app.local()
    return TestClient(api)


def upload(client, text=FASTA, name="yeast.fasta", **params):
    return client.post(
        "/upload_fasta", files={"file": (name, text.encode())}, params=params
    )


def test_health(client):
    assert client.get("/health").json() == {"status": "healthy"}


def test_upload_status_download(client, worker):
    res = upload(client)
    assert res.status_code == 200
    body = res.json()
    run_id = body["run_id"]
    uuid.UUID(run_id)
    assert body["status"] == "submitted"
    assert body["job_ids"] == ["fc-1"]
    assert body["status_url"] == f"/status/{run_id}"
    assert body["download_url"] is None

    # Longest sequences are spawned first.
    assert [si["sequence"].id for si in worker.batches[0]] == [
        "Drs1_IDR_2_[660-680]",
        "Drs1_IDR_1_[0-10]",
        "Drs1_IDR_3_[0-5]",
    ]

    status = client.get(f"/status/{run_id}").json()
    assert status["status"] == "pending"
    assert sorted(status["pending_sequences"]) == [
        "Drs1_IDR_1_[0-10]",
        "Drs1_IDR_2_[660-680]",
        "Drs1_IDR_3_[0-5]",
    ]
    assert client.get(f"/download/{run_id}").status_code == 404

    worker.finish_all()
    assert client.get(f"/status/{run_id}").json() == {
        "run_id": run_id,
        "status": "complete",
        "pending_sequences": [],
        "download_url": f"/download/{run_id}",
    }

    res = client.get(f"/download/{run_id}")
    assert res.status_code == 200
    assert res.headers["content-type"] == "application/zip"
    assert (
        res.headers["content-disposition"] == "attachment; filename=yeast_results.zip"
    )
    with zipfile.ZipFile(io.BytesIO(res.content)) as zf:
        names = zf.namelist()
        master = zf.read("master_sequences.tsv").decode().splitlines()
    assert "regular-Drs1_IDR_3_[0-5].png" in names
    assert len(master) == 4

    res = client.get(f"/download/{run_id}", params={"output_filename": "x.fasta"})
    assert res.headers["content-disposition"] == "attachment; filename=x.zip"


def test_reupload_is_served_from_cache(client, worker):
    first = upload(client).json()
    worker.finish_all()
    client.get(f"/status/{first['run_id']}")

    second = upload(client, name="again.fa").json()
    assert second["status"] == "ready"
    assert second["job_ids"] == []
    assert second["download_url"] == f"/download/{second['run_id']}"
    assert len(worker.batches) == 1
    assert client.get(f"/status/{second['run_id']}").json()["status"] == "complete"
    res = client.get(f"/download/{second['run_id']}")
    assert res.status_code == 200
    assert (
        res.headers["content-disposition"] == "attachment; filename=again_results.zip"
    )


def test_batches_of_16(client, worker):
    text = "".join(f">s{i}\n{'A' * (i + 1)}\n" for i in range(40))
    body = upload(client, text=text).json()
    assert [len(b) for b in worker.batches] == [16, 16, 8]
    assert body["job_ids"] == ["fc-1", "fc-2", "fc-3"]


def test_output_filename_param(client, worker):
    run_id = upload(client, output_filename="custom.fasta").json()["run_id"]
    worker.finish_all()
    res = client.get(f"/download/{run_id}")
    assert res.headers["content-disposition"] == "attachment; filename=custom.zip"


@pytest.mark.parametrize(
    "name, text, code",
    [
        ("notes.txt", FASTA, 400),
        ("empty.fasta", "", 400),
        ("big.fasta", ">a\n" + "A" * (backend.MAX_FILE_SIZE + 1), 413),
    ],
)
def test_upload_rejects(client, name, text, code):
    assert upload(client, text=text, name=name).status_code == code


def test_bad_run_ids(client):
    assert client.get("/status/not-a-uuid").status_code == 400
    assert client.get("/download/not-a-uuid").status_code == 400
    assert client.get(f"/status/{uuid.uuid4()}").status_code == 404
    assert client.get(f"/download/{uuid.uuid4()}").status_code == 404


class FakeSingleWorker:
    def __init__(self):
        self.inputs = []

    def spawn(self, sequence_input):
        self.inputs.append(sequence_input)
        return FakeCall(f"retry-{len(self.inputs)}")


def test_retry_spawns_the_retry_function(client, monkeypatch):
    retry = FakeSingleWorker()
    monkeypatch.setattr(backend, "retry_pending_sequences", retry)
    run_id = upload(client).json()["run_id"]

    assert client.get(f"/retry/{run_id}").json() == {
        "run_id": run_id,
        "status": "retry_submitted",
    }
    assert retry.inputs == [run_id]
    assert client.get("/retry/not-a-uuid").status_code == 400
    assert client.get(f"/retry/{uuid.uuid4()}").status_code == 404
    assert retry.inputs == [run_id]


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_retry_respawns_only_unfinished_sequences(client, worker, monkeypatch):
    single = FakeSingleWorker()
    monkeypatch.setattr(backend, "process_single_sequence", single)
    run_id = upload(client).json()["run_id"]
    done = worker.batches[0][0]
    write_nardini_zip(
        get_zip_by_idr_dir() / f"{done['seq_uuid']}.zip",
        done["sequence"].id,
        str(done["sequence"].seq),
    )

    job_ids = backend.retry_pending_sequences.local(run_id)

    assert sorted(job_ids) == ["fc-retry-1", "fc-retry-2"]
    retried = {str(si["sequence"].seq): si for si in single.inputs}
    all_seqs = {"MVVGTKKYSN", "MLKHKKEIQARPRRTWFQSE", "MKVLL"}
    assert set(retried) == all_seqs - {str(done["sequence"].seq)}
    sequences = get_run_metadata(run_id)["sequences"]
    assert sequences[str(done["sequence"].seq)]["status"] == "complete"
    for seq, si in retried.items():
        assert si["sequence"].id == sequences[seq]["sequence_id"]
        assert si["seq_uuid"] == sequences[seq]["seq_uuid"]
        assert sequences[seq]["job_id"] in job_ids

    worker.finish_all()
    assert backend.retry_pending_sequences.local(run_id) == []
