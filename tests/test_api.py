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
from shared_utils.file_utils import get_zip_by_idr_dir

FASTA = (
    ">Drs1_IDR_1_[0-10]\nMVVGTKKYSN\n"
    ">Drs1_IDR_2_[660-680]\nMLKHKKEIQARPRRTWFQSE\n"
    ">Drs1_IDR_3_[0-5]\nMKVLL\n"
)


class FakeVolume:
    def commit(self):
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
    assert set(body) == {"run_id", "status", "message", "job_ids"}
    assert body["status"] == "submitted"
    assert body["job_ids"] == ["fc-1"]
    run_id = body["run_id"]
    uuid.UUID(run_id)

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


def test_retry_acknowledges(client):
    # Current behavior: /retry only acknowledges; it does not respawn jobs.
    run_id = str(uuid.uuid4())
    assert client.get(f"/retry/{run_id}").json() == {
        "run_id": run_id,
        "status": "retry_submitted",
    }
