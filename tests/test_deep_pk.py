import json

from app import app
from chemo_suite.apps.nitro_ra import deep_pk


def test_submit_uses_canonical_smiles_and_returns_job(monkeypatch):
    captured = {}

    def fake_request(method, *, form=None, query=None):
        captured.update(method=method, form=form, query=query)
        return {"job_id": "metabolism_123.45"}

    monkeypatch.setattr(deep_pk, "_request_json", fake_request)
    result = deep_pk.submit_deep_pk_metabolism("C1CCCCC1")

    assert result["status"] == "running"
    assert result["job_id"] == "metabolism_123.45"
    assert captured["method"] == "POST"
    assert captured["form"] == {"smiles": "C1CCCCC1", "pred_type": "metabolism"}


def test_completed_job_normalizes_cyp_substrate_and_inhibitor_endpoints(monkeypatch):
    payload = {
        "0": {
            "SMILES": "O=NN1CCCCC1",
            "[Metabolism/CYP 3A4 Substrate] Predictions": "Non-Substrate",
            "[Metabolism/CYP 3A4 Substrate] Probability": 0.001,
            "[Metabolism/CYP 3A4 Substrate] Interpretation": "Non-Substrate (High Confidence)",
            "[Metabolism/CYP 3A4 Inhibitor] Predictions": "Non-Inhibitor",
            "[Metabolism/CYP 3A4 Inhibitor] Probability": 0.064,
            "[Metabolism/CYP 3A4 Inhibitor] Interpretation": "Non-Inhibitor (High Confidence)",
        }
    }
    monkeypatch.setattr(deep_pk, "_request_json", lambda method, **kwargs: payload)

    result = deep_pk.get_deep_pk_metabolism("metabolism_123.45")
    cyp3a4 = next(item for item in result["isoforms"] if item["isoform"] == "CYP3A4")

    assert result["status"] == "ok"
    assert result["smiles"] == "O=NN1CCCCC1"
    assert cyp3a4["substrate"]["prediction"] == "Non-Substrate"
    assert cyp3a4["substrate"]["probability"] == 0.001
    assert cyp3a4["inhibitor"]["prediction"] == "Non-Inhibitor"
    assert cyp3a4["inhibitor"]["probability"] == 0.064


def test_running_job_is_preserved(monkeypatch):
    monkeypatch.setattr(deep_pk, "_request_json", lambda method, **kwargs: {"status": "running"})

    result = deep_pk.get_deep_pk_metabolism("metabolism_123.45")

    assert result["status"] == "running"
    assert result["job_id"] == "metabolism_123.45"
    assert result["isoforms"] == []


def test_network_failure_is_explicit_and_non_fatal(monkeypatch):
    def fail_request(method, **kwargs):
        raise deep_pk.DeepPkError("serviço indisponível", status_code=503)

    monkeypatch.setattr(deep_pk, "_request_json", fail_request)
    result = deep_pk.submit_deep_pk_metabolism("O=NN1CCCCC1")

    assert result["status"] == "deep_pk_unavailable"
    assert result["http_status"] == 503
    assert result["isoforms"] == []


def test_invalid_job_id_is_rejected_without_network_call():
    result = deep_pk.get_deep_pk_metabolism("../../etc/passwd")

    assert result["status"] == "invalid_job_id"
    assert result["isoforms"] == []


def test_submit_route_returns_job_payload(monkeypatch):
    monkeypatch.setattr(
        "app.submit_deep_pk_metabolism",
        lambda smiles: {"module": "deep_pk_metabolism", "status": "running", "job_id": "metabolism_1"},
    )
    client = app.test_client()

    response = client.post("/nitro-ra/deep-pk", json={"smiles": "O=NN1CCCCC1"})

    assert response.status_code == 200
    assert response.get_json()["job_id"] == "metabolism_1"


def test_status_route_returns_normalized_payload(monkeypatch):
    expected = {"module": "deep_pk_metabolism", "status": "ok", "isoforms": []}
    monkeypatch.setattr("app.get_deep_pk_metabolism", lambda job_id: expected)
    client = app.test_client()

    response = client.get("/nitro-ra/deep-pk/metabolism_1")

    assert response.status_code == 200
    assert response.get_json() == expected


def test_submit_route_rejects_invalid_smiles(monkeypatch):
    monkeypatch.setattr(
        "app.submit_deep_pk_metabolism",
        lambda smiles: {"module": "deep_pk_metabolism", "status": "invalid_smiles", "isoforms": []},
    )
    client = app.test_client()

    response = client.post("/nitro-ra/deep-pk", json={"smiles": "not-a-smiles"})

    assert response.status_code == 400
    assert response.get_json()["status"] == "invalid_smiles"


def test_multipart_builder_contains_form_fields():
    body, content_type = deep_pk._multipart_form({"smiles": "O=NN1CCCCC1", "pred_type": "metabolism"})

    assert content_type.startswith("multipart/form-data; boundary=")
    assert b'name="smiles"' in body
    assert b"O=NN1CCCCC1" in body
    assert b'name="pred_type"' in body
