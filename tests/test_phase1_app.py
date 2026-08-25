import base64
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app import app


VALID_PNG_BASE64 = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO5mY8kAAAAASUVORK5CYII="


@pytest.fixture()
def client():
    app.config["TESTING"] = True
    with app.test_client() as test_client:
        yield test_client


def _assert_invalid_request_contract(response):
    body = response.get_json()
    assert isinstance(body, dict)
    assert "error" in body
    assert body["error"]["code"] == "invalid_request"
    assert "message" in body["error"]


@pytest.mark.parametrize("route", ["/compare", "/bulk-compare", "/report-preview", "/export-pdf"])
def test_requires_application_json(route, client):
    response = client.post(route, data="{}", content_type="text/plain")
    assert response.status_code == 415
    _assert_invalid_request_contract(response)


@pytest.mark.parametrize("route", ["/compare", "/bulk-compare", "/report-preview", "/export-pdf"])
def test_invalid_json_never_returns_500(route, client):
    response = client.post(route, data="{", content_type="application/json")
    assert response.status_code == 400
    _assert_invalid_request_contract(response)


def test_compare_rejects_non_object_payload(client):
    response = client.post("/compare", json=["not-an-object"])
    assert response.status_code == 400
    _assert_invalid_request_contract(response)


def test_compare_rejects_wrong_smiles_type(client):
    response = client.post(
        "/compare",
        json={"smiles_ref": 123, "smiles_test": "CCO", "metric": "Tanimoto", "fp_type": "Morgan2"},
    )
    assert response.status_code == 400
    _assert_invalid_request_contract(response)


def test_compare_accepts_legacy_show_logd_bool(client):
    response = client.post(
        "/compare",
        json={"smiles_ref": "CCO", "smiles_test": "CCN", "metric": "Tanimoto", "fp_type": "Morgan2", "show_logd": True},
    )
    assert response.status_code != 500


def test_lookup_name_rejects_non_string_smiles(client):
    response = client.post("/lookup-name", json={"smiles": ["CCO"]})
    assert response.status_code == 400
    _assert_invalid_request_contract(response)


def test_lookup_name_returns_empty_name_for_blank_smiles(client):
    response = client.post("/lookup-name", json={"smiles": "   "})
    assert response.status_code == 200
    assert response.get_json() == {"name": None}


def test_compare_rejects_invalid_metric(client):
    response = client.post(
        "/compare",
        json={"smiles_ref": "CCO", "smiles_test": "CCN", "metric": "Cosine", "fp_type": "Morgan2"},
    )
    assert response.status_code == 400
    _assert_invalid_request_contract(response)


def test_bulk_compare_names_length_must_match(client):
    response = client.post(
        "/bulk-compare",
        json={
            "ref_smiles": "CCO",
            "smiles_list": ["CCN", "CCC"],
            "names_list": ["ApenasUm"],
            "metric": "Tanimoto",
            "fp_type": "Morgan2",
        },
    )
    assert response.status_code == 400
    _assert_invalid_request_contract(response)


def test_export_pdf_rejects_non_finite_similarity(client):
    response = client.post("/export-pdf", json={"similarity": float("inf")})
    assert response.status_code == 400
    _assert_invalid_request_contract(response)


def test_export_pdf_rejects_invalid_base64(client):
    response = client.post("/export-pdf", json={"png_ref": "%%%%"})
    assert response.status_code == 400
    _assert_invalid_request_contract(response)


def test_export_pdf_rejects_non_png_base64(client):
    non_png = base64.b64encode(b"not-a-png").decode("ascii")
    response = client.post("/export-pdf", json={"png_ref": non_png})
    assert response.status_code == 400
    _assert_invalid_request_contract(response)


def test_export_pdf_accepts_png_data_uri(client):
    response = client.post("/export-pdf", json={"png_ref": f"data:image/png;base64,{VALID_PNG_BASE64}"})
    assert response.status_code == 200


def test_healthz_is_alive(client):
    response = client.get("/healthz")
    assert response.status_code == 200
    assert response.get_json() == {"status": "ok"}


def test_healthz_accepts_head(client):
    response = client.head("/healthz")
    assert response.status_code == 200
    assert response.data == b""


def test_compare_produces_valid_png_fingerprints():
    from analysis import compare

    result, error = compare(
        "CCO",
        "CCN",
        "Referência",
        "Teste",
        "Morgan2",
        "Tanimoto",
        show_map=False,
    )
    assert error is None
    assert result["fingerprint_ref_png"]
    assert result["fingerprint_test_png"]
    ref_png = base64.b64decode(result["fingerprint_ref_png"])
    test_png = base64.b64decode(result["fingerprint_test_png"])
    assert ref_png.startswith(b"\x89PNG\r\n\x1a\n")
    assert test_png.startswith(b"\x89PNG\r\n\x1a\n")
    assert all("Ligações Rotacionais" not in prop["Propriedade"] for prop in result["properties"])


def test_method_not_allowed_has_structured_error(client):
    response = client.get("/compare")
    assert response.status_code == 405
    body = response.get_json()
    assert body["error"]["code"] == "method_not_allowed"
    assert "message" in body["error"]


def test_cors_default_allows_127_localhost(client):
    response = client.post(
        "/report-preview",
        json={"a": 1},
        headers={"Origin": "http://127.0.0.1:3000"},
    )
    assert response.status_code == 200
    assert response.headers.get("Access-Control-Allow-Origin") == "http://127.0.0.1:3000"


def test_cors_default_allows_localhost(client):
    response = client.post(
        "/report-preview",
        json={"a": 1},
        headers={"Origin": "http://localhost:3000"},
    )
    assert response.status_code == 200
    assert response.headers.get("Access-Control-Allow-Origin") == "http://localhost:3000"


def test_cors_default_blocks_other_origins(client):
    response = client.post(
        "/report-preview",
        json={"a": 1},
        headers={"Origin": "http://evil.example"},
    )
    assert response.status_code == 200
    assert response.headers.get("Access-Control-Allow-Origin") is None


def test_security_headers_are_present(client):
    response = client.get("/healthz")
    assert response.status_code == 200
    assert response.headers.get("X-Content-Type-Options") == "nosniff"
    assert response.headers.get("X-Frame-Options") == "DENY"
    assert response.headers.get("Referrer-Policy") == "strict-origin-when-cross-origin"


def test_request_too_large_has_structured_error(client):
    old_limit = app.config.get("MAX_CONTENT_LENGTH")
    app.config["MAX_CONTENT_LENGTH"] = 5
    try:
        response = client.post("/compare", data="123456", content_type="application/json")
        assert response.status_code == 413
        _assert_invalid_request_contract(response)
    finally:
        app.config["MAX_CONTENT_LENGTH"] = old_limit


def test_pubchem_lookup_uses_cache_and_handles_failures(monkeypatch):
    import app as app_module

    app_module.PUBCHEM_CACHE.clear()

    class DummyResponse:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self):
            return b'{"PropertyTable": {"Properties": [{"Title": "Ethanol"}]}}'

    calls = {"count": 0}

    def fake_urlopen(req, timeout=None):
        calls["count"] += 1
        return DummyResponse()

    monkeypatch.setattr(app_module.urllib_request, "urlopen", fake_urlopen)
    assert app_module.get_pubchem_name_for_smiles("CCO") == "Ethanol"
    assert app_module.get_pubchem_name_for_smiles("CCO") == "Ethanol"
    assert calls["count"] == 1

    def fake_urlopen_failure(req, timeout=None):
        raise ConnectionError("boom")

    monkeypatch.setattr(app_module.urllib_request, "urlopen", fake_urlopen_failure)
    assert app_module.get_pubchem_name_for_smiles("CCN") is None
