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
    assert cyp3a4["substrate"]["probability_percent"] == 0.1
    assert cyp3a4["inhibitor"]["prediction"] == "Non-Inhibitor"
    assert cyp3a4["inhibitor"]["probability"] == 0.064
    assert cyp3a4["inhibitor"]["probability_percent"] == 6.4


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


def test_main_analyze_route_includes_deep_pk_submission(monkeypatch):
    monkeypatch.setattr("app.evaluate_metabolism", lambda smiles: {"module": "metabolism", "status": "ok"})
    monkeypatch.setattr("app.submit_deep_pk_metabolism", lambda smiles: {"module": "deep_pk_metabolism", "status": "running", "job_id": "metabolism_1"})
    client = app.test_client()

    response = client.post("/nitro-ra/analyze", json={"smiles": "O=NN1CCCCC1", "modules": ["metabolism", "deep_pk"]})

    assert response.status_code == 200
    payload = response.get_json()
    assert payload["modules"] == ["metabolism", "deep_pk"]
    assert payload["results"]["deep_pk"]["job_id"] == "metabolism_1"


def test_main_analyze_route_requires_metabolism_for_deep_pk():
    client = app.test_client()

    response = client.post("/nitro-ra/analyze", json={"smiles": "O=NN1CCCCC1", "modules": ["deep_pk"]})

    assert response.status_code == 400
    assert "Metabolism" in str(response.get_json())


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


def test_deep_pk_report_layout_keeps_external_results_separate():
    payload = {
        "module": "nitro_ra",
        "mode": "nitro",
        "status": "ok",
        "smiles": "O=NN1CCCCC1",
        "name": "N-nitroso-piperidina",
        "modules": ["metabolism", "deep_pk"],
        "results": {
            "metabolism": {"status": "ok", "summary": {"alpha_sites": 0, "metabolites": 0, "reactive_intermediates": 0}},
            "deep_pk": {
                "module": "deep_pk_metabolism",
                "provider": "Deep-PK",
                "prediction_type": "metabolism",
                "status": "ok",
                "message": "Endpoints recuperados.",
                "smiles": "O=NN1CCCCC1",
                "job_id": "metabolism_123",
                "isoforms": [{
                    "isoform": "CYP3A4",
                    "substrate": {"prediction": "Non-Substrate", "probability": 0.001, "interpretation": "Alta confiança"},
                    "inhibitor": {"prediction": "Non-Inhibitor", "probability": 0.064, "interpretation": "Alta confiança"},
                }],
                "warnings": [],
            },
        },
    }
    response = app.test_client().post("/report-preview", json=payload)

    assert response.status_code == 200
    html = response.get_data(as_text=True)
    deep_pk_html = html.split('id="report-deep-pk"', 1)[1].split('id="report-', 1)[0]
    assert "nitro-deep-pk-summary" in deep_pk_html
    assert "Escopo da análise" in deep_pk_html
    assert "Identificação da consulta" in deep_pk_html
    assert "Isoforma CYP" in deep_pk_html
    assert "Substrato" in deep_pk_html
    assert "Inibidor" in deep_pk_html
    assert "metabolism_123" in deep_pk_html
    assert "complemento externo" in deep_pk_html
    assert 'class="nitro-deep-pk-section"' in html
    assert html.index("Resultado Deep-PK") < html.index("Ativação Enzimática")
    assert html.index("Resultado Deep-PK") < html.index("Regra e contexto enzimático")


def test_deep_pk_report_layout_handles_running_state_without_table():
    payload = {
        "module": "nitro_ra",
        "mode": "nitro",
        "status": "ok",
        "smiles": "O=NN1CCCCC1",
        "modules": ["metabolism", "deep_pk"],
        "results": {
            "metabolism": {"status": "ok", "summary": {"alpha_sites": 0, "metabolites": 0, "reactive_intermediates": 0}, "alpha_sites": []},
            "deep_pk": {
                "module": "deep_pk_metabolism",
                "provider": "Deep-PK",
                "status": "running",
                "message": "Consulta enviada; aguardando os endpoints CYP.",
                "job_id": "metabolism_456",
                "isoforms": [],
                "warnings": [],
            },
        },
    }
    response = app.test_client().post("/report-preview", json=payload)

    assert response.status_code == 200
    deep_pk_html = response.get_data(as_text=True).split('id="report-deep-pk"', 1)[1].split('id="report-', 1)[0]
    assert "Consulta enviada; aguardando os endpoints CYP." in deep_pk_html
    assert "nitro-deep-pk-table" not in deep_pk_html
    assert "metabolism_456" in deep_pk_html


def test_success_presentation_uses_result_available_without_completed_label():
    presentation = deep_pk.normalize_deep_pk_message({"status": "ok"})["presentation"]
    assert presentation["label"] == "Resultados disponíveis"
    assert presentation["label"] != "Concluído"


def test_presentation_normalizer_maps_states_to_actionable_messages():
    expected = {
        "running": ("Processando", "info"),
        "deep_pk_unavailable": ("Indisponível", "error"),
        "deep_pk_error": ("Resposta não processada", "error"),
        "deep_pk_timeout": ("Tempo excedido", "warning"),
        "invalid_smiles": ("SMILES inválido", "error"),
        "invalid_job_id": ("Consulta inválida", "error"),
    }
    for status, (label, tone) in expected.items():
        presentation = deep_pk.normalize_deep_pk_message({"status": status})["presentation"]
        assert presentation["label"] == label
        assert presentation["tone"] == tone
        assert presentation["message"]
        assert presentation["action"]


def test_base_result_contains_presentation_contract():
    result = deep_pk._base_result("O=NN1CCCCC1", "deep_pk_timeout", "mensagem técnica")
    assert result["presentation"]["status"] == "deep_pk_timeout"
    assert result["presentation"]["label"] == "Tempo excedido"
    assert "mensagem técnica" == result["message"]
