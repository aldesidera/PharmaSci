import pytest

from app import app
from chemo_suite.apps.nitro_ra.cpca import calculate_cpca, evaluate_cpca


@pytest.mark.parametrize(
    ("smiles", "category", "ai_ng_day"),
    [
        ("O=NN1CCCC1", 4, 1500.0),  # pyrrolidine: alpha-H 2,2 + ring +3
        ("O=NN1CCCCC1", 3, 400.0),  # piperidine: alpha-H 2,2 + ring +2
        ("O=NN1CCOCC1", 2, 100.0),  # morpholine: alpha-H 2,2 + ring +1
    ],
)
def test_cpca_supported_ring_examples(smiles, category, ai_ng_day):
    result = evaluate_cpca(smiles)

    assert result["status"] == "ok"
    assert result["potency_category"] == category
    assert result["ai_ng_day"] == ai_ng_day
    assert result["center_count"] == 1
    assert result["centers"][0]["alpha_hydrogen_counts"] == [2, 2]


def test_cpca_returns_manual_review_for_unmapped_alpha_hydrogen_pair():
    result = evaluate_cpca("CN(C)N=O")

    assert result["status"] == "manual_review"
    assert result["center_count"] == 1
    assert result["centers"][0]["alpha_hydrogen_counts"] == [3, 3]
    assert "não está definido" in result["centers"][0]["message"]


def test_cpca_rejects_non_nitrosamine_and_invalid_smiles():
    assert evaluate_cpca("CCO")["status"] == "not_nitrosamine"
    assert evaluate_cpca("not-a-smiles")["status"] == "invalid_smiles"


def test_cpca_excludes_nitrosamide_like_structure():
    result = evaluate_cpca("CC(=O)N(C)N=O")

    assert result["status"] == "not_applicable"
    assert result["center_count"] == 1
    assert result["centers"][0]["excluded_reasons"]


def test_cpca_converts_ai_to_ppm_using_mdd():
    result = calculate_cpca("O=NN1CCCCC1", mdd_mg=10)

    assert result["ai_ng_day"] == 400.0
    assert result["mdd_mg"] == 10.0
    assert result["ppm_limit"] == 40.0
    assert result["ppm_formula"] == "AI (ng/dia) / dose diária máxima (mg)"


def test_cpca_api_returns_structured_result():
    client = app.test_client()
    response = client.post(
        "/nitro-ra/cpca",
        json={"smiles": "O=NN1CCCCC1", "mdd_mg": 10},
    )

    assert response.status_code == 200
    payload = response.get_json()
    assert payload["status"] == "ok"
    assert payload["potency_category"] == 3
    assert payload["ppm_limit"] == 40.0


def test_cpca_api_requires_json_content_type():
    client = app.test_client()
    response = client.post("/nitro-ra/cpca", data='{"smiles":"CCO"}')

    assert response.status_code == 415
    assert response.get_json()["error"]["code"] == "invalid_request"


def test_cpca_api_reports_invalid_smiles_as_bad_request():
    client = app.test_client()
    response = client.post("/nitro-ra/cpca", json={"smiles": "not-a-smiles"})

    assert response.status_code == 400
    assert response.get_json()["status"] == "invalid_smiles"


def test_nitro_ra_analyze_returns_one_result_per_selected_module():
    client = app.test_client()
    response = client.post(
        "/nitro-ra/analyze",
        json={
            "smiles": "O=NN1CCCCC1",
            "mdd_mg": 10,
            "modules": ["cpca", "quantum", "metabolism"],
        },
    )

    assert response.status_code == 200
    payload = response.get_json()
    assert payload["module"] == "nitro_ra"
    assert payload["modules"] == ["cpca", "quantum", "metabolism"]
    assert payload["results"]["cpca"]["potency_category"] == 3
    assert payload["results"]["quantum"]["status"] == "not_implemented"
    assert payload["results"]["metabolism"]["status"] == "ok"
    assert payload["results"]["metabolism"]["summary"]["alpha_sites"] >= 1
    assert payload["results"]["metabolism"]["summary"]["metabolites"] >= 1


def test_nitro_ra_analyze_rejects_missing_or_invalid_modules():
    client = app.test_client()
    missing = client.post("/nitro-ra/analyze", json={"smiles": "O=NN1CCCCC1"})
    invalid = client.post(
        "/nitro-ra/analyze",
        json={"smiles": "O=NN1CCCCC1", "modules": ["unknown"]},
    )

    assert missing.status_code == 400
    assert invalid.status_code == 400


def test_cpca_includes_structure_svg_and_ema_appendix_match():
    result = evaluate_cpca("O=NN1CCCCC1", mdd_mg=10)

    assert result["structure_svg"].startswith("<?xml")
    assert result["canonical_smiles"] == "O=NN1CCCCC1"
    assert result["ema"]["listed"] is True
    assert result["ema"]["status"] == "listed"
    assert result["ema"]["name"] == "N-nitroso-piperidine"
    assert result["ema"]["ai_ng_day"] == 1300
    assert result["ema"]["ppm_limit"] == 130.0
    assert result["ema"]["reference_number"] == "EMA/42261/2025 Rev.13"


def test_cpca_reports_when_structure_is_not_in_ema_appendix():
    result = evaluate_cpca("O=NN1CCCCC1C")

    assert result["ema"]["listed"] is False
    assert result["ema"]["status"] == "not_listed"
    assert result["ema"].get("ai_ng_day") is None
    assert result["ema"]["message"].startswith("Nitrosamina não listada no Apêndice I da EMA")


def test_nitro_ra_analyze_dispatches_nitrosamine_space(monkeypatch):
    captured = {}

    def fake_space(smiles):
        captured["smiles"] = smiles
        return {"module": "nitrosamine_space", "status": "no_nitrosamines", "candidates": [], "points": []}

    monkeypatch.setattr("app.search_nitrosamine_space", fake_space)
    client = app.test_client()
    response = client.post(
        "/nitro-ra/analyze",
        json={"smiles": "O=NN1CCCCC1", "modules": ["nitrosamine_space"]},
    )

    assert response.status_code == 200
    assert captured["smiles"] == "O=NN1CCCCC1"
    payload = response.get_json()
    assert payload["modules"] == ["nitrosamine_space"]
    assert payload["results"]["nitrosamine_space"]["status"] == "no_nitrosamines"


def test_nitro_ra_analyze_preserves_pubchem_unavailable_state(monkeypatch):
    def fake_space(smiles):
        return {
            "module": "nitrosamine_space",
            "status": "pubchem_unavailable",
            "message": "PubChem retornou HTTP 503.",
            "search": {"retrieved_cids": 0, "n_nitroso_candidates": 0, "selected_candidates": 0},
            "candidates": [],
            "points": [],
        }

    monkeypatch.setattr("app.search_nitrosamine_space", fake_space)
    client = app.test_client()
    response = client.post(
        "/nitro-ra/analyze",
        json={"smiles": "O=NN1CCCCC1", "modules": ["cpca", "nitrosamine_space"]},
    )

    assert response.status_code == 200
    payload = response.get_json()
    space_result = payload["results"]["nitrosamine_space"]
    assert space_result["status"] == "pubchem_unavailable"
    assert space_result["search"]["retrieved_cids"] == 0
    assert payload["results"]["cpca"]["status"] == "ok"
