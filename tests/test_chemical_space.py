import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from analysis import build_chemical_space, bulk_compare
from app import app


def test_build_chemical_space_returns_reference_and_coordinates():
    results, error = bulk_compare("CCO", ["CCN", "CCCl"], ["Amina", "Clorado"], "Morgan2", "Tanimoto")
    assert error is None
    assert "Ligações rotacionáveis (RotB)" in results[0]["properties"]
    space = build_chemical_space("CCO", ["CCN", "CCCl"], ["Amina", "Clorado"], results, "Morgan2", "Tanimoto")
    assert len(space["points"]) == 3
    assert any(point["role"] == "reference" for point in space["points"])
    assert all("x" in point and "y" in point for point in space["points"])
    assert all("physicochemical_distance" in point and "global_distance" in point for point in space["points"])
    assert "RotB" in space["descriptors"]
    assert space["weights"] == {"structural": 0.6, "physicochemical": 0.4}


def test_bulk_endpoint_conditionally_returns_chemical_space():
    client = app.test_client()
    response = client.post("/bulk-compare", json={
        "ref_smiles": "CCO",
        "smiles_list": ["CCN", "CCCl"],
        "names_list": ["Amina", "Clorado"],
        "fp_type": "Morgan2",
        "metric": "Tanimoto",
        "show_chemical_space": True,
    })
    assert response.status_code == 200, response.get_json()
    payload = response.get_json()
    assert len(payload["chemical_space"]["points"]) == 3

    legacy = client.post("/bulk-compare", json={
        "ref_smiles": "CCO",
        "smiles_list": ["CCN"],
        "names_list": ["Amina"],
        "fp_type": "Morgan2",
        "metric": "Tanimoto",
    })
    assert legacy.status_code == 200
    assert "chemical_space" not in legacy.get_json()


def test_batch_report_preview_includes_chemical_space_section():
    results, error = bulk_compare("CCO", ["CCN", "CCCl"], ["Amina", "Clorado"], "MACCS", "Tanimoto")
    assert error is None
    space = build_chemical_space("CCO", ["CCN", "CCCl"], ["Amina", "Clorado"], results, "MACCS", "Tanimoto")
    client = app.test_client()
    response = client.post("/report-preview", json={
        "mode": "batch",
        "ref_name": "Referência",
        "ref_smiles": "CCO",
        "results": results,
        "chemical_space": space,
        "fp_type": "MACCS",
        "metric": "Tanimoto",
    })
    assert response.status_code == 200
    html = response.get_data(as_text=True)
    assert "Espaço Químico 2D" in html
    assert "60% de distância estrutural + 40% de Dist.FQ normalizada" in html
    assert "RotB" in html
    assert "chemical-space-report-svg" in html


def test_build_chemical_space_displays_at_most_ten_nearest_neighbors():
    smiles_list = ["C" * length for length in range(2, 14)]
    names = [f"Mol_{index}" for index in range(len(smiles_list))]
    results, error = bulk_compare("CCO", smiles_list, names, "MACCS", "Tanimoto")
    assert error is None

    space = build_chemical_space("CCO", smiles_list, names, results, "MACCS", "Tanimoto")

    assert space["display_limit"] == 10
    assert space["total_valid_points"] == len(results) + 1
    assert len(space["points"]) == 11
    assert space["points"][0]["role"] == "reference"
    assert space["displayed_candidates"] == 10


def test_pair_endpoint_does_not_return_batch_chemical_space():
    client = app.test_client()
    response = client.post("/compare", json={
        "smiles_ref": "CCO",
        "smiles_test": "CCN",
        "name_ref": "Etanol",
        "name_test": "Etilamina",
        "fp_type": "MACCS",
        "metric": "Tanimoto",
        "show_similarity_map": False,
    })
    assert response.status_code == 200, response.get_json()
    payload = response.get_json()
    assert "chemical_space" not in payload
    assert "ema_space" not in payload
    assert "nitro_ra" not in payload


def test_batch_chemical_space_is_user_batch_only_contract():
    client = app.test_client()
    response = client.post("/bulk-compare", json={
        "ref_smiles": "CCO",
        "smiles_list": ["CCN", "CCCl"],
        "names_list": ["Amina", "Clorado"],
        "fp_type": "MACCS",
        "metric": "Tanimoto",
        "show_chemical_space": True,
    })
    assert response.status_code == 200, response.get_json()
    space = response.get_json()["chemical_space"]
    assert space["reference_included"] is True
    assert space["displayed_candidates"] == 2
    assert "PubChem" not in space.get("method", "")
    assert "EMA" not in space.get("method", "")
    assert "source" not in space or space.get("source") in (None, "")
