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
