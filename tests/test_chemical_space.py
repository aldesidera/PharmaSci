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
    assert all("structural_distance" in point and "physicochemical_distance" in point and "global_distance" in point for point in space["points"])
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
    assert payload.get("reference_png")

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
        "reference_png": "cGxhY2Vob2xkZXI=",
        "chemical_space": space,
        "fp_type": "MACCS",
        "metric": "Tanimoto",
    })
    assert response.status_code == 200
    html = response.get_data(as_text=True)
    assert "Espaço Químico 2D" in html
    assert "reference-card--with-structure" in html
    assert "chemical-space-report-layout" in html
    assert "chemical-space-report-table--side" in html
    assert ".batch-similarity-table thead { display: table-header-group; }" in html
    assert "body.batch-report .batch-similarity-table { break-inside: auto; }" in html
    assert 'text-anchor="start"' in html
    assert 'style="fill:#1e293b;"' in html
    assert html.index("Similaridade e Propriedades Físico-Químicas") < html.index("Espaço Químico 2D")
    assert "Similaridade estrutural" in html
    assert "Propriedades físico-químicas" in html
    assert "Descritores:</strong> Massa Molecular (MW), LogP, TPSA, HBD, HBA e RotB." in html
    assert "Distância de similaridade" in html
    assert "Distância físico-química" in html
    assert "Dist.FQ normalizada" not in html
    assert "chemical-space-report-svg--batch" in html


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



def test_nitro_report_preview_is_modular_and_preserves_quantum_state():
    client = app.test_client()
    response = client.post("/report-preview", json={
        "mode": "nitro",
        "module": "nitro_ra",
        "status": "ok",
        "smiles": "O=NN1CCCCC1",
        "modules": ["quantum"],
        "results": {
            "quantum": {
                "module": "quantum",
                "status": "not_implemented",
                "message": "Cálculo de HOMO/LUMO será integrado na próxima etapa do Nitro.RA.",
                "smiles": "O=NN1CCCCC1",
            }
        },
    })
    assert response.status_code == 200
    html = response.get_data(as_text=True)
    assert "Relatório de Análise Nitro.RA" in html
    assert "Quantum" in html
    assert "O modelo Quantum ainda não foi implementado nesta versão do aplicativo." in html
    assert "Cálculo de HOMO/LUMO será integrado" not in html
    assert "Escopo atual" not in html
    assert "HOMO" not in html
    assert "LUMO" not in html
    assert "cPCA" not in html
    assert "Metabolism" not in html
    assert "PubChem" not in html
    assert "1300 ng/dia" not in html



def test_nitro_space_report_page_two_uses_compact_vertical_layout():
    target = {"name": "Molécula alvo", "is_target": True, "x": 0.0, "y": 0.0}
    candidate = {
        "name": "N-nitrosopiperidina",
        "source": "PubChem",
        "global_distance": 0.123456,
        "similarity": 0.876543,
    }
    space = {
        "module": "nitrosamine_space",
        "status": "ok",
        "points": [target, {"name": candidate["name"], "x": 0.2, "y": -0.1, "is_target": False}],
        "candidates": [candidate],
        "search": {
            "scored_candidates": 40,
            "selected_candidates": 1,
            "mds_stress": 0.123456,
            "selection_method": "texto de método que não deve aparecer no relatório",
        },
        "warnings": [
            "A busca depende da disponibilidade e da cobertura do PubChem; ausência no lote não prova ausência no universo químico.",
            "A similaridade 2D, os descritores e a projeção MDS são ferramentas de triagem.",
        ],
        "ema_space": {
            "points": [target],
            "candidates": [],
            "search": {"scored_candidates": 0, "selected_candidates": 0, "mds_stress": 0.0},
            "warnings": ["Aviso interno que não deve aparecer abaixo da tabela EMA"],
        },
    }
    response = app.test_client().post("/report-preview", json={
        "mode": "nitro",
        "module": "nitro_ra",
        "status": "ok",
        "smiles": "O=NN1CCCCC1",
        "name": "N-nitroso-piperidina",
        "modules": ["nitrosamine_space"],
        "results": {"nitrosamine_space": space},
    })
    assert response.status_code == 200
    html = response.get_data(as_text=True)
    assert "PubChem - Lote de 40 CIDs" in html
    assert "N-nitroso-piperidina" in html
    assert "Molécula alvo" not in html
    assert "Estruturas Avaliadas" in html
    assert "Estruturas Consideradas" in html
    assert "space-meta div { display: flex; min-height: 58px; flex-direction: column; align-items: center; justify-content: center; text-align: center; }" in html
    assert "space-meta span, .nitro-space-pair .space-meta strong { width: 100%; text-align: center; }" in html
    assert "space-metric-help { display: block; width: 100%;" in html
    assert "0.1235" in html
    assert "Diferença entre as distâncias originais e a projeção 2D" in html
    assert "texto de método que não deve aparecer no relatório" not in html
    assert "Método:" not in html
    assert "A busca depende da disponibilidade e da cobertura do PubChem" not in html
    assert "A similaridade 2D, os descritores e a projeção MDS" not in html
    assert "Aviso interno que não deve aparecer abaixo da tabela EMA" not in html
    assert "grid-template-columns: minmax(0, 1fr); gap: 8px;" in html
    assert "height: 340px" in html
    assert "MDS 1" in html
    assert "MDS 2" in html
    assert 'text-anchor="' in html
    assert 'font-size="12" fill="#1e293b"' in html
    assert "Distância de similaridade" in html
    assert "Distância físico-química" in html
    assert "Distância global" in html
    assert "AI / metadado" not in html
    assert "MACCS/Tanimoto" not in html
    assert "font-size=\"13\" font-weight=\"700\" fill=\"#1e293b\"" not in html
