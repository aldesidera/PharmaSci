from app import app
from analysis import get_mol, get_properties


def test_properties_include_estimated_water_solubility_at_ph7():
    mol, error = get_mol("CCO")
    properties = get_properties(mol)

    assert error is None
    assert properties is not None
    key = "Solubilidade em água (pH 7, estimada) (mg/L)"
    assert key in properties
    assert isinstance(properties[key], float)
    assert properties[key] > 0


def test_pair_report_uses_pairwise_layout_without_batch_space():
    logd = [{"pH": pH, "LogD": round(-0.1 - abs(7 - pH) * 0.12, 3)} for pH in range(15)]
    payload = {
        "mode": "pair",
        "name_ref": "Etanol",
        "name_test": "Etilamina",
        "smiles_ref": "CCO",
        "smiles_test": "CCN",
        "similarity": 0.3571,
        "classification": "Baixa",
        "fp_type": "MACCS",
        "metric": "Tanimoto",
        "similarity_map_fingerprint": "Morgan2",
        "similarity_map_metric": "Tanimoto",
        "properties": [
                {"Propriedade": "Massa Molecular (g/mol)", "Referência": 46.07, "Teste": 45.08, "Diferença": 0.99},
            {"Propriedade": "Solubilidade em água (pH 7, estimada) (mg/L)", "Referência": 1000.0, "Teste": 2000.0, "Diferença": 1000.0},
        ],
        "logd_ref": logd,
        "logd_test": logd,
        "similarity_map": '<svg preserveAspectRatio="xMidYMid meet" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"><circle cx="5" cy="5" r="3"/></svg>',
        "fingerprint_ref_png": None,
        "fingerprint_test_png": None,
        "warnings": [],
    }
    response = app.test_client().post("/report-preview", json=payload)

    assert response.status_code == 200
    html = response.get_data(as_text=True)
    assert 'class="pair-report"' in html
    assert "pair-summary-grid" in html
    assert "pair-structure-grid" in html
    assert "Mapa de similaridade" in html
    assert "Métrica adotada:" in html
    assert "Morgan2" in html
    assert "Tanimoto" in html
    assert "preserveAspectRatio=\"xMidYMid meet\"" in html
    assert "Propriedades Físico-Químicas" in html
    assert "Indicadores" not in html[html.index('<body'):]
    assert "pair-fingerprint-inline" not in html[html.index('<body'):]
    assert "pka_heuristic_estimate" not in html[html.index('<body'):]
    assert "prova regulatória" not in html[html.index('<body'):]
    assert "pair-properties-section" in html
    assert "@page { size: A4 portrait; margin: 11mm; }" in html
    body_html = html[html.index('<body'):]
    assert body_html.count("CCO") == 1
    assert body_html.count("CCN") == 1
    assert "chemical-space-report-card" not in body_html
    assert html.index("Resumo") < html.index("Estruturas Moleculares") < html.index("Similaridade Estrutural") < html.index("Propriedades Físico-Químicas")
