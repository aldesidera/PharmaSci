from app import app
from chemo_suite.apps.nitro_ra.cpca import calculate_cpca
from chemo_suite.apps.nitro_ra.metabolism import evaluate_metabolism


def test_nitro_report_layout_renders_cpca_ema_and_metabolism_contracts():
    smiles = "O=NN1CCCCC1"
    payload = {
        "module": "nitro_ra",
        "mode": "nitro",
        "status": "ok",
        "smiles": smiles,
        "name": "N-nitroso-piperidina",
        "generated_at": "27/08/2026 12:30",
        "modules": ["cpca", "metabolism"],
        "results": {
            "cpca": calculate_cpca(smiles, mdd_mg=10),
            "metabolism": evaluate_metabolism(smiles),
        },
    }

    response = app.test_client().post("/report-preview", json=payload)

    assert response.status_code == 200
    html = response.get_data(as_text=True)
    assert "Relatório de Análise Nitro.RA" in html
    assert "nitro-target-card" in html
    assert html.index(">Molécula analisada<") < html.index(">Módulos selecionados<")
    assert "Resultado estrutural cPCA" in html
    assert "Limite de ingestão aceitável — Apêndice I da EMA" in html
    assert "nitro-evidence-table" in html
    assert "nitro-metabolism-table" in html
    assert "nitro-simulation-grid" in html
    assert "Simulação 1 · α-hidroxilação" in html
    assert "Simulação 2 · intermediário diazônio" in html
    assert "Metabólitos hipotéticos gerados pela adição de OH" in html
    assert "Surrogates mecanísticos hipotéticos" in html
    metabolism_html = html.split('id="report-metabolism"', 1)[1].split('id="report-', 1)[0]
    assert "nitro-product-grid" not in metabolism_html
    assert "Sítios alfa elegíveis" in html
    assert "CYP2E1" in html
    assert "CYP3A4" in html
    assert "nitro-metabolism-summary-grid" in html
    assert "nitro-metabolism-sites" in html
    assert "mechanistic_status" not in html
    assert "Quantum" not in html
    assert "Carcinogenic Potency Categorisation Approach (CPCA)" in html
    assert "Apêndice 2 da EMA" in html
    assert "https://www.ema.europa.eu/en/documents/other/appendix-2-carcinogenic-potency-categorisation-approach-n-nitrosamines_en.pdf" in html
    assert "Estado da análise" not in html
    cpca_html = html.split('id="report-cpca"', 1)[1].split('id="report-', 1)[0]
    assert "Avaliação estrutural CPCA concluída." not in cpca_html
    assert '<span class="status-badge' not in cpca_html
    assert "body.nitro-report .summary-grid.nitro-summary-grid { grid-template-columns: minmax(0, 1fr) minmax(7.25rem, 8.5rem);" in html
    assert "nitro-modules-box { align-items: center; justify-self: end; width: 100%; max-width: 8.5rem; text-align: center; }" in html
    assert "nitro-analyzed-box { align-items: flex-start; text-align: left; }" in html
    assert "min-height: 48px; height: 100%" in html
    assert "Potency Score" not in html
    assert "AI FDA cPCA" not in html
    assert "Fonte estrutural" not in html
    assert "A estrutura é apresentada como referência comum" not in html
    assert "Deep-PK" not in html
    assert html.index("id=\"report-cpca\"") < html.index("id=\"report-metabolism\"")
