from pathlib import Path


TEMPLATE = Path(__file__).parents[1] / "templates" / "index.html"


def test_nitro_module_order_and_labels_are_stable():
    html = TEMPLATE.read_text(encoding="utf-8")
    checks = html[html.index('<div class="nitro-module-checks">'):]
    tabs = html[html.index('<div class="nitro-tabs"'):]

    check_order = [
        'data-nitro-module="cpca"',
        'data-nitro-module="nitrosamine_space"',
        'data-nitro-module="quantum"',
        'data-nitro-module="metabolism"',
    ]
    tab_order = [
        'data-nitro-tab="cpca"',
        'data-nitro-tab="nitrosamine_space"',
        'data-nitro-tab="quantum"',
        'data-nitro-tab="metabolism"',
    ]

    assert [checks.index(value) for value in check_order] == sorted(checks.index(value) for value in check_order)
    assert [tabs.index(value) for value in tab_order] == sorted(tabs.index(value) for value in tab_order)
    assert "Espaço Químico</strong><small>PubChem e similaridade" in checks
    assert "Espaço Químico</button>" in tabs
    assert "<th>RotB</th>" in html
    assert 'id="nitro-space-chart"' in html
    assert 'id="nitro-ema-space-chart"' in html
    assert "PubChem — vizinhança relativa" in html
    assert "EMA — Apêndice I · referência" in html
    assert "10 menores distâncias exibidas" in html


def test_module_switch_has_local_display_fallback():
    html = TEMPLATE.read_text(encoding="utf-8")

    assert 'id="nitro-modules-card" class="card-glass p-8 shadow-sm h-auto hidden" style="display: none;"' in html
    assert "molsimControls.style.display = activeApp === 'molsim' ? '' : 'none'" in html
    assert "nitroModulesCard.style.display = activeApp === 'nitro' ? '' : 'none'" in html


def test_molsim_chemical_space_is_explicitly_batch_only_and_independent():
    html = TEMPLATE.read_text(encoding="utf-8")
    batch_section = html[html.index('id="batch-mode"'):html.index('id="nitro-modules-card"')]
    pair_section = html[html.index('id="pair-mode"'):html.index('id="batch-mode"')]
    nitro_section = html[html.index('id="nitro-modules-card"'):html.index('id="result-content"')]

    assert 'id="batch-show-chemical-space"' in batch_section
    assert 'id="batch-show-chemical-space"' not in pair_section
    assert 'id="batch-show-chemical-space"' not in nitro_section
    assert 'data-nitro-module="cpca"' in nitro_section
    assert 'data-nitro-module="quantum"' in nitro_section
    assert 'data-nitro-module="metabolism"' in nitro_section
    assert 'id="nitro-module-results"' in html


def test_molsim_space_text_excludes_external_sources():
    html = TEMPLATE.read_text(encoding="utf-8")
    space_start = html.index('id="batch-chemical-space-panel"')
    space_end = html.index('id="batch-properties-panel"')
    space_section = html[space_start:space_end]
    assert "não consulta PubChem nem EMA" in space_section


def test_logo_switch_has_explicit_active_and_inactive_visibility_rules():
    html = TEMPLATE.read_text(encoding="utf-8")
    assert ".brand-logo-image" in html
    assert "display: none;" in html
    assert ".brand-visual.is-molsim .brand-logo-image--molsim" in html
    assert ".brand-visual.is-nitro .brand-logo-image--nitro" in html
    assert ".brand-visual.is-molsim .brand-logo-image--nitro" in html
    assert ".brand-visual.is-nitro .brand-logo-image--molsim" in html
    assert "if (brandVisual)" in html
    assert "brandVisual.classList.toggle('is-nitro', activeApp === 'nitro')" in html
