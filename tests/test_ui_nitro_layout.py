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


def test_module_switch_has_local_display_fallback():
    html = TEMPLATE.read_text(encoding="utf-8")

    assert 'id="nitro-modules-card" class="card-glass p-8 shadow-sm h-auto hidden" style="display: none;"' in html
    assert "molsimControls.style.display = activeApp === 'molsim' ? '' : 'none'" in html
    assert "nitroModulesCard.style.display = activeApp === 'nitro' ? '' : 'none'" in html
