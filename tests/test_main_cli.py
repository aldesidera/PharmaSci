import json

from main import main


def test_main_cli_cpca_isolated(capsys):
    exit_code = main([
        "--app",
        "nitro_ra",
        "--module",
        "cpca",
        "--smiles",
        "O=NN1CCCCC1",
        "--mdd-mg",
        "10",
    ])

    assert exit_code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "ok"
    assert payload["potency_category"] == 3
    assert payload["ai_ng_day"] == 400.0
    assert payload["ppm_limit"] == 40.0
