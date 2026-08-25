from chemo_suite.apps.mol_sim.batch import run_batch_compare
from chemo_suite.apps.mol_sim.pairwise import run_pairwise_compare
from chemo_suite.apps.nitro_ra.cpca import evaluate_cpca
from chemo_suite.apps.nitro_ra.metabolism import evaluate_metabolism
from chemo_suite.apps.nitro_ra.quantum import evaluate_quantum
from chemo_suite.core.parser import sanitize_smiles, sanitize_smiles_list


def test_pairwise_module_preserves_show_logd_fallback():
    captured = {}

    def fake_compare(smiles_ref, smiles_test, name_ref, name_test, fp_type, metric, show_map=True):
        captured["show_map"] = show_map
        return {"ok": True}, None

    data = {
        "smiles_ref": "CCO",
        "smiles_test": "CCN",
        "name_ref": "A",
        "name_test": "B",
        "metric": "Tanimoto",
        "fp_type": "Morgan2",
        "show_logd": False,
    }
    result, error = run_pairwise_compare(data, fake_compare)
    assert error is None
    assert result == {"ok": True}
    assert captured["show_map"] is False


def test_batch_module_resolves_blank_names():
    def fake_bulk(ref_smiles, smiles_list, names_list, fp_type, metric):
        return [{"name": names_list[0]}], None

    def fake_lookup(smiles):
        return "Resolved"

    data = {
        "ref_smiles": "CCO",
        "smiles_list": ["CCN"],
        "names_list": [""],
        "fp_type": "Morgan2",
        "metric": "Tanimoto",
    }
    result, error = run_batch_compare(data, fake_bulk, fake_lookup)
    assert error is None
    assert result[0]["name"] == "Resolved"


def test_parser_sanitize_smiles():
    smiles, err = sanitize_smiles(" CCO ")
    assert err is None
    assert smiles == "CCO"


def test_parser_sanitize_smiles_list():
    smiles_list, err = sanitize_smiles_list([" CCO ", "CCN"])
    assert err is None
    assert smiles_list == ["CCO", "CCN"]


def test_nitro_ra_scaffold_contracts():
    cpca = evaluate_cpca("CCO")
    quantum = evaluate_quantum("CCO")
    metabolism = evaluate_metabolism("CCO")

    assert cpca["status"] == "not_nitrosamine"
    assert cpca["module"] == "cpca"
    assert quantum["status"] == "not_implemented"
    assert metabolism["status"] == "not_implemented"

