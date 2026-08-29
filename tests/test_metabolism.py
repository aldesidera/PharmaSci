from chemo_suite.apps.nitro_ra.metabolism import evaluate_metabolism, predict_cyp450_metabolism


def test_public_predictor_alias_uses_the_same_contract():
    result = predict_cyp450_metabolism("O=NN1CCCCC1")

    assert result["status"] == "ok"
    assert result["module"] == "metabolism"


def test_predicts_alpha_sites_metabolites_and_diazonium_surrogates():
    result = evaluate_metabolism("O=NN1CCCCC1")

    assert result["module"] == "metabolism"
    assert result["status"] == "ok"
    assert result["prediction_mode"] == "rule_based"
    assert result["rule_id"] == "CYP450_ALPHA_HYDROXYLATION_N_NITROSO"
    assert result["enzyme_context"] == ["CYP2E1", "CYP3A4"]
    assert result["summary"]["alpha_sites"] == 2
    assert result["summary"]["metabolites"] == 2
    assert result["summary"]["reactive_intermediates"] == 2
    assert len(result["alpha_sites"]) == 2

    for site in result["alpha_sites"]:
        assert site["atom_index_display"] > 0
        assert site["confidence"] == "rule_supported"
        assert site["rule_match"] is True
        assert site["rule_match_label"] == "Sítio identificado"
        assert site["metabolite"]["kind"] == "alpha_hydroxynitrosamine"
        assert site["metabolite"]["smiles"]
        assert site["reactive_intermediate"]["kind"] == "alkyl_diazonium_surrogate"
        assert site["reactive_intermediate"]["mechanistic_status"] == "hypothetical"
        assert site["reactive_intermediate"]["smiles"]

    assert "experimental" in result["warnings"][0]
    assert "hipotética" in result["warnings"][1]
    assert result["disclaimer"]


def test_non_nitrosamine_is_explicitly_out_of_scope():
    result = evaluate_metabolism("CCO")

    assert result["status"] == "not_nitrosamine"
    assert result["alpha_sites"] == []
    assert result["metabolites"] == []
    assert result["reactive_intermediates"] == []
    assert result["target"]["canonical_smiles"] == "CCO"


def test_invalid_smiles_is_safe_and_structured():
    result = evaluate_metabolism("not-a-smiles")

    assert result["status"] == "invalid_smiles"
    assert result["alpha_sites"] == []
    assert result["metabolites"] == []
    assert result["reactive_intermediates"] == []
    assert result["warnings"]


def test_aromatic_n_nitroso_imidazole_is_recognized_without_alpha_sp3_sites():
    result = evaluate_metabolism("C1=CN(C=N1)N=O")

    assert result["status"] == "no_alpha_sites"
    assert result["canonical_smiles"] == "O=Nn1ccnc1"
    assert result["nitrosamine_centers"] == 1
    assert result["alpha_sites"] == []
    assert result["metabolites"] == []
    assert result["reactive_intermediates"] == []


def test_nitrosamine_with_heteroatom_in_ring_keeps_alpha_sites():
    result = evaluate_metabolism("O=NN1CCOCC1")

    assert result["status"] == "ok"
    assert result["summary"]["alpha_sites"] == 2
    assert all(site["metabolite"]["smiles"] for site in result["alpha_sites"])


def test_asymmetric_nitrosamine_keeps_site_specific_metabolites():
    result = evaluate_metabolism("O=NN(CC)CCO")

    assert result["status"] == "ok"
    products = [site["metabolite"]["smiles"] for site in result["alpha_sites"]]
    assert len(products) == 2
    assert len(set(products)) == 2
