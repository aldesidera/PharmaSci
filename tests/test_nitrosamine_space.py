from chemo_suite.apps.nitro_ra import nitrosamine_space


def _property_item(cid, smiles, title):
    return {
        "CID": cid,
        "Title": title,
        "IUPACName": title,
        "SMILES": smiles,
        "MolecularFormula": "C6H12N2O",
        "MolecularWeight": 128.17,
        "XLogP": 1.1,
        "TPSA": 32.7,
        "HBondDonorCount": 0,
        "HBondAcceptorCount": 2,
        "RotatableBondCount": 1,
    }


def test_search_filters_n_nitroso_and_builds_space(monkeypatch):
    def fake_get_json(url):
        if "/cids/JSON" in url:
            return {"IdentifierList": {"CID": [101, 102, 103]}}, None
        return {
            "PropertyTable": {
                "Properties": [
                    _property_item(101, "O=NN1CCOCC1", "N-nitrosomorpholine"),
                    _property_item(102, "CCN(CC)N=O", "N-nitrosodiethylamine"),
                    _property_item(103, "CCO", "Ethanol"),
                ]
            }
        }, None

    monkeypatch.setattr(nitrosamine_space, "_get_json", fake_get_json)
    fingerprint_types = []
    original_get_fingerprint = nitrosamine_space.get_fingerprint

    def tracked_get_fingerprint(mol, fp_type):
        fingerprint_types.append(fp_type)
        return original_get_fingerprint(mol, fp_type)

    monkeypatch.setattr(nitrosamine_space, "get_fingerprint", tracked_get_fingerprint)
    result = nitrosamine_space.search_nitrosamine_space(
        "O=NN1CCCCC1", threshold=70, max_records=50, max_candidates=10
    )

    assert result["status"] == "ok"
    assert result["search"]["retrieved_cids"] == 3
    assert result["search"]["n_nitroso_candidates"] == 2
    assert result["search"]["selected_candidates"] == 2
    assert len(result["candidates"]) == 2
    assert all(candidate["is_n_nitroso"] for candidate in result["candidates"])
    assert all(candidate["similarity"] >= 0 for candidate in result["candidates"])
    assert all(candidate["global_distance"] is not None for candidate in result["candidates"])
    assert result["points"][0]["is_target"] is True
    assert result["target"]["svg"].startswith("<?xml")
    assert fingerprint_types and set(fingerprint_types) == {"MACCS"}
    assert result["search"]["fingerprint"] == "MACCS"
    assert all(candidate["fingerprint"] == "MACCS" for candidate in result["candidates"])
    assert result["search"]["selection_method"].startswith("Filtro SMARTS [N;X3][N;X2]=O + MACCS/Tanimoto")


def test_search_handles_pubchem_failure(monkeypatch):
    monkeypatch.setattr(
        nitrosamine_space,
        "_get_json",
        lambda url: (None, {"status": "network_error", "message": "PubChem indisponível"}),
    )

    result = nitrosamine_space.search_nitrosamine_space("O=NN1CCCCC1")

    assert result["status"] == "pubchem_unavailable"
    assert result["candidates"] == []
    assert result["points"] == []
    assert result["search"]["retrieved_cids"] == 0
    assert result["search"]["n_nitroso_candidates"] == 0
    assert result["search"]["selected_candidates"] == 0
    assert result["message"] == "PubChem indisponível"


def test_search_reports_no_nitrosamine_after_filter(monkeypatch):
    def fake_get_json(url):
        if "/cids/JSON" in url:
            return {"IdentifierList": {"CID": [201]}}, None
        return {"PropertyTable": {"Properties": [_property_item(201, "CCO", "Ethanol")]}}, None

    monkeypatch.setattr(nitrosamine_space, "_get_json", fake_get_json)
    result = nitrosamine_space.search_nitrosamine_space("O=NN1CCCCC1")

    assert result["status"] == "no_nitrosamines"
    assert result["search"]["retrieved_cids"] == 1
    assert result["search"]["n_nitroso_candidates"] == 0
    assert result["search"]["selected_candidates"] == 0
    assert result["candidates"] == []
    assert len(result["points"]) == 1
    assert result["points"][0]["is_target"] is True
    assert "Nenhuma nitrosamina" in result["message"]
