from app import app


def test_mol_fq_returns_selected_static_sections_and_structure():
    response = app.test_client().post('/mol-fq/analyze', json={
        'smiles': 'CNCCC1=CC=CC=N1',
        'name': 'betahistine',
        'analyses': ['summary', 'structural', 'identifiers'],
    })
    assert response.status_code == 200
    payload = response.get_json()
    assert payload['module'] == 'mol_fq'
    assert payload['status'] == 'ok'
    assert payload['name'] == 'betahistine'
    assert payload['structure_svg'].startswith('<?xml')
    assert set(payload['sections']) == {'summary', 'structural', 'identifiers'}
    assert payload['sections']['structural']['values']['tpsa'] == 24.92


def test_mol_fq_ionization_exposes_explicit_method_and_states():
    response = app.test_client().post('/mol-fq/analyze', json={
        'smiles': 'CNCCC1=CC=CC=N1',
        'analyses': ['ionization'],
    })
    assert response.status_code == 200
    ionization = response.get_json()['sections']['ionization']
    assert ionization['status'] == 'estimated'
    assert 'Dimorphite-DL' in ionization['method']
    assert ionization['values']['enumerated_states']


def test_mol_fq_rejects_invalid_smiles():
    response = app.test_client().post('/mol-fq/analyze', json={
        'smiles': 'not-a-smiles',
        'analyses': ['summary'],
    })
    assert response.status_code == 400
    assert response.get_json()['status'] == 'invalid_smiles'


def test_mol_fq_ui_contains_green_module_and_analysis_checkboxes():
    html = app.test_client().get('/').get_data(as_text=True)
    assert 'data-app="fq"' in html
    assert 'theme-fq' in html
    assert 'data-fq-analysis="summary"' in html
    assert 'data-fq-analysis="nmr"' in html
    assert '/mol-fq/analyze' in html
