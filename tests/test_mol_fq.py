from app import app


def test_mol_fq_returns_selected_static_sections_and_structure():
    response = app.test_client().post('/mol-fq/analyze', json={
        'smiles': 'CNCCC1=CC=CC=N1',
        'name': 'betahistine',
        'analyses': ['summary', 'structural', 'identifiers', 'geometry'],
    })
    assert response.status_code == 200
    payload = response.get_json()
    assert payload['module'] == 'mol_fq'
    assert payload['status'] == 'ok'
    assert payload['name'] == 'betahistine'
    assert payload['structure_svg'].startswith('<?xml')
    assert set(payload['sections']) == {'summary', 'structural'}
    assert 'identifiers' not in payload['sections']
    assert 'geometry' not in payload['sections']
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
    assert 'data-fq-analysis="identifiers"' not in html
    assert 'data-fq-analysis="geometry"' not in html
    assert '/mol-fq/analyze' in html


def test_mol_fq_nmr_preserves_external_provenance(monkeypatch):
    from chemo_suite.apps.mol_fq import engine
    from rmn_suite.nmrshiftdb import RmnPeak, RmnPrediction

    def fake_predict(smiles, nucleus='13C', **kwargs):
        return RmnPrediction(
            smiles_input=smiles,
            smiles_canonical='CCO',
            nucleus=nucleus,
            provenance='hose-predicted',
            source_url='https://example.test/nmr',
            spectrum_id=None,
            molecule_id=None,
            solvent=None,
            temperature_kelvin=None,
            field_mhz=None,
            peaks=[RmnPeak(ppm=1.23, intensity=1.0, multiplicity='s', atom_refs=['a1'])],
            warnings=['Resultado de teste: previsão HOSE.'],
        )

    monkeypatch.setattr(engine, 'predict_from_smiles', fake_predict)
    result = engine.analyze_mol_fq('CCO', selected=['nmr'])
    section = result['sections']['nmr']
    assert result['status'] == 'ok'
    assert section['status'] == 'external'
    assert section['values']['13C']['provenance'] == 'hose-predicted'
    assert section['values']['1H']['peaks'][0]['ppm'] == 1.23
    assert any('NMRShiftDB2' in warning for warning in result['warnings'])


def test_mol_fq_report_template_contains_nmr_presentation():
    html = app.test_client().get('/').get_data(as_text=True)
    assert 'fq-nmr-provenance' in html
    assert 'Previsão HOSE' in html


def test_rmn_cache_persists_and_returns_fresh_entry(monkeypatch, tmp_path):
    import rmn_suite.nmrshiftdb as provider
    from rmn_suite.nmrshiftdb import RmnPrediction

    cache_path = tmp_path / 'nmr-cache.json'
    monkeypatch.setenv('PHARMASCI_NMR_CACHE_PATH', str(cache_path))
    provider.clear_cache()
    calls = {'count': 0}

    class Response:
        text = '<cml:cml xmlns:cml="http://www.xml-cml.org/schema"><cml:spectrum id="nmrshiftdb1"><cml:peak xValue="1.2" peakMultiplicity="q"/></cml:spectrum></cml:cml>'
        def raise_for_status(self):
            return None

    def fake_get(*args, **kwargs):
        calls['count'] += 1
        return Response()

    monkeypatch.setattr(provider.requests, 'get', fake_get)
    first = provider.predict_from_smiles('CCO', nucleus='1H', cache_ttl=900)
    provider.clear_cache()
    second = provider.predict_from_smiles('CCO', nucleus='1H', cache_ttl=900)
    assert calls['count'] == 1
    assert first.cache_status == 'live'
    assert second.cache_status == 'fresh-cache'
    assert second.peaks[0].ppm == 1.2


def test_rmn_cache_uses_stale_entry_when_provider_is_unavailable(monkeypatch, tmp_path):
    import rmn_suite.nmrshiftdb as provider

    cache_path = tmp_path / 'nmr-cache.json'
    monkeypatch.setenv('PHARMASCI_NMR_CACHE_PATH', str(cache_path))
    provider.clear_cache()

    class Response:
        text = '<cml:cml xmlns:cml="http://www.xml-cml.org/schema"><cml:spectrum id="nmrshiftdb1"><cml:peak xValue="18.6" peakMultiplicity="Q"/></cml:spectrum></cml:cml>'
        def raise_for_status(self):
            return None

    monkeypatch.setattr(provider.requests, 'get', lambda *args, **kwargs: Response())
    provider.predict_from_smiles('CCO', nucleus='13C', cache_ttl=900)
    provider.clear_cache()

    def unavailable(*args, **kwargs):
        raise provider.requests.RequestException('offline')

    monkeypatch.setattr(provider.requests, 'get', unavailable)
    stale = provider.predict_from_smiles('CCO', nucleus='13C', cache_ttl=0, stale_max_age=900)
    assert stale.cache_status == 'stale-cache'
    assert stale.peaks[0].ppm == 18.6
    assert any('cache local' in warning for warning in stale.warnings)


def test_mol_fq_geometry_3d_returns_provenance_and_coordinates(monkeypatch, tmp_path):
    monkeypatch.setenv('PHARMASCI_3D_CACHE_PATH', str(tmp_path / 'geometry.json'))
    from chemo_suite.apps.mol_fq.engine import analyze_mol_fq

    result = analyze_mol_fq('CCO', selected=['summary', 'geometry_3d'], name='Etanol')
    geometry = result['sections']['geometry_3d']
    assert result['provenance']['requested_mode'] == 'enriched_3d'
    assert result['provenance']['calculated_at'].endswith('+00:00')
    assert geometry['method'] == 'RDKit ETKDGv3 + UFF'
    assert geometry['coordinates']
    assert geometry['cache_status'] == 'live'

    cached = analyze_mol_fq('CCO', selected=['geometry_3d'])['sections']['geometry_3d']
    assert cached['cache_status'] == 'fresh-cache'


def test_mol_fq_ui_contains_geometry_3d_option():
    html = app.test_client().get('/').get_data(as_text=True)
    assert 'data-fq-analysis="geometry_3d"' in html
    assert 'ETKDGv3/UFF' in html
