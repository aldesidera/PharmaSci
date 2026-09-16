from analysis import get_mol
from chemo_suite.core.molecule_standardization import POLICY_VERSION, standardize_molecule
from chemo_suite.apps.nitro_ra.cpca import evaluate_cpca
from chemo_suite.apps.nitro_ra.metabolism import evaluate_metabolism
from chemo_suite.apps.mol_fq.engine import analyze_mol_fq


def test_standardization_is_conservative_and_reproducible():
    mol, identity, error = standardize_molecule(' CCO ')
    assert error is None
    assert mol is not None
    assert identity['input_smiles'] == 'CCO'
    assert identity['canonical_smiles'] == 'CCO'
    assert identity['policy_version'] == POLICY_VERSION
    assert identity['salt_removal_applied'] is False
    assert identity['tautomer_normalization_applied'] is False


def test_shared_loader_attaches_identity_without_changing_molecule():
    mol, error = get_mol(' CCO ')
    assert error is None
    assert mol is not None
    assert mol.GetProp('_pharmasci_policy_version') == POLICY_VERSION
    assert mol.GetProp('_pharmasci_canonical_smiles') == 'CCO'


def test_molfq_reports_standardization_provenance():
    result = analyze_mol_fq(' CCO ', selected=['summary'])
    assert result['provenance']['standardization']['canonical_smiles'] == 'CCO'
    assert result['provenance']['standardization']['salt_removal_applied'] is False


def test_nitrora_engines_report_standardization_provenance():
    cpca = evaluate_cpca('CCN(CC)N=O')
    metabolism = evaluate_metabolism('CCN(CC)N=O')
    assert 'standardization' in cpca
    assert cpca['standardization']['policy_version'] == POLICY_VERSION
    assert 'standardization' in metabolism
    assert metabolism['standardization']['policy_version'] == POLICY_VERSION
