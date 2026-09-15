"""Motor inicial do Mol.FQ.

A primeira versão separa descritores estruturais calculados de propriedades
estimadas ou ainda não disponíveis. Isso evita apresentar heurísticas como se
fossem equivalentes a motores proprietários de previsão físico-química.
"""
from __future__ import annotations

from typing import Any, Dict, Iterable, Optional

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

try:
    from dimorphite_dl import protonate_smiles
except ImportError:  # pragma: no cover - recurso opcional de ambiente
    protonate_smiles = None

from analysis import calc_logd_vs_ph, get_mol, get_properties, mol_to_svg

ANALYSIS_KEYS = (
    "summary",
    "structural",
    "identifiers",
    "ionization",
    "lipophilicity",
    "solubility",
    "geometry",
    "nmr",
)

LABELS = {
    "summary": "Resumo molecular",
    "structural": "Descritores estruturais",
    "identifiers": "Identificadores e composição",
    "ionization": "Ionização e pH",
    "lipophilicity": "Lipofilicidade",
    "solubility": "Solubilidade",
    "geometry": "Geometria molecular",
    "nmr": "H-NMR predito",
}


def _selected_set(selected: Optional[Iterable[str]]) -> set[str]:
    values = {str(item) for item in (selected or ANALYSIS_KEYS)}
    return values.intersection(ANALYSIS_KEYS) or {"summary", "structural"}


def _composition(mol: Chem.Mol) -> Dict[str, Any]:
    formula = rdMolDescriptors.CalcMolFormula(mol)
    counts: Dict[str, int] = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1
    return {"formula": formula, "element_counts": counts}


def _identifiers(mol: Chem.Mol, smiles: str, name: Optional[str]) -> Dict[str, Any]:
    return {
        "name": name or None,
        "smiles": smiles,
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "inchi": Chem.MolToInchi(mol) if hasattr(Chem, "MolToInchi") else None,
        "inchikey": Chem.InchiToInchiKey(Chem.MolToInchi(mol)) if hasattr(Chem, "MolToInchi") else None,
        "method": "RDKit",
    }


def _structural(mol: Chem.Mol) -> Dict[str, Any]:
    return {
        "atom_count_heavy": mol.GetNumHeavyAtoms(),
        "atom_count_total_explicit_h": Chem.AddHs(mol).GetNumAtoms(),
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "ring_count": rdMolDescriptors.CalcNumRings(mol),
        "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(mol),
        "hetero_ring_count": rdMolDescriptors.CalcNumHeterocycles(mol),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 4),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "formal_charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        "tpsa": round(Descriptors.TPSA(mol), 3),
        "molar_refractivity": round(Crippen.MolMR(mol), 3),
    }


def _status(value: Any, method: str, unit: Optional[str] = None) -> Dict[str, Any]:
    return {"value": value, "unit": unit, "method": method, "status": "calculated" if value is not None else "not_available"}


def analyze_mol_fq(smiles: str, selected: Optional[Iterable[str]] = None, name: Optional[str] = None) -> Dict[str, Any]:
    selected_keys = _selected_set(selected)
    mol, error = get_mol(smiles)
    if mol is None:
        return {
            "module": "mol_fq",
            "status": "invalid_smiles",
            "smiles": smiles,
            "selected": sorted(selected_keys),
            "error": error or "SMILES inválido.",
        }

    canonical = Chem.MolToSmiles(mol, canonical=True)
    composition = _composition(mol)
    structure = _structural(mol)
    properties = get_properties(mol) or {}
    result: Dict[str, Any] = {
        "module": "mol_fq",
        "status": "ok",
        "smiles": smiles,
        "canonical_smiles": canonical,
        "name": name or None,
        "selected": sorted(selected_keys),
        "labels": {key: LABELS[key] for key in sorted(selected_keys)},
        "structure_svg": mol_to_svg(mol, size=460),
        "sections": {},
        "warnings": [],
    }

    if "summary" in selected_keys:
        result["sections"]["summary"] = {
            "molar_mass": _status(round(Descriptors.MolWt(mol), 3), "RDKit Descriptors.MolWt", "g/mol"),
            "exact_mass": _status(round(Descriptors.ExactMolWt(mol), 9), "RDKit Descriptors.ExactMolWt", "Da"),
            "formula": _status(composition["formula"], "RDKit CalcMolFormula", None),
            "formal_charge": _status(structure["formal_charge"], "RDKit atom formal charges", "e"),
            "lipinski_rule_of_five": {"status": "calculated", "method": "RDKit descriptor thresholds", "passed": None},
        }

    if "structural" in selected_keys:
        result["sections"]["structural"] = {"status": "calculated", "method": "RDKit", "values": structure}

    if "identifiers" in selected_keys:
        result["sections"]["identifiers"] = {
            "status": "calculated",
            "values": _identifiers(mol, smiles, name),
            "composition": composition,
        }

    if "ionization" in selected_keys:
        ionized_states = []
        if protonate_smiles is not None:
            try:
                ionized_states = list(protonate_smiles(canonical, ph_min=1.7, ph_max=8.0, precision=1.0, max_variants=32))
            except Exception:
                ionized_states = []
        result["sections"]["ionization"] = {
            "status": "estimated",
            "method": "Dimorphite-DL 2.0.2 para enumeração de estados + heuristic_rdkit_v1 para pKa indicativo",
            "values": {
                "pka_acidic": properties.get("pKa ácido"),
                "pka_basic": properties.get("pKa básico"),
                "enumerated_states": ionized_states,
                "ph_range": [1.7, 8.0],
            },
        }
        result["warnings"].append("Os estados ionizados são enumerados pelo Dimorphite-DL; os valores de pKa ainda são indicativos e não substituem validação experimental ou um preditor calibrado.")

    if "lipophilicity" in selected_keys:
        result["sections"]["lipophilicity"] = {
            "status": "estimated",
            "method": "RDKit Crippen MolLogP + curva LogD heurística",
            "logp": _status(properties.get("Coeficiente de Partição (LogP)"), "RDKit Crippen.MolLogP", None),
            "logd_vs_ph": calc_logd_vs_ph(mol, [0, 1.7, 4.6, 6.5, 7.0, 7.4, 8.0, 10.0, 14]),
        }
        result["warnings"].append("LogD versus pH é uma aproximação interna e não reproduz automaticamente o Chemicalize.")

    if "solubility" in selected_keys:
        result["sections"]["solubility"] = {
            "status": "estimated",
            "method": "estimativa interna de solubilidade em água em pH 7",
            "intrinsic_solubility": {"value": None, "unit": "mg/L", "status": "not_available"},
            "ph7": _status(properties.get("Solubilidade em água (pH 7, estimada) (mg/L)"), "PharmaSci heuristic solubility", "mg/L"),
        }
        result["warnings"].append("Solubilidade é estimada em pH 7; a curva completa e a solubilidade intrínseca ainda não estão implementadas.")

    if "geometry" in selected_keys:
        result["sections"]["geometry"] = {
            "status": "partial",
            "method": "RDKit MMFF conformer disponível; métricas Chemicalize ainda não implementadas",
            "values": {"van_der_waals_volume": None, "van_der_waals_surface": None, "sasa": None},
        }
        result["warnings"].append("Métricas geométricas avançadas requerem implementação adicional e não são inferidas nesta versão.")

    if "nmr" in selected_keys:
        result["sections"]["nmr"] = {
            "status": "not_implemented",
            "method": None,
            "values": None,
            "message": "Predição H-NMR ainda não está disponível no motor local.",
        }

    return result
