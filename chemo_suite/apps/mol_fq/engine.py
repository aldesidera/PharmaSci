"""Motor inicial do Mol.FQ.

A primeira versão separa descritores estruturais calculados de propriedades
estimadas ou ainda não disponíveis. Isso evita apresentar heurísticas como se
fossem equivalentes a motores proprietários de previsão físico-química.
"""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, Iterable, Optional

from chemo_suite.core.conformer import generate_3d_conformer
from chemo_suite.core.molecule_standardization import molecule_identity

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

try:
    from dimorphite_dl import protonate_smiles
except ImportError:  # pragma: no cover - recurso opcional de ambiente
    protonate_smiles = None

from analysis import calc_logd_vs_ph, get_mol, get_properties, mol_to_svg
from rmn_suite.nmrshiftdb import RmnProviderError, predict_from_smiles
from chemo_suite.apps.mol_fq.solubility import build_solubility

ANALYSIS_KEYS = (
    "summary",
    "structural",
    "ionization",
    "lipophilicity",
    "solubility",
    "nmr",
    "geometry_3d",
)
DEFAULT_ANALYSIS_KEYS = ANALYSIS_KEYS[:-1]

LABELS = {
    "summary": "Resumo molecular",
    "structural": "Descritores estruturais",
    "ionization": "Ionização e pH",
    "lipophilicity": "Lipofilicidade",
    "solubility": "Solubilidade",
    "nmr": "H-NMR predito",
    "geometry_3d": "Geometria 3D enriquecida",
}


def _selected_set(selected: Optional[Iterable[str]]) -> set[str]:
    values = {str(item) for item in (selected or DEFAULT_ANALYSIS_KEYS)}
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


def _spectrum_payload(prediction: Any) -> Dict[str, Any]:
    peaks = list(getattr(prediction, "peaks", []) or [])
    if not peaks:
        return {"min_ppm": None, "max_ppm": None, "peaks": [], "mode": "stick"}
    ppm_values = [float(peak.ppm) for peak in peaks]
    min_ppm, max_ppm = min(ppm_values), max(ppm_values)
    span = max(max_ppm - min_ppm, 1.0)
    raw_intensities = [abs(float(peak.intensity)) if peak.intensity is not None else 1.0 for peak in peaks]
    maximum = max(raw_intensities) or 1.0
    rendered = []
    for peak, raw in zip(peaks, raw_intensities):
        x_pct = 50.0 if max_ppm == min_ppm else (max_ppm - float(peak.ppm)) / span * 100.0
        rendered.append({"ppm": float(peak.ppm), "height": round(max(0.08, raw / maximum), 4), "x_pct": round(x_pct, 4)})
    return {"min_ppm": min_ppm, "max_ppm": max_ppm, "peaks": rendered, "mode": "stick"}


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

    identity = molecule_identity(mol, smiles)
    canonical = identity["canonical_smiles"]
    properties = get_properties(mol) or {}
    structure = _structural(mol)
    composition = _composition(mol)
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
        "provenance": {
            "input_smiles": smiles,
            "canonical_smiles": canonical,
            "normalization": "RDKit canonical SMILES",
            "engine_version": "mol_fq_v4_enriched_1",
            "calculated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "requested_mode": "enriched_3d" if "geometry_3d" in selected_keys else "standard_2d",
            "standardization": identity,
        },
    }

    geometry_payload = None
    geometry_error = None
    if "geometry_3d" in selected_keys:
        geometry_payload, geometry_error = generate_3d_conformer(canonical)

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
        pka_acidic = properties.get("pKa ácido")
        pka_basic = properties.get("pKa básico")
        solubility = build_solubility(mol, canonical, pka_acidic=pka_acidic, pka_basic=pka_basic, geometry=geometry_payload)
        legacy_ph7 = properties.get("Solubilidade em água (pH 7, estimada) (mg/L)")
        solubility["legacy_ph7_mg_l"] = _status(legacy_ph7, "PharmaSci heuristic legacy", "mg/L")
        result["sections"]["solubility"] = solubility
        result["warnings"].append("A curva pH–solubilidade é uma estimativa baseada em ESOL/GSE e pKa indicativos; não representa dado experimental.")
        if solubility["calibrated"].get("status") != "calibrated":
            result["warnings"].append("Nenhum artefato calibrado aprovado está configurado; ESOL permanece como baseline exploratório.")

    if "geometry_3d" in selected_keys:
        if geometry_payload is None:
            result["sections"]["geometry_3d"] = {
                "status": "unavailable",
                "method": "RDKit ETKDGv3 + UFF",
                "message": geometry_error or "Não foi possível gerar a geometria 3D.",
            }
            result["warnings"].append(geometry_error or "Geometria 3D indisponível.")
        else:
            result["sections"]["geometry_3d"] = geometry_payload
            result["warnings"].append(geometry_payload.get("warning", "A geometria 3D é estimativa."))

    if "nmr" in selected_keys:
        nmr_results: Dict[str, Any] = {}
        nmr_warnings: list[str] = []
        for nucleus in ("1H", "13C"):
            try:
                prediction = predict_from_smiles(smiles, nucleus=nucleus)
                nmr_results[nucleus] = prediction.to_dict()
                nmr_results[nucleus]["spectrum"] = _spectrum_payload(prediction)
                nmr_warnings.extend(prediction.warnings)
            except RmnProviderError as exc:
                nmr_results[nucleus] = {
                    "status": "unavailable",
                    "provenance": "unavailable",
                    "peaks": [],
                    "warnings": [str(exc)],
                }
            except Exception as exc:  # pragma: no cover - proteção de integração externa
                nmr_results[nucleus] = {
                    "status": "unavailable",
                    "provenance": "unavailable",
                    "peaks": [],
                    "warnings": [f"Falha inesperada no provedor externo: {exc}"],
                }
        result["sections"]["nmr"] = {
            "status": "external",
            "method": "NMRShiftDB2 — espectro medido quando disponível ou previsão HOSE quando não há registro",
            "values": nmr_results,
            "message": "Predição obtida por provedor externo; não representa dado experimental quando a proveniência for hose-predicted.",
        }
        result["warnings"].append("A análise RMN envia a estrutura ao NMRShiftDB2 e depende de rede; verifique a proveniência de cada núcleo.")
        result["warnings"].extend(list(dict.fromkeys(nmr_warnings)))

    return result
