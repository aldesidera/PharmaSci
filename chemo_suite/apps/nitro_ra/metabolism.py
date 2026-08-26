"""Rule-based CYP450 metabolism hypotheses for N-nitrosamines.

This module intentionally produces transparent, hypothesis-generating results. It
is not a validated toxicology predictor and does not claim enzyme-specific
regioselectivity or experimental confirmation.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdchem import Mol

from analysis import get_properties, mol_to_svg


NITROSAMINE_SMARTS = "[N;X3]-[N;X2]=O"
NITROSAMINE_QUERY = Chem.MolFromSmarts(NITROSAMINE_SMARTS)
RULE_ID = "CYP450_ALPHA_HYDROXYLATION_N_NITROSO"
ALPHA_HYDROXYLATION_SMIRKS = "[N;X2:10](=O)[N;X3:1][C;X4;H1,H2:2]>>[N;X2:10](=O)[N:1][C:2][O;H1]"
ALPHA_HYDROXYLATION_REACTION = rdChemReactions.ReactionFromSmarts(ALPHA_HYDROXYLATION_SMIRKS)
ENZYME_CONTEXT = ["CYP2E1", "CYP3A4"]
BASE_WARNING = (
    "Predição baseada em regra estrutural de triagem; a enzima predominante, "
    "a regioseletividade e o produto principal requerem confirmação experimental."
)
DIAZONIUM_WARNING = (
    "O SMILES do íon diazônio é uma representação mecanística hipotética do "
    "fragmento alquilante e não comprova a formação experimental da espécie."
)


def _empty_result(smiles: str, status: str, message: str) -> Dict[str, Any]:
    return {
        "module": "metabolism",
        "status": status,
        "message": message,
        "smiles": smiles,
        "engine": "RDKit rule-based CYP450 alpha-hydroxylation",
        "prediction_mode": "rule_based",
        "rule_id": RULE_ID,
        "reaction_smarts": ALPHA_HYDROXYLATION_SMIRKS,
        "enzyme_context": ENZYME_CONTEXT,
        "alpha_sites": [],
        "metabolites": [],
        "reactive_intermediates": [],
        "summary": {"alpha_sites": 0, "metabolites": 0, "reactive_intermediates": 0},
        "warnings": [],
    }


def _canonical_smiles(mol: Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _safe_svg(mol: Optional[Mol]) -> Optional[str]:
    if mol is None:
        return None
    try:
        return mol_to_svg(mol)
    except Exception:
        return None


def _target_payload(mol: Mol) -> Dict[str, Any]:
    canonical = _canonical_smiles(mol)
    return {
        "canonical_smiles": canonical,
        "svg": _safe_svg(mol),
        "properties": get_properties(mol) or {},
    }


def _find_nitrosamine_matches(mol: Mol) -> Tuple[Tuple[int, int, int], ...]:
    if NITROSAMINE_QUERY is None:
        return ()
    return mol.GetSubstructMatches(NITROSAMINE_QUERY, uniquify=True)


def _find_alpha_sites(mol: Mol, matches: Tuple[Tuple[int, int, int], ...]) -> List[Dict[str, int]]:
    sites: List[Dict[str, int]] = []
    seen: Set[Tuple[int, int, int]] = set()
    for amine_n_idx, nitroso_n_idx, oxygen_idx in matches:
        amine_n = mol.GetAtomWithIdx(amine_n_idx)
        for neighbor in amine_n.GetNeighbors():
            alpha_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() != 6:
                continue
            if neighbor.GetIsAromatic():
                continue
            if neighbor.GetHybridization() != Chem.HybridizationType.SP3:
                continue
            if neighbor.GetTotalNumHs() < 1:
                continue
            key = (alpha_idx, amine_n_idx, nitroso_n_idx)
            if key in seen:
                continue
            seen.add(key)
            sites.append(
                {
                    "atom_index": alpha_idx,
                    "amine_n_index": amine_n_idx,
                    "nitroso_n_index": nitroso_n_idx,
                    "nitroso_o_index": oxygen_idx,
                }
            )
    return sites


def _add_hydroxyl(mol: Mol, alpha_idx: int) -> Optional[Mol]:
    """Apply the alpha-hydroxylation SMIRKS semantics to one selected atom.

    A localized graph edit is used after the rule is selected because RDKit
    reaction matching otherwise returns every eligible alpha carbon at once,
    which is ambiguous for asymmetric nitrosamines. The public SMIRKS remains
    attached to the payload as the auditable transformation definition.
    """
    try:
        editable = Chem.RWMol(Chem.Mol(mol))
        hydroxyl = Chem.Atom(8)
        hydroxyl.SetFormalCharge(0)
        hydroxyl.SetNumExplicitHs(1)
        hydroxyl.SetNoImplicit(True)
        oxygen_idx = editable.AddAtom(hydroxyl)
        editable.AddBond(alpha_idx, oxygen_idx, Chem.BondType.SINGLE)
        product = editable.GetMol()
        Chem.SanitizeMol(product)
        return product
    except Exception:
        return None


def _remapped_index(old_idx: int, removed: Set[int]) -> int:
    return old_idx - sum(1 for index in removed if index < old_idx)


def _build_diazonium_surrogate(
    mol: Mol,
    alpha_idx: int,
    amine_n_idx: int,
    nitroso_n_idx: int,
    nitroso_o_idx: int,
) -> Optional[Mol]:
    """Build a transparent alkyl-diazonium surrogate by replacing N-nitroso N.

    The original N-nitroso nitrogen, its oxygen, and the amine nitrogen are
    removed from the parent graph. The alpha carbon is then connected to a
    formal diazonium group [N+]#N. This is deliberately labelled as a
    mechanistic surrogate because complete fragmentation can be substrate
    dependent.
    """

    removed = {amine_n_idx, nitroso_n_idx, nitroso_o_idx}
    if alpha_idx in removed:
        return None
    try:
        editable = Chem.RWMol(Chem.Mol(mol))
        for index in sorted(removed, reverse=True):
            editable.RemoveAtom(index)
        alpha_new_idx = _remapped_index(alpha_idx, removed)

        diazonium_n = Chem.Atom(7)
        diazonium_n.SetFormalCharge(1)
        diazonium_n.SetNoImplicit(True)
        terminal_n = Chem.Atom(7)
        terminal_n.SetFormalCharge(0)
        terminal_n.SetNoImplicit(True)
        diazonium_idx = editable.AddAtom(diazonium_n)
        terminal_idx = editable.AddAtom(terminal_n)
        editable.AddBond(alpha_new_idx, diazonium_idx, Chem.BondType.SINGLE)
        editable.AddBond(diazonium_idx, terminal_idx, Chem.BondType.TRIPLE)
        intermediate = editable.GetMol()
        Chem.SanitizeMol(intermediate)
        return intermediate
    except Exception:
        return None


def _site_record(mol: Mol, site: Dict[str, int], rank: int) -> Dict[str, Any]:
    alpha_idx = site["atom_index"]
    atom = mol.GetAtomWithIdx(alpha_idx)
    product = _add_hydroxyl(mol, alpha_idx)
    product_smiles = _canonical_smiles(product) if product is not None else None
    intermediate = _build_diazonium_surrogate(
        mol,
        alpha_idx,
        site["amine_n_index"],
        site["nitroso_n_index"],
        site["nitroso_o_index"],
    )
    intermediate_smiles = _canonical_smiles(intermediate) if intermediate is not None else None
    product_record: Dict[str, Any] = {
        "id": f"alpha_site_{rank}",
        "kind": "alpha_hydroxynitrosamine",
        "reaction": "alpha-hydroxylation",
        "rule_id": RULE_ID,
        "smiles": product_smiles,
        "svg": _safe_svg(product),
        "properties": get_properties(product) if product is not None else {},
    }
    intermediate_record: Dict[str, Any] = {
        "id": f"diazonium_{rank}",
        "kind": "alkyl_diazonium_surrogate",
        "smiles": intermediate_smiles,
        "svg": _safe_svg(intermediate),
        "mechanistic_status": "hypothetical",
    }
    return {
        "id": f"alpha_site_{rank}",
        "atom_index": alpha_idx,
        "atom_index_display": alpha_idx + 1,
        "carbon_hydrogens": int(atom.GetTotalNumHs()),
        "is_in_ring": bool(atom.IsInRing()),
        "carbon_degree": int(atom.GetDegree()),
        "enzyme_hypotheses": ENZYME_CONTEXT,
        "confidence": "rule_supported",
        "rule": "Carbono sp3 com hidrogênio diretamente adjacente ao N-nitroso; hipótese de α-hidroxilação CYP450.",
        "metabolite": product_record,
        "reactive_intermediate": intermediate_record,
    }


def evaluate_metabolism(smiles: str) -> Dict[str, Any]:
    """Predict CYP450 alpha-hydroxylation hypotheses for an N-nitrosamine SMILES."""

    raw_smiles = smiles if isinstance(smiles, str) else ""
    mol = Chem.MolFromSmiles(raw_smiles.strip()) if raw_smiles.strip() else None
    if mol is None:
        result = _empty_result(raw_smiles, "invalid_smiles", "SMILES inválido ou não sanitizável pelo RDKit.")
        result["warnings"] = ["Não foi possível construir a estrutura para a predição metabólica."]
        return result

    target = _target_payload(mol)
    matches = _find_nitrosamine_matches(mol)
    if not matches:
        result = _empty_result(
            raw_smiles,
            "not_nitrosamine",
            "A estrutura não contém o grupo N-nitroso necessário para esta regra de α-hidroxilação.",
        )
        result["canonical_smiles"] = target["canonical_smiles"]
        result["target"] = target
        return result

    sites = _find_alpha_sites(mol, matches)
    if not sites:
        result = _empty_result(
            raw_smiles,
            "no_alpha_sites",
            "O grupo N-nitroso foi identificado, mas nenhum carbono sp3 α com hidrogênio foi elegível para a regra.",
        )
        result["canonical_smiles"] = target["canonical_smiles"]
        result["target"] = target
        result["nitrosamine_centers"] = len(matches)
        result["warnings"] = [BASE_WARNING]
        return result

    site_records = [_site_record(mol, site, rank) for rank, site in enumerate(sites, start=1)]
    metabolites = [record["metabolite"] for record in site_records if record["metabolite"]["smiles"]]
    intermediates = [record["reactive_intermediate"] for record in site_records if record["reactive_intermediate"]["smiles"]]
    warnings = [BASE_WARNING, DIAZONIUM_WARNING]
    if len(metabolites) < len(site_records):
        warnings.append("Uma ou mais estruturas de produto não puderam ser sanitizadas e foram omitidas.")
    if len(intermediates) < len(site_records):
        warnings.append("Um ou mais fragmentos diazônio não puderam ser representados como SMILES válido.")

    return {
        "module": "metabolism",
        "status": "ok",
        "message": "Sítios α candidatos e produtos de α-hidroxilação gerados por regra estrutural.",
        "smiles": raw_smiles,
        "canonical_smiles": target["canonical_smiles"],
        "engine": "RDKit rule-based CYP450 alpha-hydroxylation",
        "prediction_mode": "rule_based",
        "rule_id": RULE_ID,
        "reaction_smarts": ALPHA_HYDROXYLATION_SMIRKS,
        "enzyme_context": ENZYME_CONTEXT,
        "target": target,
        "nitrosamine_centers": len(matches),
        "alpha_sites": site_records,
        "metabolites": metabolites,
        "reactive_intermediates": intermediates,
        "summary": {
            "alpha_sites": len(site_records),
            "metabolites": len(metabolites),
            "reactive_intermediates": len(intermediates),
        },
        "warnings": warnings,
        "disclaimer": "Resultado in silico de triagem. Não substitui avaliação toxicológica, confirmação analítica ou decisão regulatória.",
    }


def predict_cyp450_metabolism(smiles: str) -> Dict[str, Any]:
    """Public name for CYP450 metabolism prediction in the Nitro.RA module."""

    return evaluate_metabolism(smiles)
