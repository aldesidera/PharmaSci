"""Transparent FDA cPCA screening for nitrosamine structures.

This module implements the structural screening flow described in the FDA
Recommended Acceptable Intake Limits for NDSRIs guidance, Appendix A and
Figure 1. It is intentionally conservative: structures outside the supported
scope return ``manual_review`` instead of an invented category.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


CPCA_RULE_VERSION = "FDA-RAIL-2023-Appendix-A"
CPCA_SOURCE_URL = "https://www.fda.gov/media/170794/download"
EMA_APPENDIX_VERSION = "EMA/42261/2025 Rev.13"
EMA_APPENDIX_UPDATED = "2026-06-24"
EMA_APPENDIX_URL = "https://www.ema.europa.eu/en/documents/other/appendix-1-acceptable-intakes-established-n-nitrosamines_en.xlsx"
EMA_APPENDIX_PATH = Path(__file__).resolve().parent / "data" / "ema_appendix1.json"

AI_LIMITS_NG_DAY: Dict[int, float] = {
    1: 26.5,
    2: 100.0,
    3: 400.0,
    4: 1500.0,
    5: 1500.0,
}

# Table A: the lower alpha-hydrogen count is listed first.
ALPHA_HYDROGEN_SCORES: Dict[Tuple[int, int], int] = {
    (0, 2): 3,
    (0, 3): 2,
    (1, 2): 3,
    (1, 3): 3,
    (2, 2): 1,
    (2, 3): 1,
}


@lru_cache(maxsize=1)
def _ema_index() -> Dict[str, Dict[str, Any]]:
    try:
        payload = json.loads(EMA_APPENDIX_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    index: Dict[str, Dict[str, Any]] = {}
    for record in payload.get("records", []):
        canonical = record.get("canonical_smiles")
        if canonical and canonical not in index:
            index[canonical] = record
    return index


def _canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _structure_svg(mol: Chem.Mol, width: int = 460, height: int = 280) -> str:
    try:
        mol_no_h = Chem.RemoveHs(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        options = drawer.drawOptions()
        options.addStereoAnnotation = True
        options.prepareMolsBeforeDrawing = True
        options.bondLineWidth = 2.2
        options.minFontSize = 14
        options.annotationFontScale = 0.85
        drawer.DrawMolecule(mol_no_h)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace(f'width="{width}px"', 'width="100%"').replace(f'height="{height}px"', 'height="100%"')
    except Exception:
        return ""


def _ema_lookup(mol: Chem.Mol, mdd_mg: Optional[float] = None) -> Dict[str, Any]:
    canonical = _canonical_smiles(mol)
    record = _ema_index().get(canonical)
    base = {
        "listed": bool(record),
        "status": "listed" if record else "not_listed",
        "reference_number": EMA_APPENDIX_VERSION,
        "last_updated": EMA_APPENDIX_UPDATED,
        "source_url": EMA_APPENDIX_URL,
        "canonical_smiles": canonical,
        "ai_unit": "ng/day",
    }
    if not record:
        base["message"] = "A estrutura não foi localizada no snapshot do Appendix 1 da EMA; nenhum AI EMA foi inferido."
        return base

    base.update(
        {
            "sheet": record.get("sheet"),
            "name": record.get("name"),
            "iupac_name": record.get("iupac_name"),
            "listed_smiles": record.get("smiles"),
            "cas_rn": record.get("cas_rn"),
            "synonym_acronym": record.get("synonym_acronym"),
            "source": record.get("source"),
            "cpca_category": record.get("cpca_category"),
            "ai_ng_day": record.get("ai_ng_day"),
            "ai_ng_day_raw": record.get("ai_ng_day_raw"),
            "note": record.get("note"),
            "publication_date": record.get("publication_date"),
            "message": "A estrutura foi localizada no snapshot do Appendix 1 da EMA.",
        }
    )
    if record.get("ai_ng_day") is None:
        base["status"] = "listed_no_numeric_ai"
    if mdd_mg is not None and record.get("ai_ng_day") is not None:
        base["ppm_limit"] = round(float(record["ai_ng_day"]) / float(mdd_mg), 6)
        base["ppm_formula"] = "AI EMA (ng/dia) / dose diária máxima (mg)"
    return base


def _add_valid_structure_fields(result: Dict[str, Any], mol: Chem.Mol, mdd_mg: Optional[float]) -> Dict[str, Any]:
    result["canonical_smiles"] = _canonical_smiles(mol)
    result["structure_svg"] = _structure_svg(mol)
    result["ema"] = _ema_lookup(mol, mdd_mg=mdd_mg)
    return result


def _base_result(smiles: str, status: str, message: str) -> Dict[str, Any]:
    return {
        "module": "cpca",
        "status": status,
        "message": message,
        "smiles": smiles,
        "rule_version": CPCA_RULE_VERSION,
        "regulatory_basis": "FDA RAIL Guidance, Appendix A and Figure 1",
        "source_url": CPCA_SOURCE_URL,
        "ai_unit": "ng/day",
        "warnings": [
            "Resultado de triagem estrutural cPCA; não substitui avaliação toxicológica, dados específicos do composto, read-across ou decisão regulatória.",
            "A versão da regra e a jurisdição devem ser verificadas antes de qualquer uso regulatório.",
        ],
    }


def _bond_is(atom_a: Chem.Atom, atom_b: Chem.Atom, bond_type: Chem.BondType) -> bool:
    bond = atom_a.GetOwningMol().GetBondBetweenAtoms(atom_a.GetIdx(), atom_b.GetIdx())
    return bond is not None and bond.GetBondType() == bond_type


def _find_n_nitroso_centers(mol: Chem.Mol) -> List[Tuple[int, int]]:
    """Return ``(amine_nitrogen_idx, nitroso_nitrogen_idx)`` pairs."""
    centers: List[Tuple[int, int]] = []
    for nitroso_n in mol.GetAtoms():
        if nitroso_n.GetAtomicNum() != 7:
            continue
        has_nitroso_oxygen = any(
            neighbor.GetAtomicNum() == 8
            and _bond_is(nitroso_n, neighbor, Chem.BondType.DOUBLE)
            for neighbor in nitroso_n.GetNeighbors()
        )
        if not has_nitroso_oxygen:
            continue
        for neighbor in nitroso_n.GetNeighbors():
            if neighbor.GetAtomicNum() != 7:
                continue
            if _bond_is(nitroso_n, neighbor, Chem.BondType.SINGLE):
                pair = (neighbor.GetIdx(), nitroso_n.GetIdx())
                if pair not in centers:
                    centers.append(pair)
    return centers


def _heavy_atom_paths(
    mol: Chem.Mol,
    start_idx: int,
    blocked: Iterable[int],
    max_depth: int = 9,
) -> Iterable[List[int]]:
    blocked_set = set(blocked)

    def visit(path: List[int]) -> Iterable[List[int]]:
        yield path
        if len(path) >= max_depth:
            return
        current = mol.GetAtomWithIdx(path[-1])
        for neighbor in current.GetNeighbors():
            idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 1 or idx in blocked_set or idx in path:
                continue
            yield from visit(path + [idx])

    return visit([start_idx])


def _path_respects_ring_limit(mol: Chem.Mol, path: Sequence[int]) -> bool:
    for ring in mol.GetRingInfo().AtomRings():
        if sum(idx in ring for idx in path) > 4:
            return False
    return True


def _has_long_chain(mol: Chem.Mol, alpha_idx: int, blocked: Iterable[int]) -> bool:
    return any(
        len(path) >= 5 and _path_respects_ring_limit(mol, path)
        for path in _heavy_atom_paths(mol, alpha_idx, blocked)
    )


def _smallest_ring_containing(mol: Chem.Mol, atom_idx: int) -> Optional[Tuple[int, ...]]:
    rings = [ring for ring in mol.GetRingInfo().AtomRings() if atom_idx in ring]
    if not rings:
        return None
    return min(rings, key=len)


def _ring_feature(mol: Chem.Mol, amine_n_idx: int) -> Optional[Dict[str, Any]]:
    ring = _smallest_ring_containing(mol, amine_n_idx)
    if ring is None:
        return None

    size = len(ring)
    symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
    hetero = {symbol for symbol in symbols if symbol not in {"C", "H"}}

    if size == 5 and symbols.count("N") == 1 and not (hetero - {"N"}):
        return {
            "name": "pyrrolidine_ring",
            "type": "deactivating",
            "score": 3,
            "detail": "Grupo N-nitroso em anel de pirrolidina de cinco membros.",
        }
    if size == 6 and "S" in symbols:
        return {
            "name": "six_membered_sulfur_ring",
            "type": "deactivating",
            "score": 3,
            "detail": "Grupo N-nitroso em anel de seis membros contendo enxofre.",
        }
    if size == 6 and symbols.count("N") == 1 and symbols.count("O") == 1:
        return {
            "name": "morpholine_ring",
            "type": "deactivating",
            "score": 1,
            "detail": "N-nitroso group in a morpholine ring.",
        }
    if size in {5, 6}:
        return {
            "name": "five_or_six_membered_ring",
            "type": "deactivating",
            "score": 2,
            "detail": f"Grupo N-nitroso em anel de {size} membros.",
        }
    if size == 7:
        return {
            "name": "seven_membered_ring",
            "type": "deactivating",
            "score": 1,
            "detail": "N-nitroso group in a seven-membered ring.",
        }
    return None


def _has_carboxylic_acid(mol: Chem.Mol) -> bool:
    pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    return bool(pattern and mol.HasSubstructMatch(pattern))


def _has_direct_double_bond_to_heteroatom(atom: Chem.Atom) -> bool:
    return any(
        neighbor.GetAtomicNum() not in {1, 6}
        and _bond_is(atom, neighbor, Chem.BondType.DOUBLE)
        for neighbor in atom.GetNeighbors()
    )


def _is_sp3(atom: Chem.Atom) -> bool:
    return atom.GetHybridization() == Chem.HybridizationType.SP3


def _is_tertiary_alpha(atom: Chem.Atom) -> bool:
    carbon_neighbors = sum(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors())
    return atom.GetAtomicNum() == 6 and _is_sp3(atom) and carbon_neighbors == 3


def _is_part_of_ethyl_group(mol: Chem.Mol, alpha: Chem.Atom, amine_n_idx: int) -> bool:
    if alpha.GetTotalNumHs() != 2:
        return False
    return any(
        neighbor.GetIdx() != amine_n_idx
        and neighbor.GetAtomicNum() == 6
        and neighbor.GetTotalNumHs() == 3
        for neighbor in alpha.GetNeighbors()
    )


def _has_beta_hydroxyl(alpha: Chem.Atom, amine_n_idx: int) -> bool:
    for beta in alpha.GetNeighbors():
        if beta.GetIdx() == amine_n_idx or beta.GetAtomicNum() != 6:
            continue
        if not _is_sp3(beta):
            continue
        for neighbor in beta.GetNeighbors():
            if neighbor.GetAtomicNum() != 8:
                continue
            if _bond_is(beta, neighbor, Chem.BondType.SINGLE) and neighbor.GetTotalNumHs() > 0:
                return True
    return False


def _has_aryl_substituent(alpha: Chem.Atom, amine_n_idx: int) -> bool:
    return any(
        neighbor.GetIdx() != amine_n_idx and neighbor.GetIsAromatic()
        for neighbor in alpha.GetNeighbors()
    )


def _has_methyl_on_beta(alpha: Chem.Atom, amine_n_idx: int) -> bool:
    for beta in alpha.GetNeighbors():
        if beta.GetIdx() == amine_n_idx or beta.GetAtomicNum() != 6 or not _is_sp3(beta):
            continue
        for neighbor in beta.GetNeighbors():
            if neighbor.GetIdx() == alpha.GetIdx() or neighbor.GetAtomicNum() != 6:
                continue
            if neighbor.GetTotalNumHs() == 3:
                return True
    return False


def _is_supported_ewg_substituent(atom: Chem.Atom) -> Optional[str]:
    """Return a conservative EWG label for supported Appendix A patterns."""
    symbol = atom.GetSymbol()
    if symbol == "C":
        has_nitrile = any(
            neighbor.GetAtomicNum() == 7
            and _bond_is(atom, neighbor, Chem.BondType.TRIPLE)
            for neighbor in atom.GetNeighbors()
        )
        return "nitrile" if has_nitrile else None
    if symbol == "S":
        double_oxygen_count = sum(
            neighbor.GetAtomicNum() == 8
            and _bond_is(atom, neighbor, Chem.BondType.DOUBLE)
            for neighbor in atom.GetNeighbors()
        )
        return "sulfonyl" if double_oxygen_count >= 2 else None
    if symbol == "N":
        has_nitro_oxygen = any(
            neighbor.GetAtomicNum() == 8
            and _bond_is(atom, neighbor, Chem.BondType.DOUBLE)
            for neighbor in atom.GetNeighbors()
        )
        return "nitro" if has_nitro_oxygen else None
    if symbol in {"F", "Cl", "Br", "I"}:
        return "halogen"
    return None


def _ewg_on_alpha(alpha: Chem.Atom, amine_n_idx: int) -> Optional[str]:
    for neighbor in alpha.GetNeighbors():
        if neighbor.GetIdx() == amine_n_idx:
            continue
        label = _is_supported_ewg_substituent(neighbor)
        if label:
            return label
    return None


def _alpha_summary(alpha: Chem.Atom, amine_n_idx: int) -> Dict[str, Any]:
    return {
        "atom_index": alpha.GetIdx(),
        "hydrogens": int(alpha.GetTotalNumHs()),
        "hybridization": str(alpha.GetHybridization()),
        "is_sp3": _is_sp3(alpha),
        "is_aromatic": bool(alpha.GetIsAromatic()),
        "is_tertiary": _is_tertiary_alpha(alpha),
        "part_of_ethyl_group": _is_part_of_ethyl_group(alpha.GetOwningMol(), alpha, amine_n_idx),
    }


def _classify_score(score: int) -> Tuple[int, float]:
    if score >= 4:
        return 4, AI_LIMITS_NG_DAY[4]
    if score == 3:
        return 3, AI_LIMITS_NG_DAY[3]
    if score == 2:
        return 2, AI_LIMITS_NG_DAY[2]
    return 1, AI_LIMITS_NG_DAY[1]


def _add_feature(features: List[Dict[str, Any]], feature: Optional[Dict[str, Any]]) -> None:
    if feature is None:
        return
    if any(item["name"] == feature["name"] for item in features):
        return
    features.append(feature)


def _analyze_center(mol: Chem.Mol, amine_n_idx: int, nitroso_n_idx: int) -> Dict[str, Any]:
    amine_n = mol.GetAtomWithIdx(amine_n_idx)
    nitroso_n = mol.GetAtomWithIdx(nitroso_n_idx)
    alpha_atoms = [
        neighbor
        for neighbor in amine_n.GetNeighbors()
        if neighbor.GetIdx() != nitroso_n_idx and neighbor.GetAtomicNum() == 6
    ]

    center: Dict[str, Any] = {
        "amine_nitrogen_index": amine_n_idx,
        "nitroso_nitrogen_index": nitroso_n_idx,
        "alpha_carbons": [_alpha_summary(atom, amine_n_idx) for atom in alpha_atoms],
        "features": [],
        "excluded_reasons": [],
        "status": "ok",
    }

    if amine_n.GetIsAromatic() or nitroso_n.GetIsAromatic():
        center["excluded_reasons"].append("N-nitroso group is within an aromatic system.")
    if len(alpha_atoms) != 2:
        center["status"] = "manual_review"
        center["message"] = "O fluxo cPCA suportado requer dois carbonos alfa diretamente ligados."
        return center
    if any(_has_direct_double_bond_to_heteroatom(atom) for atom in alpha_atoms):
        center["excluded_reasons"].append(
            "An alpha carbon is directly double-bonded to a heteroatom; the FDA cPCA scope excludes this case."
        )

    ring_feature = _ring_feature(mol, amine_n_idx)
    _add_feature(center["features"], ring_feature)
    if _has_carboxylic_acid(mol):
        _add_feature(
            center["features"],
            {
                "name": "carboxylic_acid",
                "type": "deactivating",
                "score": 3,
                "detail": "Grupo ácido carboxílico presente em qualquer parte da molécula.",
            },
        )

    ewg_labels = [_ewg_on_alpha(atom, amine_n_idx) for atom in alpha_atoms]
    ewg_labels = [label for label in ewg_labels if label]
    if ewg_labels:
        score = 2 if len(ewg_labels) == 2 else 1
        _add_feature(
            center["features"],
            {
                "name": "alpha_electron_withdrawing_group",
                "type": "deactivating",
                "score": score,
                "detail": f"Padrão EWG alfa suportado: {', '.join(ewg_labels)} em {len(ewg_labels)} lado(s).",
            },
        )

    beta_hydroxyl_sides = [
        _has_beta_hydroxyl(atom, amine_n_idx) for atom in alpha_atoms
    ]
    if any(beta_hydroxyl_sides):
        score = 2 if all(beta_hydroxyl_sides) else 1
        _add_feature(
            center["features"],
            {
                "name": "beta_hydroxyl",
                "type": "deactivating",
                "score": score,
                "detail": f"Hidroxila beta detectada em {sum(beta_hydroxyl_sides)} lado(s).",
            },
        )

    if amine_n.GetIsAromatic() is False:
        long_chain_sides = [
            _has_long_chain(mol, atom.GetIdx(), {amine_n_idx, nitroso_n_idx})
            for atom in alpha_atoms
        ]
        if all(long_chain_sides):
            _add_feature(
                center["features"],
                {
                    "name": "long_chains_both_sides",
                    "type": "deactivating",
                    "score": 1,
                    "detail": "Foram encontradas pelo menos cinco átomos pesados consecutivos nos dois lados, respeitando a restrição de percurso em anel do Appendix A.",
                },
            )

    if any(_has_aryl_substituent(atom, amine_n_idx) for atom in alpha_atoms):
        _add_feature(
            center["features"],
            {
                "name": "alpha_aryl",
                "type": "activating",
                "score": -1,
                "detail": "Substituinte aril ligado a um carbono alfa.",
            },
        )
    if any(_has_methyl_on_beta(atom, amine_n_idx) for atom in alpha_atoms):
        _add_feature(
            center["features"],
            {
                "name": "beta_methyl",
                "type": "activating",
                "score": -1,
                "detail": "Substituinte metil ligado a um carbono beta.",
            },
        )

    alpha_hydrogens = sorted(int(atom.GetTotalNumHs()) for atom in alpha_atoms)
    center["alpha_hydrogen_counts"] = alpha_hydrogens
    center["has_tertiary_alpha"] = any(_is_tertiary_alpha(atom) for atom in alpha_atoms)
    center["features"] = list(center["features"])

    if center["excluded_reasons"]:
        center["status"] = "not_applicable"
        center["message"] = "A estrutura está fora do escopo FDA cPCA suportado."
        return center

    if not any(alpha_hydrogens):
        center.update(
            {
                "potency_category": 5,
                "ai_ng_day": AI_LIMITS_NG_DAY[5],
                "category_basis": "Sem hidrogênios alfa; resultado direto de Categoria 5 da Figura 1.",
                "potency_score": None,
                "alpha_hydrogen_score": None,
                "deactivating_feature_score": 0,
                "activating_feature_score": 0,
            }
        )
        return center

    if not any(count > 1 for count in alpha_hydrogens):
        center.update(
            {
                "potency_category": 5,
                "ai_ng_day": AI_LIMITS_NG_DAY[5],
                "category_basis": "Nenhum lado possui mais de um hidrogênio alfa; resultado direto de Categoria 5 da Figura 1.",
                "potency_score": None,
                "alpha_hydrogen_score": None,
                "deactivating_feature_score": 0,
                "activating_feature_score": 0,
            }
        )
        return center

    if center["has_tertiary_alpha"]:
        center.update(
            {
                "potency_category": 5,
                "ai_ng_day": AI_LIMITS_NG_DAY[5],
                "category_basis": "Carbono alfa terciário; resultado direto de Categoria 5 da Figura 1.",
                "potency_score": None,
                "alpha_hydrogen_score": None,
                "deactivating_feature_score": 0,
                "activating_feature_score": 0,
            }
        )
        return center

    pair = tuple(alpha_hydrogens)
    alpha_score = ALPHA_HYDROGEN_SCORES.get(pair)
    if alpha_score is None:
        center["status"] = "manual_review"
        center["message"] = f"O par de hidrogênios alfa {pair} não está definido na Table A do Appendix A da FDA."
        return center

    if pair == (0, 2):
        methylene = next(atom for atom in alpha_atoms if atom.GetTotalNumHs() == 2)
        if _is_part_of_ethyl_group(mol, methylene, amine_n_idx):
            alpha_score = 2
            center["alpha_hydrogen_note"] = "Methylene alpha carbon is part of an ethyl group; FDA note applies score 2."

    deactivating_score = sum(
        int(feature["score"])
        for feature in center["features"]
        if feature["type"] == "deactivating"
    )
    activating_score = sum(
        int(feature["score"])
        for feature in center["features"]
        if feature["type"] == "activating"
    )
    potency_score = alpha_score + deactivating_score + activating_score
    category, ai_ng_day = _classify_score(potency_score)
    center.update(
        {
            "alpha_hydrogen_score": alpha_score,
            "deactivating_feature_score": deactivating_score,
            "activating_feature_score": activating_score,
            "potency_score": potency_score,
            "potency_category": category,
            "ai_ng_day": ai_ng_day,
            "category_basis": "Potency Score mapeado pela Figura 1 da FDA.",
        }
    )
    return center


def evaluate_cpca(smiles: str, mdd_mg: Optional[float] = None) -> Dict[str, Any]:
    """Evaluate a nitrosamine SMILES with a conservative FDA cPCA screen.

    Args:
        smiles: Candidate nitrosamine SMILES.
        mdd_mg: Optional maximum daily dose in mg for ppm conversion.

    Returns:
        A JSON-serializable result containing structural evidence, score,
        potency category and recommended AI when the structure is supported.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return _base_result(str(smiles), "invalid_input", "SMILES deve ser uma string não vazia.")
    if mdd_mg is not None:
        if isinstance(mdd_mg, bool) or not isinstance(mdd_mg, (int, float)) or mdd_mg <= 0:
            return _base_result(smiles, "invalid_input", "mdd_mg deve ser um número positivo quando informado.")

    normalized_smiles = smiles.strip()
    try:
        mol = Chem.MolFromSmiles(normalized_smiles)
    except Exception:
        mol = None
    if mol is None:
        return _base_result(normalized_smiles, "invalid_smiles", "O SMILES não pôde ser interpretado pelo RDKit.")

    centers = _find_n_nitroso_centers(mol)
    if not centers:
        result = _base_result(
            normalized_smiles,
            "not_nitrosamine",
            "Nenhum centro N-nitroso suportado foi detectado no SMILES.",
        )
        result.update({"center_count": 0, "centers": []})
        return _add_valid_structure_fields(result, mol, mdd_mg)

    if len(centers) > 2:
        result = _base_result(
            normalized_smiles,
            "manual_review",
            "Foram detectados mais de dois grupos N-nitroso; a Figura 1 da FDA orienta buscar orientação adicional.",
        )
        result.update({"center_count": len(centers), "centers": []})
        return _add_valid_structure_fields(result, mol, mdd_mg)

    center_results = [_analyze_center(mol, *center) for center in centers]
    unsupported = [center for center in center_results if center.get("status") != "ok"]
    if unsupported:
        result = _base_result(
            normalized_smiles,
            unsupported[0].get("status", "manual_review"),
            unsupported[0].get("message", "Pelo menos um centro nitroso requer revisão manual."),
        )
        result.update({"center_count": len(center_results), "centers": center_results})
        return _add_valid_structure_fields(result, mol, mdd_mg)

    selected = max(
        center_results,
        key=lambda center: (
            int(center.get("potency_category", 0)),
            float(center.get("ai_ng_day", 0.0)),
        ),
    )
    result = _base_result(
        normalized_smiles,
        "ok",
        "Triagem estrutural FDA cPCA concluída.",
    )
    result.update(
        {
            "center_count": len(center_results),
            "centers": center_results,
            "potency_category": selected["potency_category"],
            "potency_score": selected.get("potency_score"),
            "ai_ng_day": selected["ai_ng_day"],
            "category_basis": selected.get("category_basis"),
        }
    )
    if mdd_mg is not None:
        result["mdd_mg"] = float(mdd_mg)
        result["ppm_limit"] = round(float(selected["ai_ng_day"]) / float(mdd_mg), 6)
        result["ppm_formula"] = "AI (ng/dia) / dose diária máxima (mg)"
    return _add_valid_structure_fields(result, mol, mdd_mg)


def calculate_cpca(smiles: str, mdd_mg: Optional[float] = None) -> Dict[str, Any]:
    """Compatibility alias matching the public specification name."""
    return evaluate_cpca(smiles, mdd_mg=mdd_mg)
