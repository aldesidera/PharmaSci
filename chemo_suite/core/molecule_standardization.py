"""Política compartilhada de identidade e padronização molecular.

A política é deliberadamente conservadora. Ela corrige/sanitiza a estrutura e
produz identificadores reprodutíveis, mas não remove sais, não escolhe
sistematicamente tautômeros e não neutraliza cargas sem solicitação explícita.
Isso evita alterar silenciosamente a molécula usada em comparações ou avaliações
regulatórias.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import inchi

POLICY_VERSION = "conservative_rdkit_v1"


def _canonical(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _inchikey(mol: Chem.Mol) -> Optional[str]:
    try:
        value = inchi.MolToInchi(mol)
        return inchi.InchiToInchiKey(value) if value else None
    except Exception:
        return None


def standardize_molecule(smiles: str, *, preserve_stereo: bool = True) -> Tuple[Optional[Chem.Mol], Dict[str, Any], Optional[str]]:
    """Parseia, sanitiza e retorna identidade/proveniência sem transformação destrutiva."""
    original = (smiles or "").strip() if isinstance(smiles, str) else ""
    if not original:
        return None, {"input_smiles": original, "policy_version": POLICY_VERSION}, "SMILES vazio."
    try:
        mol = Chem.MolFromSmiles(original, sanitize=True)
        if mol is None:
            return None, {"input_smiles": original, "policy_version": POLICY_VERSION}, "Não foi possível interpretar o SMILES."
        Chem.SanitizeMol(mol)
        canonical = _canonical(mol)
        canonical_no_stereo = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        components = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True))
        metadata: Dict[str, Any] = {
            "input_smiles": original,
            "canonical_smiles": canonical,
            "canonical_smiles_no_stereo": canonical_no_stereo,
            "inchikey": _inchikey(mol),
            "policy_version": POLICY_VERSION,
            "stereo_preserved": bool(preserve_stereo),
            "fragment_count": len(components),
            "has_multiple_fragments": len(components) > 1,
            "formal_charge": int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
            "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
            "transformations": ["trim_whitespace", "rdkit_parse", "rdkit_sanitize", "canonicalize_identity"],
            "salt_removal_applied": False,
            "tautomer_normalization_applied": False,
            "charge_normalization_applied": False,
        }
        return mol, metadata, None
    except Exception as exc:
        return None, {"input_smiles": original, "policy_version": POLICY_VERSION}, f"Falha na padronização conservadora: {exc}"


def molecule_identity(mol: Chem.Mol, input_smiles: Optional[str] = None) -> Dict[str, Any]:
    """Produz identidade para uma molécula já carregada."""
    canonical = _canonical(mol)
    fragments = Chem.GetMolFrags(mol)
    return {
        "input_smiles": input_smiles,
        "canonical_smiles": canonical,
        "canonical_smiles_no_stereo": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        "inchikey": _inchikey(mol),
        "policy_version": POLICY_VERSION,
        "fragment_count": len(fragments),
        "has_multiple_fragments": len(fragments) > 1,
        "formal_charge": int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
        "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
        "salt_removal_applied": False,
        "tautomer_normalization_applied": False,
        "charge_normalization_applied": False,
    }
