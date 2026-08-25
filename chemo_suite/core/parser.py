from typing import Any, Dict, List, Optional, Tuple

try:
    from rdkit import Chem
except Exception:  # pragma: no cover - optional dependency guard
    Chem = None


def sanitize_smiles(value: Any, max_length: int = 4096) -> Tuple[Optional[str], Optional[str]]:
    if not isinstance(value, str):
        return None, "SMILES deve ser string."
    normalized = value.strip()
    if not normalized:
        return None, "SMILES não pode ser vazio."
    if len(normalized) > max_length:
        return None, f"SMILES excede o limite de {max_length} caracteres."
    return normalized, None


def sanitize_smiles_list(values: Any, max_items: int = 100, max_length: int = 4096) -> Tuple[Optional[List[str]], Optional[str]]:
    if not isinstance(values, list):
        return None, "smiles_list deve ser uma lista."
    if not values:
        return None, "smiles_list deve ser uma lista não vazia."
    if len(values) > max_items:
        return None, f"smiles_list excede o limite de {max_items} itens."
    sanitized: List[str] = []
    for index, item in enumerate(values):
        normalized, err = sanitize_smiles(item, max_length=max_length)
        if err:
            return None, f"Item {index} inválido em smiles_list: {err}"
        sanitized.append(normalized)
    return sanitized, None


def to_rdkit_mol(smiles: str):
    if Chem is None:
        return None
    return Chem.MolFromSmiles(smiles)


def parse_smiles_payload(payload: Dict[str, Any], field: str, max_length: int = 4096) -> Tuple[Optional[str], Optional[str]]:
    value = payload.get(field)
    return sanitize_smiles(value, max_length=max_length)

