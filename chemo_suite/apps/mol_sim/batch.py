from typing import Any, Callable, Dict, List, Optional, Tuple


def _resolve_names(
    names_list: Optional[List[Any]],
    smiles_list: List[str],
    lookup_name_fn: Callable[[Optional[str]], Optional[str]],
) -> Optional[List[str]]:
    if not isinstance(names_list, list):
        return names_list
    resolved_names: List[str] = []
    for index, name in enumerate(names_list):
        candidate = name if isinstance(name, str) and name.strip() else None
        if candidate:
            resolved_names.append(candidate)
            continue
        smiles = smiles_list[index] if index < len(smiles_list) else None
        pubchem_name = lookup_name_fn(smiles)
        resolved_names.append(pubchem_name or f"Mol_{index + 1}")
    return resolved_names


def run_batch_compare(
    data: Dict[str, Any],
    bulk_compare_fn: Callable[[str, List[str], Optional[List[str]], str, str], Tuple[Optional[List[Dict[str, Any]]], Optional[str]]],
    lookup_name_fn: Callable[[Optional[str]], Optional[str]],
) -> Tuple[Optional[List[Dict[str, Any]]], Optional[str]]:
    names_list = _resolve_names(data.get("names_list"), data["smiles_list"], lookup_name_fn)
    data["names_list"] = names_list
    return bulk_compare_fn(
        data["ref_smiles"],
        data["smiles_list"],
        data.get("names_list"),
        data.get("fp_type", "Morgan2"),
        data["metric"],
    )

