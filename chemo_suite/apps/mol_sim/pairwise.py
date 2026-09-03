from typing import Any, Callable, Dict, Optional, Tuple


def run_pairwise_compare(
    data: Dict[str, Any],
    compare_fn: Callable[[str, str, str, str, str, str, bool], Tuple[Optional[Dict[str, Any]], Optional[str]]],
) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    fp_type = data.get("fp_type", "Morgan2")
    show_similarity_map = data.get("show_similarity_map")
    show_logd = data.get("show_logd")
    # O heatmap é parte obrigatória da comparação par a par; flags legadas não podem suprimi-lo.
    show_map = True

    return compare_fn(
        data["smiles_ref"],
        data["smiles_test"],
        data.get("name_ref") or "Molécula Referência",
        data.get("name_test") or "Molécula Teste",
        fp_type,
        data["metric"],
        show_map=show_map,
    )

