from typing import Dict


def evaluate_metabolism(smiles: str) -> Dict[str, str]:
    return {
        "module": "metabolism",
        "status": "not_implemented",
        "message": "Predição CYP450 será integrada na próxima etapa do Nitro.RA.",
        "smiles": smiles,
    }

