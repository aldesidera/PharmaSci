from typing import Dict


def evaluate_quantum(smiles: str) -> Dict[str, str]:
    return {
        "module": "quantum",
        "status": "not_implemented",
        "message": "Cálculo de HOMO/LUMO será integrado na próxima etapa do Nitro.RA.",
        "smiles": smiles,
    }

