from typing import Dict, Optional, Tuple

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except Exception:  # pragma: no cover - optional dependency guard
    Chem = None
    AllChem = None


def generate_3d_conformer(smiles: str, max_iters: int = 500) -> Tuple[Optional[Dict[str, float]], Optional[str]]:
    if Chem is None or AllChem is None:
        return None, "RDKit não está disponível para geração conformacional."
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "SMILES inválido para geração conformacional."
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        return None, "Falha ao gerar conformação 3D."
    ff_status = AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    return {"embed_status": float(status), "uff_status": float(ff_status)}, None

