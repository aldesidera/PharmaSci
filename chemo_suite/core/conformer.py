"""Geração 3D leve e auditável para o modo enriquecido do Mol.FQ."""
from __future__ import annotations

import hashlib
import json
import math
import os
import time
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
except Exception:  # pragma: no cover - dependência opcional
    Chem = None
    AllChem = None
    Descriptors = None
    rdMolDescriptors = None


def _finite(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def generate_3d_conformer(smiles: str, max_iters: int = 500) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """Gera um conformero ETKDGv3, otimiza com UFF e retorna coordenadas auditáveis.

    O método é deliberadamente leve: não representa cálculo quântico e não deve ser
    interpretado como geometria otimizada por DFT. Os status antigos continuam no
    payload para preservar compatibilidade com consumidores existentes.
    """
    if Chem is None or AllChem is None:
        return None, "RDKit não está disponível para geração conformacional."
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "SMILES inválido para geração conformacional."
    canonical = Chem.MolToSmiles(mol, canonical=True)
    cache_path = Path(os.getenv("PHARMASCI_3D_CACHE_PATH", ".cache/pharmasci_3d_geometry.json"))
    ttl = _finite(os.getenv("PHARMASCI_3D_CACHE_TTL_SECONDS", "604800")) or 604800.0
    cache_key = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    try:
        cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}
        cached = cache.get(cache_key)
        if isinstance(cached, dict) and (time.time() - float(cached.get("cached_at", 0))) <= ttl:
            payload = dict(cached.get("payload") or {})
            payload["cache_status"] = "fresh-cache"
            payload["cache_key"] = cache_key
            return payload, None
    except Exception:
        cache = {}
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20240916
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        return None, "Falha ao gerar conformação 3D."
    try:
        ff_status = AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    except Exception as exc:  # pragma: no cover - varia por molécula/ambiente
        ff_status = -1
        optimization_error = str(exc)
    else:
        optimization_error = None

    conformer = mol.GetConformer()
    coordinates = []
    for atom in mol.GetAtoms():
        point = conformer.GetAtomPosition(atom.GetIdx())
        coordinates.append({
            "atom_index": atom.GetIdx(),
            "element": atom.GetSymbol(),
            "x": round(float(point.x), 5),
            "y": round(float(point.y), 5),
            "z": round(float(point.z), 5),
        })

    distances = []
    for left in range(mol.GetNumAtoms()):
        for right in range(left + 1, mol.GetNumAtoms()):
            distance = ((coordinates[left]["x"] - coordinates[right]["x"]) ** 2 +
                        (coordinates[left]["y"] - coordinates[right]["y"]) ** 2 +
                        (coordinates[left]["z"] - coordinates[right]["z"]) ** 2) ** 0.5
            if distance <= 4.0:
                distances.append(float(distance))

    payload = {
        "embed_status": float(status),
        "uff_status": float(ff_status),
        "optimization_error": optimization_error,
        "method": "RDKit ETKDGv3 + UFF",
        "status": "calculated" if ff_status == 0 else "estimated",
        "cache_status": "live",
        "cache_key": cache_key,
        "atom_count_3d": mol.GetNumAtoms(),
        "heavy_atom_count_3d": sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() != "H"),
        "conformer_id": int(conformer.GetId()),
        "coordinate_unit": "angstrom",
        "coordinates": coordinates,
        "distance_summary": {
            "min_angstrom": round(min(distances), 5) if distances else None,
            "max_angstrom": round(max(distances), 5) if distances else None,
            "mean_angstrom": round(sum(distances) / len(distances), 5) if distances else None,
        },
        "warning": "Geometria 3D gerada por ETKDGv3/UFF; não é cálculo quântico nem substitui otimização QM.",
    }
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}
        cache[cache_key] = {"cached_at": time.time(), "payload": payload}
        cache_path.write_text(json.dumps(cache, ensure_ascii=False), encoding="utf-8")
    except Exception:
        pass
    return payload, None
