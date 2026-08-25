"""Optional PubChem-backed nitrosamine chemical-space search for Nitro.RA."""

from __future__ import annotations

import json
import logging
import math
import os
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple
from urllib import error as urllib_error
from urllib import parse as urllib_parse
from urllib import request as urllib_request

import numpy as np
from rdkit import Chem

from analysis import calc_similarity, classify_similarity, get_fingerprint, get_properties, mol_to_svg

logger = logging.getLogger(__name__)

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_SOURCE_URL = "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest"
DEFAULT_THRESHOLD = int(os.getenv("MOLSIM_PUBCHEM_SPACE_THRESHOLD", "50"))
DEFAULT_MAX_RECORDS = int(os.getenv("MOLSIM_PUBCHEM_SPACE_MAX_RECORDS", "100"))
DEFAULT_MAX_CANDIDATES = int(os.getenv("MOLSIM_PUBCHEM_SPACE_MAX_CANDIDATES", "10"))
DEFAULT_TIMEOUT = float(os.getenv("MOLSIM_PUBCHEM_SPACE_TIMEOUT", "12"))
DEFAULT_CACHE_TTL = float(os.getenv("MOLSIM_PUBCHEM_SPACE_CACHE_TTL_SECONDS", "3600"))

N_NITROSO_PATTERN = Chem.MolFromSmarts("[N;X3][N;X2]=O")
SPACE_FINGERPRINT = "MACCS"
DESCRIPTOR_KEYS = (
    "Massa Molecular (g/mol)",
    "Coeficiente de Partição (LogP)",
    "Área de Superfície Polar (Å²)",
    "Doadores de H (HBD)",
    "Receptores de H (HBA)",
)

_HTTP_CACHE: Dict[str, Tuple[float, Any]] = {}


def _cache_get(url: str) -> Any:
    cached = _HTTP_CACHE.get(url)
    if not cached:
        return None
    created, value = cached
    if time.monotonic() - created > DEFAULT_CACHE_TTL:
        _HTTP_CACHE.pop(url, None)
        return None
    return value


def _cache_put(url: str, value: Any) -> None:
    _HTTP_CACHE[url] = (time.monotonic(), value)


def _get_json(url: str) -> Tuple[Optional[Any], Optional[Dict[str, Any]]]:
    cached = _cache_get(url)
    if cached is not None:
        return cached, None
    request = urllib_request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "PharmaSci-NitroRA/1.0 (PubChem chemical-space feature)",
        },
    )
    try:
        with urllib_request.urlopen(request, timeout=DEFAULT_TIMEOUT) as response:
            payload = json.loads(response.read().decode("utf-8", errors="replace"))
        _cache_put(url, payload)
        return payload, None
    except urllib_error.HTTPError as exc:
        return None, {"status": "http_error", "http_status": exc.code, "message": f"PubChem retornou HTTP {exc.code}."}
    except (urllib_error.URLError, TimeoutError) as exc:
        return None, {"status": "network_error", "message": f"Não foi possível consultar o PubChem: {exc.reason if hasattr(exc, 'reason') else exc}."}
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        return None, {"status": "invalid_response", "message": f"A resposta do PubChem não era JSON válido: {exc}."}
    except Exception as exc:
        logger.warning("Falha inesperada no PubChem: %s", exc, exc_info=True)
        return None, {"status": "network_error", "message": "Falha inesperada ao consultar o PubChem."}


def _base_result(smiles: str, status: str, message: str) -> Dict[str, Any]:
    return {
        "module": "nitrosamine_space",
        "status": status,
        "message": message,
        "smiles": smiles,
        "source_url": PUBCHEM_SOURCE_URL,
        "warnings": [
            "A busca depende da disponibilidade e da cobertura do PubChem; ausência no lote não prova ausência no universo químico.",
            "A similaridade 2D e os descritores são ferramentas de triagem e não substituem avaliação toxicológica ou regulatória.",
        ],
    }


def _canonical(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _has_n_nitroso(mol: Chem.Mol) -> bool:
    return bool(N_NITROSO_PATTERN and mol.HasSubstructMatch(N_NITROSO_PATTERN))


def _descriptor_vector(properties: Optional[Dict[str, Any]]) -> Optional[np.ndarray]:
    if not properties:
        return None
    values: List[float] = []
    for key in DESCRIPTOR_KEYS:
        value = properties.get(key)
        if isinstance(value, bool):
            return None
        try:
            value_float = float(value)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(value_float):
            return None
        values.append(value_float)
    return np.array(values, dtype=float)


def _pca_and_distances(target_properties: Dict[str, Any], candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    target_vector = _descriptor_vector(target_properties)
    vectors: List[np.ndarray] = []
    valid_candidates: List[Dict[str, Any]] = []
    if target_vector is None:
        return candidates
    vectors.append(target_vector)
    for candidate in candidates:
        vector = _descriptor_vector(candidate.get("properties"))
        if vector is None:
            continue
        vectors.append(vector)
        valid_candidates.append(candidate)
    if not valid_candidates:
        return candidates

    matrix = np.vstack(vectors)
    scale = matrix.max(axis=0) - matrix.min(axis=0)
    scale[scale == 0] = 1.0
    normalized = (matrix - matrix.min(axis=0)) / scale
    distances = np.linalg.norm(normalized[1:] - normalized[0], axis=1)

    centered = normalized - normalized.mean(axis=0, keepdims=True)
    if centered.shape[0] >= 2 and np.any(np.abs(centered) > 1e-12):
        u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
        coordinates = u[:, :2] * singular_values[:2]
    else:
        coordinates = np.zeros((centered.shape[0], 2), dtype=float)
    if coordinates.shape[1] == 1:
        coordinates = np.column_stack([coordinates[:, 0], np.zeros(coordinates.shape[0])])

    for index, candidate in enumerate(valid_candidates):
        candidate["global_distance"] = round(float(distances[index]), 6)
        candidate["pca"] = {
            "x": round(float(coordinates[index + 1, 0] - coordinates[0, 0]), 6),
            "y": round(float(coordinates[index + 1, 1] - coordinates[0, 1]), 6),
        }
    for candidate in candidates:
        if "global_distance" not in candidate:
            candidate["global_distance"] = None
            candidate["pca"] = {"x": 0.0, "y": 0.0}
    return candidates


def _build_points(target: Dict[str, Any], candidates: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    points = [{
        "id": "target",
        "label": "Molécula alvo",
        "name": "Molécula alvo",
        "is_target": True,
        "x": 0.0,
        "y": 0.0,
        "similarity": 1.0,
        "global_distance": 0.0,
    }]
    for index, candidate in enumerate(candidates, start=1):
        pca = candidate.get("pca") or {"x": 0.0, "y": 0.0}
        points.append({
            "id": f"cid-{candidate.get('cid', index)}",
            "label": candidate.get("name") or f"CID {candidate.get('cid', index)}",
            "name": candidate.get("name") or f"CID {candidate.get('cid', index)}",
            "is_target": False,
            "x": pca.get("x", 0.0),
            "y": pca.get("y", 0.0),
            "similarity": candidate.get("similarity"),
            "global_distance": candidate.get("global_distance"),
        })
    return points


def _property_rows(properties: Dict[str, Any]) -> Dict[str, Any]:
    return {key: properties.get(key) for key in DESCRIPTOR_KEYS}


def search_nitrosamine_space(
    smiles: str,
    *,
    threshold: int = DEFAULT_THRESHOLD,
    max_records: int = DEFAULT_MAX_RECORDS,
    max_candidates: int = DEFAULT_MAX_CANDIDATES,
) -> Dict[str, Any]:
    normalized = smiles.strip() if isinstance(smiles, str) else ""
    if not normalized:
        return _base_result(normalized, "invalid_input", "SMILES deve ser uma string não vazia.")
    try:
        target_mol = Chem.MolFromSmiles(normalized)
    except Exception:
        target_mol = None
    if target_mol is None:
        return _base_result(normalized, "invalid_smiles", "O SMILES alvo não pôde ser interpretado pelo RDKit.")

    target_canonical = _canonical(target_mol)
    target_properties = get_properties(target_mol) or {}
    target = {
        "smiles": normalized,
        "canonical_smiles": target_canonical,
        "svg": mol_to_svg(target_mol, size=320),
        "properties": _property_rows(target_properties),
    }
    threshold = max(1, min(int(threshold), 99))
    max_records = max(1, min(int(max_records), 100))
    max_candidates = max(1, min(int(max_candidates), 10))
    encoded = urllib_parse.quote(normalized, safe="")
    search_url = f"{PUBCHEM_BASE_URL}/compound/fastsimilarity_2d/smiles/{encoded}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}"
    search_payload, search_error = _get_json(search_url)
    if search_error:
        result = _base_result(normalized, "pubchem_unavailable", search_error["message"])
        result.update({"target": target, "search": {"threshold": threshold, "max_records": max_records, "retrieved_cids": 0, "n_nitroso_candidates": 0, "selected_candidates": 0}, "candidates": [], "points": []})
        return result

    cids = (((search_payload or {}).get("IdentifierList") or {}).get("CID") or [])
    cids = [int(cid) for cid in cids if str(cid).isdigit()][:max_records]
    if not cids:
        result = _base_result(normalized, "no_nitrosamines", "Nenhuma nitrosamina foi encontrada no lote amostrado do PubChem.")
        result.update({"target": target, "search": {"threshold": threshold, "max_records": max_records, "retrieved_cids": 0, "n_nitroso_candidates": 0, "selected_candidates": 0}, "candidates": [], "points": []})
        return result

    cid_path = ",".join(str(cid) for cid in cids)
    properties_url = (
        f"{PUBCHEM_BASE_URL}/compound/cid/{cid_path}/property/Title,IUPACName,SMILES,MolecularFormula,"
        "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
    )
    properties_payload, properties_error = _get_json(properties_url)
    if properties_error:
        result = _base_result(normalized, "pubchem_unavailable", properties_error["message"])
        result.update({"target": target, "search": {"threshold": threshold, "max_records": max_records, "retrieved_cids": len(cids), "n_nitroso_candidates": 0, "selected_candidates": 0}, "candidates": [], "points": []})
        return result

    target_fp = get_fingerprint(target_mol, SPACE_FINGERPRINT)
    candidates: List[Dict[str, Any]] = []
    raw_properties = ((properties_payload or {}).get("PropertyTable") or {}).get("Properties") or []
    for item in raw_properties:
        if not isinstance(item, dict):
            continue
        candidate_smiles = item.get("SMILES")
        if not isinstance(candidate_smiles, str) or not candidate_smiles.strip():
            continue
        try:
            candidate_mol = Chem.MolFromSmiles(candidate_smiles)
        except Exception:
            candidate_mol = None
        if candidate_mol is None or not _has_n_nitroso(candidate_mol):
            continue
        canonical = _canonical(candidate_mol)
        if canonical == target_canonical:
            continue
        candidate_fp = get_fingerprint(candidate_mol, SPACE_FINGERPRINT)
        if target_fp is None or candidate_fp is None:
            continue
        similarity = calc_similarity(target_fp, candidate_fp, "Tanimoto")
        properties = get_properties(candidate_mol) or {}
        name = item.get("Title") or item.get("IUPACName") or f"CID {item.get('CID')}"
        candidate = {
            "cid": item.get("CID"),
            "name": name,
            "iupac_name": item.get("IUPACName"),
            "smiles": candidate_smiles,
            "canonical_smiles": canonical,
            "is_n_nitroso": True,
            "fingerprint": SPACE_FINGERPRINT,
            "similarity": similarity,
            "classification": classify_similarity(similarity, "Tanimoto"),
            "svg": mol_to_svg(candidate_mol, size=220),
            "properties": _property_rows(properties),
            "pubchem_properties": {
                "molecular_formula": item.get("MolecularFormula"),
                "molecular_weight": item.get("MolecularWeight"),
                "xlogp": item.get("XLogP"),
                "tpsa": item.get("TPSA"),
                "hbd": item.get("HBondDonorCount"),
                "hba": item.get("HBondAcceptorCount"),
                "rotatable_bonds": item.get("RotatableBondCount"),
            },
        }
        candidates.append(candidate)

    candidates.sort(key=lambda item: (-float(item.get("similarity", 0.0)), int(item.get("cid") or 0)))
    candidates = candidates[:max_candidates]
    candidates = _pca_and_distances(target_properties, candidates)
    candidates.sort(key=lambda item: (float(item.get("global_distance") or 999999), -float(item.get("similarity", 0.0))))
    points = _build_points(target, candidates)
    result = _base_result(normalized, "ok" if candidates else "no_nitrosamines", "Busca e mapeamento de nitrosaminas concluídos." if candidates else "Nenhuma nitrosamina foi encontrada no lote amostrado do PubChem.")
    result.update({
        "target": target,
        "search": {
            "threshold": threshold,
            "max_records": max_records,
            "retrieved_cids": len(cids),
            "n_nitroso_candidates": len(candidates),
            "selected_candidates": len(candidates),
            "search_url": search_url,
            "fingerprint": SPACE_FINGERPRINT,
            "selection_method": "Filtro SMARTS [N;X3][N;X2]=O + MACCS/Tanimoto + distância euclidiana normalizada dos descritores",
        },
        "descriptor_keys": list(DESCRIPTOR_KEYS),
        "candidates": candidates,
        "points": points,
    })
    return result
