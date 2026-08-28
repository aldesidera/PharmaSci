"""Nitro.RA chemical-space search and EMA reference-space mapping."""

from __future__ import annotations

import json
import logging
import os
import time
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
from urllib import error as urllib_error
from urllib import parse as urllib_parse
from urllib import request as urllib_request

import numpy as np
from rdkit import Chem

from analysis import classify_similarity, get_fingerprint, get_properties, mol_to_svg
from chemo_suite.core.chemical_space import (
    DESCRIPTOR_KEYS,
    calculate_multimodal_space,
    classical_mds,
    descriptor_vector,
    normalized_stress,
    select_nearest_indices,
)
from .cpca import EMA_APPENDIX_UPDATED, EMA_APPENDIX_URL, EMA_APPENDIX_VERSION, _ema_index

logger = logging.getLogger(__name__)

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_SOURCE_URL = "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest"
EMA_PROFILE_PATH = Path(__file__).resolve().parent / "data" / "ema_chemical_space_profile.json"
DEFAULT_THRESHOLD = int(os.getenv("MOLSIM_PUBCHEM_SPACE_THRESHOLD", "50"))
DEFAULT_MAX_RECORDS = min(int(os.getenv("MOLSIM_PUBCHEM_SPACE_MAX_RECORDS", "40")), 40)
DEFAULT_MAX_CANDIDATES = int(os.getenv("MOLSIM_PUBCHEM_SPACE_MAX_CANDIDATES", "10"))
DEFAULT_TIMEOUT = float(os.getenv("MOLSIM_PUBCHEM_SPACE_TIMEOUT", "12"))
DEFAULT_CACHE_TTL = float(os.getenv("MOLSIM_PUBCHEM_SPACE_CACHE_TTL_SECONDS", "3600"))
SPACE_FINGERPRINT = "MACCS"
SPACE_METRIC = "Tanimoto"
PUBCHEM_BATCH_LIMIT = 40
DISPLAY_LIMIT = 10

N_NITROSO_PATTERN = Chem.MolFromSmarts("[N;X3][N;X2]=O")
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
            "A similaridade 2D, os descritores e a projeção MDS são ferramentas de triagem e não substituem avaliação toxicológica ou regulatória.",
        ],
    }


def _canonical(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _has_n_nitroso(mol: Chem.Mol) -> bool:
    return bool(N_NITROSO_PATTERN and mol.HasSubstructMatch(N_NITROSO_PATTERN))


def _descriptor_vector(properties: Optional[Dict[str, Any]]) -> Optional[np.ndarray]:
    return descriptor_vector(properties)


def _property_rows(properties: Dict[str, Any]) -> Dict[str, Any]:
    return {key: properties.get(key) for key in DESCRIPTOR_KEYS}


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
            "source": candidate.get("source"),
        })
    return points


def _rank_and_project(
    target_mol: Chem.Mol,
    target_properties: Dict[str, Any],
    candidates: List[Dict[str, Any]],
    candidate_molecules: Sequence[Chem.Mol],
    *,
    profile: Optional[Dict[str, Any]] = None,
    fq_normalizer: Optional[float] = None,
    display_limit: int = DISPLAY_LIMIT,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]:
    target_fp = get_fingerprint(target_mol, SPACE_FINGERPRINT)
    target_vector = _descriptor_vector(target_properties)
    if target_fp is None or target_vector is None:
        return [], _build_points({}, []), {"scored_candidates": 0, "fq_normalizer": None, "mds_stress": 0.0}

    fingerprints = [target_fp]
    vectors = [target_vector]
    valid_candidates: List[Dict[str, Any]] = []
    valid_molecules: List[Chem.Mol] = []
    for candidate, molecule in zip(candidates, candidate_molecules):
        vector = _descriptor_vector(candidate.get("properties"))
        fingerprint = get_fingerprint(molecule, SPACE_FINGERPRINT)
        if vector is None or fingerprint is None:
            continue
        fingerprints.append(fingerprint)
        vectors.append(vector)
        valid_candidates.append(candidate)
        valid_molecules.append(molecule)

    if not valid_candidates:
        return [], _build_points({}, []), {"scored_candidates": 0, "fq_normalizer": None, "mds_stress": 0.0}

    metrics = calculate_multimodal_space(
        fingerprints,
        vectors,
        SPACE_METRIC,
        profile=profile,
        fq_normalizer=fq_normalizer,
    )
    selected_indices = select_nearest_indices(
        metrics["reference_global_distances"],
        display_limit=min(max(0, int(display_limit)), DISPLAY_LIMIT),
        reference_index=0,
    )
    selected_distances = metrics["global_distances"][np.ix_(selected_indices, selected_indices)]
    coordinates = classical_mds(selected_distances)
    coordinates = coordinates - coordinates[0]
    selected_candidates: List[Dict[str, Any]] = []
    for local_index, original_index in enumerate(selected_indices[1:], start=1):
        candidate_index = original_index - 1
        candidate = valid_candidates[candidate_index]
        candidate["similarity"] = round(float(1.0 - metrics["reference_structural_distances"][original_index]), 6)
        candidate["classification"] = classify_similarity(candidate["similarity"], SPACE_METRIC)
        candidate["physicochemical_distance"] = round(float(metrics["normalized_physicochemical_distances"][0, original_index]), 6)
        candidate["global_distance"] = round(float(metrics["reference_global_distances"][original_index]), 6)
        candidate["global_similarity"] = round(float(1.0 - metrics["reference_global_distances"][original_index]), 6)
        candidate["pca"] = {
            "x": round(float(coordinates[local_index, 0]), 6),
            "y": round(float(coordinates[local_index, 1]), 6),
        }
        if not candidate.get("svg"):
            candidate["svg"] = mol_to_svg(valid_molecules[candidate_index], size=220)
        selected_candidates.append(candidate)

    # A maior similaridade global equivale à menor distância global; manter essa
    # ordem explícita garante que tabelas e gráficos iniciem pelo vizinho mais similar.
    selected_candidates.sort(key=lambda item: (-float(item.get("global_similarity", 0.0)), float(item.get("global_distance", 1.0))))
    target = {"properties": _property_rows(target_properties)}
    points = _build_points(target, selected_candidates)
    return selected_candidates, points, {
        "scored_candidates": len(valid_candidates),
        "fq_normalizer": round(float(metrics["fq_normalizer"]), 6),
        "mds_stress": round(float(normalized_stress(selected_distances, coordinates)), 6),
    }


@lru_cache(maxsize=1)
def _ema_profile() -> Optional[Dict[str, Any]]:
    try:
        payload = json.loads(EMA_PROFILE_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if payload.get("sheet") != "N-nitrosamines":
        return None
    return payload


def _ema_candidates(target_canonical: Optional[str] = None) -> Tuple[Optional[Dict[str, Any]], List[Dict[str, Any]], List[Chem.Mol], bool]:
    profile = _ema_profile()
    if not profile:
        return None, [], [], False
    index = _ema_index()
    candidates: List[Dict[str, Any]] = []
    molecules: List[Chem.Mol] = []
    target_excluded = False
    for canonical in profile.get("canonical_smiles", []):
        record = index.get(canonical)
        if not record:
            continue
        molecule = Chem.MolFromSmiles(canonical)
        if molecule is None:
            continue
        if target_canonical and _canonical(molecule) == target_canonical:
            target_excluded = True
            continue
        properties = get_properties(molecule) or {}
        if _descriptor_vector(properties) is None:
            continue
        candidates.append({
            "name": record.get("name") or record.get("iupac_name") or canonical,
            "smiles": record.get("smiles") or canonical,
            "canonical_smiles": canonical,
            "source": "EMA Appendix 1",
            "sheet": record.get("sheet"),
            "cas_rn": record.get("cas_rn"),
            "ai_ng_day": record.get("ai_ng_day"),
            "ai_ng_day_raw": record.get("ai_ng_day_raw"),
            "cpca_category": record.get("cpca_category"),
            "note": record.get("note"),
            "svg": mol_to_svg(molecule, size=220),
            "properties": _property_rows(properties),
            "fingerprint": SPACE_FINGERPRINT,
        })
        molecules.append(molecule)
    return profile, candidates, molecules, target_excluded


def _ema_space(
    target_mol: Chem.Mol,
    target_properties: Dict[str, Any],
    *,
    display_limit: int = DISPLAY_LIMIT,
) -> Dict[str, Any]:
    target_canonical = _canonical(target_mol)
    profile, candidates, molecules, target_excluded = _ema_candidates(target_canonical)
    if not profile:
        return {
            "module": "ema_chemical_space",
            "status": "ema_reference_unavailable",
            "message": "O perfil local da biblioteca EMA não está disponível.",
            "candidates": [],
            "points": [],
            "search": {"library_size": 0, "source_library_size": 0, "selected_candidates": 0, "target_excluded": False},
            "warnings": ["A biblioteca EMA local não foi carregada; nenhum valor foi inferido."],
        }
    fixed_normalizer = profile.get("fq_normalizer")
    selected, points, diagnostics = _rank_and_project(
        target_mol,
        target_properties,
        candidates,
        molecules,
        profile=profile,
        fq_normalizer=fixed_normalizer,
        display_limit=display_limit,
    )
    message = "Mapeamento EMA concluído." if selected else "Nenhuma estrutura EMA elegível foi localizada no perfil."
    return {
        "module": "ema_chemical_space",
        "status": "ok" if selected else "no_ema_candidates",
        "message": message,
        "source": "EMA Appendix 1",
        "source_url": profile.get("source_url") or EMA_APPENDIX_URL,
        "reference_number": profile.get("reference_number") or EMA_APPENDIX_VERSION,
        "last_updated": profile.get("last_updated") or EMA_APPENDIX_UPDATED,
        "sheet": profile.get("sheet"),
        "target": {"properties": _property_rows(target_properties)},
        "search": {
            "library_size": len(candidates),
            "source_library_size": len(profile.get("canonical_smiles", [])),
            "target_excluded": target_excluded,
            "scored_candidates": diagnostics["scored_candidates"],
            "selected_candidates": len(selected),
            "selection_method": "MACCS/Tanimoto + distância multimodal quadrática com z-score fixo da biblioteca EMA; MDS clássico; exibição das 10 menores distâncias",
            "profile_n_structures": profile.get("n_structures"),
            "fq_normalizer": diagnostics["fq_normalizer"],
            "mds_stress": diagnostics["mds_stress"],
        },
        "descriptor_keys": list(DESCRIPTOR_KEYS),
        "fingerprint": SPACE_FINGERPRINT,
        "candidates": selected,
        "points": points,
        "warnings": [
            "A biblioteca EMA é um conjunto de referência regulatória e não representa o universo completo de nitrosaminas.",
            "A proximidade química não transfere automaticamente AI, cPCA ou conclusão regulatória ao alvo.",
            "O perfil z-score é fixo, versionado e não remove outliers automaticamente.",
        ],
    }


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
    max_records = max(1, min(int(max_records), PUBCHEM_BATCH_LIMIT))
    max_candidates = max(1, min(int(max_candidates), DISPLAY_LIMIT))
    encoded = urllib_parse.quote(normalized, safe="")
    search_url = f"{PUBCHEM_BASE_URL}/compound/fastsimilarity_2d/smiles/{encoded}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}"
    search_payload, search_error = _get_json(search_url)
    if search_error:
        result = _base_result(normalized, "pubchem_unavailable", search_error["message"])
        result.update({
            "target": target,
            "search": {"threshold": threshold, "max_records": max_records, "retrieved_cids": 0, "n_nitroso_candidates": 0, "selected_candidates": 0, "target_excluded": False},
            "candidates": [],
            "points": [],
            "ema_space": _ema_space(target_mol, target_properties, display_limit=DISPLAY_LIMIT),
        })
        return result

    raw_cids = (((search_payload or {}).get("IdentifierList") or {}).get("CID") or [])
    cids = []
    for cid in raw_cids[:max_records]:
        try:
            cids.append(int(cid))
        except (TypeError, ValueError):
            continue
    if not cids:
        result = _base_result(normalized, "no_nitrosamines", "Nenhuma nitrosamina foi encontrada no lote amostrado do PubChem.")
        result.update({
            "target": target,
            "search": {"threshold": threshold, "max_records": max_records, "retrieved_cids": 0, "n_nitroso_candidates": 0, "selected_candidates": 0, "target_excluded": False},
            "candidates": [],
            "points": _build_points(target, []),
            "ema_space": _ema_space(target_mol, target_properties, display_limit=DISPLAY_LIMIT),
        })
        return result

    cid_path = ",".join(str(cid) for cid in cids)
    properties_url = (
        f"{PUBCHEM_BASE_URL}/compound/cid/{cid_path}/property/Title,IUPACName,SMILES,MolecularFormula,"
        "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
    )
    properties_payload, properties_error = _get_json(properties_url)
    if properties_error:
        result = _base_result(normalized, "pubchem_unavailable", properties_error["message"])
        result.update({
            "target": target,
            "search": {"threshold": threshold, "max_records": max_records, "retrieved_cids": len(cids), "n_nitroso_candidates": 0, "selected_candidates": 0, "target_excluded": False},
            "candidates": [],
            "points": [],
            "ema_space": _ema_space(target_mol, target_properties, display_limit=DISPLAY_LIMIT),
        })
        return result

    candidates: List[Dict[str, Any]] = []
    candidate_molecules: List[Chem.Mol] = []
    target_excluded = False
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
            target_excluded = True
            continue
        properties = get_properties(candidate_mol) or {}
        if _descriptor_vector(properties) is None:
            continue
        name = item.get("Title") or item.get("IUPACName") or f"CID {item.get('CID')}"
        candidates.append({
            "cid": item.get("CID"),
            "name": name,
            "iupac_name": item.get("IUPACName"),
            "smiles": candidate_smiles,
            "canonical_smiles": canonical,
            "is_n_nitroso": True,
            "fingerprint": SPACE_FINGERPRINT,
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
        })
        candidate_molecules.append(candidate_mol)

    selected, points, diagnostics = _rank_and_project(
        target_mol,
        target_properties,
        candidates,
        candidate_molecules,
        display_limit=max_candidates,
    )
    ema_space = _ema_space(target_mol, target_properties, display_limit=DISPLAY_LIMIT)
    result = _base_result(
        normalized,
        "ok" if selected else "no_nitrosamines",
        "Busca, cálculo multimodal e mapeamento PubChem concluídos." if selected else "Nenhuma nitrosamina elegível foi encontrada no lote amostrado do PubChem.",
    )
    result.update({
        "target": target,
        "search": {
            "threshold": threshold,
            "max_records": max_records,
            "retrieved_cids": len(cids),
            "n_nitroso_candidates": len(candidates),
            "target_excluded": target_excluded,
            "scored_candidates": diagnostics["scored_candidates"],
            "selected_candidates": len(selected),
            "search_url": search_url,
            "fingerprint": SPACE_FINGERPRINT,
            "metric": SPACE_METRIC,
            "selection_method": "Filtro SMARTS [N;X3][N;X2]=O + MACCS/Tanimoto + distância multimodal quadrática em MW, LogP, TPSA, HBD, HBA e RotB; lote PubChem limitado a 40 CIDs, seleção das 10 menores distâncias e MDS clássico",
            "fq_normalizer": diagnostics["fq_normalizer"],
            "mds_stress": diagnostics["mds_stress"],
            "display_limit": DISPLAY_LIMIT,
        },
        "descriptor_keys": list(DESCRIPTOR_KEYS),
        "candidates": selected,
        "points": points if selected else _build_points(target, []),
        "ema_space": ema_space,
    })
    return result
