"""Shared, deterministic chemical-space calculations for Mol.Sim and Nitro.RA."""

from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


DESCRIPTOR_KEYS: Tuple[str, ...] = (
    "Massa Molecular (g/mol)",
    "Coeficiente de Partição (LogP)",
    "Área de Superfície Polar (Å²)",
    "Doadores de H (HBD)",
    "Receptores de H (HBA)",
    "Ligações rotacionáveis (RotB)",
)
DESCRIPTOR_LABELS: Tuple[str, ...] = ("MW", "LogP", "TPSA", "HBD", "HBA", "RotB")
STRUCTURAL_WEIGHT = 0.6
PHYSICOCHEMICAL_WEIGHT = 0.4
DEFAULT_DISPLAY_LIMIT = 10


def descriptor_vector(properties: Optional[Dict[str, Any]]) -> Optional[np.ndarray]:
    """Return the six-descriptor vector or ``None`` when a value is unusable."""
    if not properties:
        return None
    values: List[float] = []
    for key in DESCRIPTOR_KEYS:
        value = properties.get(key)
        if isinstance(value, bool):
            return None
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(number):
            return None
        values.append(number)
    return np.asarray(values, dtype=float)


def _coerce_matrix(matrix: Sequence[Sequence[float]]) -> np.ndarray:
    values = np.asarray(matrix, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(DESCRIPTOR_KEYS):
        raise ValueError(f"A matriz deve ter {len(DESCRIPTOR_KEYS)} colunas de descritores.")
    if values.shape[0] == 0:
        raise ValueError("A matriz de descritores não pode ser vazia.")
    return values.copy()


def fit_zscore_profile(fit_matrix: Sequence[Sequence[float]]) -> Dict[str, Any]:
    """Fit a population z-score profile, using medians only to fill missing values."""
    matrix = _coerce_matrix(fit_matrix)
    for column in range(matrix.shape[1]):
        values = matrix[:, column]
        finite = np.isfinite(values)
        fill = float(np.median(values[finite])) if finite.any() else 0.0
        values[~finite] = fill
    center = matrix.mean(axis=0)
    scale = matrix.std(axis=0)
    scale[scale == 0] = 1.0
    return {
        "center": [float(value) for value in center],
        "scale": [float(value) for value in scale],
        "n_structures": int(matrix.shape[0]),
        "method": "z-score populacional; valores ausentes preenchidos pela mediana do perfil",
    }


def apply_zscore_profile(matrix: Sequence[Sequence[float]], profile: Dict[str, Any]) -> np.ndarray:
    values = _coerce_matrix(matrix)
    center = np.asarray(profile.get("center", []), dtype=float)
    scale = np.asarray(profile.get("scale", []), dtype=float)
    if center.shape != (len(DESCRIPTOR_KEYS),) or scale.shape != (len(DESCRIPTOR_KEYS),):
        raise ValueError("Perfil z-score incompatível com os seis descritores.")
    scale = scale.copy()
    scale[~np.isfinite(scale) | (scale == 0)] = 1.0
    for column in range(values.shape[1]):
        finite = np.isfinite(values[:, column])
        fill = center[column]
        values[~finite, column] = fill
    return (values - center) / scale


def _pairwise_euclidean(values: np.ndarray) -> np.ndarray:
    differences = values[:, None, :] - values[None, :, :]
    distances = np.linalg.norm(differences, axis=2)
    return (distances + distances.T) / 2.0


def _pairwise_structural_distance(fingerprints: Sequence[Any], metric: str) -> np.ndarray:
    from analysis import calc_similarity

    count = len(fingerprints)
    distances = np.zeros((count, count), dtype=float)
    for left in range(count):
        for right in range(left + 1, count):
            similarity = float(calc_similarity(fingerprints[left], fingerprints[right], metric))
            structural = max(0.0, min(1.0, 1.0 - similarity))
            distances[left, right] = distances[right, left] = structural
    return distances


def _resolve_fq_normalizer(
    fq_distances: np.ndarray,
    *,
    reference_index: int,
    fq_normalizer: Optional[float],
) -> float:
    if fq_normalizer is not None:
        value = float(fq_normalizer)
    else:
        reference_distances = np.delete(fq_distances[reference_index], reference_index)
        value = float(np.max(reference_distances)) if reference_distances.size else 0.0
    return value if math.isfinite(value) and value > 0 else 1.0


def calculate_multimodal_space(
    fingerprints: Sequence[Any],
    descriptor_matrix: Sequence[Sequence[float]],
    metric: str = "Tanimoto",
    *,
    fit_matrix: Optional[Sequence[Sequence[float]]] = None,
    profile: Optional[Dict[str, Any]] = None,
    fq_normalizer: Optional[float] = None,
    reference_index: int = 0,
) -> Dict[str, Any]:
    """Calculate structural, physicochemical and quadratic global distances."""
    if not fingerprints or len(fingerprints) != len(descriptor_matrix):
        raise ValueError("Fingerprints e descritores devem ter o mesmo número de entradas.")
    matrix = _coerce_matrix(descriptor_matrix)
    if profile is None:
        profile = fit_zscore_profile(fit_matrix if fit_matrix is not None else matrix)
    standardized = apply_zscore_profile(matrix, profile)
    fq_distances = _pairwise_euclidean(standardized)
    structural_distances = _pairwise_structural_distance(fingerprints, metric)
    normalizer = _resolve_fq_normalizer(
        fq_distances,
        reference_index=reference_index,
        fq_normalizer=fq_normalizer,
    )
    normalized_fq = np.clip(fq_distances / normalizer, 0.0, 1.0)
    global_distances = np.sqrt(
        STRUCTURAL_WEIGHT * np.square(structural_distances)
        + PHYSICOCHEMICAL_WEIGHT * np.square(normalized_fq)
    )
    np.fill_diagonal(global_distances, 0.0)
    return {
        "profile": profile,
        "standardized": standardized,
        "structural_distances": structural_distances,
        "physicochemical_distances": fq_distances,
        "normalized_physicochemical_distances": normalized_fq,
        "global_distances": global_distances,
        "fq_normalizer": float(normalizer),
        "reference_global_distances": global_distances[reference_index].copy(),
        "reference_structural_distances": structural_distances[reference_index].copy(),
    }


def select_nearest_indices(
    reference_distances: Sequence[float],
    *,
    display_limit: int = DEFAULT_DISPLAY_LIMIT,
    reference_index: int = 0,
) -> List[int]:
    """Return the reference and at most ``display_limit`` nearest other entries."""
    limit = max(0, int(display_limit))
    candidates = [index for index in range(len(reference_distances)) if index != reference_index]
    candidates.sort(key=lambda index: (float(reference_distances[index]), index))
    return [reference_index] + candidates[:limit]


def classical_mds(distances: Sequence[Sequence[float]]) -> np.ndarray:
    """Return deterministic two-dimensional classical-MDS coordinates."""
    matrix = np.asarray(distances, dtype=float)
    count = matrix.shape[0]
    if count == 0:
        return np.zeros((0, 2), dtype=float)
    if count == 1:
        return np.zeros((1, 2), dtype=float)
    matrix = np.maximum((matrix + matrix.T) / 2.0, 0.0)
    np.fill_diagonal(matrix, 0.0)
    centering = np.eye(count) - np.ones((count, count)) / count
    gram = -0.5 * centering @ np.square(matrix) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh((gram + gram.T) / 2.0)
    order = np.argsort(eigenvalues)[::-1]
    positive = np.maximum(eigenvalues[order[:2]], 0.0)
    vectors = eigenvectors[:, order[:2]]
    coordinates = vectors * np.sqrt(positive)
    if coordinates.shape[1] == 1:
        coordinates = np.column_stack([coordinates[:, 0], np.zeros(count)])
    elif coordinates.shape[1] == 0:
        coordinates = np.zeros((count, 2))
    for axis in range(coordinates.shape[1]):
        anchor = coordinates[0, axis]
        if anchor < -1e-12:
            coordinates[:, axis] *= -1.0
        elif abs(anchor) <= 1e-12:
            nonzero = np.flatnonzero(np.abs(coordinates[:, axis]) > 1e-12)
            if nonzero.size and coordinates[nonzero[0], axis] < 0:
                coordinates[:, axis] *= -1.0
    return coordinates


def normalized_stress(distances: Sequence[Sequence[float]], coordinates: Sequence[Sequence[float]]) -> float:
    """Calculate normalized raw stress for the displayed MDS configuration."""
    target = np.asarray(distances, dtype=float)
    coords = np.asarray(coordinates, dtype=float)
    if target.shape[0] < 2 or coords.shape[0] != target.shape[0]:
        return 0.0
    embedded = _pairwise_euclidean(coords)
    upper = np.triu_indices(target.shape[0], k=1)
    denominator = float(np.sum(np.square(target[upper])))
    if denominator <= 0:
        return 0.0
    return float(np.sqrt(np.sum(np.square(target[upper] - embedded[upper])) / denominator))


def profile_to_json_dict(profile: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "descriptor_keys": list(DESCRIPTOR_KEYS),
        "descriptor_labels": list(DESCRIPTOR_LABELS),
        "center": [round(float(value), 12) for value in profile["center"]],
        "scale": [round(float(value), 12) for value in profile["scale"]],
        "n_structures": int(profile["n_structures"]),
        "method": profile["method"],
    }
