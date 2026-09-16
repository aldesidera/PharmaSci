"""Modelos transparentes e contrato de previsão de solubilidade do Mol.FQ.

ESOL e GSE são baselines estimativos. O modelo calibrado só é usado quando um
artefato versionado e aprovado é explicitamente configurado.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors


def _finite(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def esol(mol: Any) -> Dict[str, Any]:
    """Calcula o baseline ESOL de Delaney em log10(mol/L).

    A implementação usa os descritores RDKit e mantém a fórmula explícita no
    retorno para que o resultado seja auditável e não confundido com S0 medida.
    """
    mw = float(Descriptors.MolWt(mol))
    logp = float(Crippen.MolLogP(mol))
    rotb = float(Lipinski.NumRotatableBonds(mol))
    heavy = max(mol.GetNumHeavyAtoms(), 1)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    aromatic_proportion = aromatic_atoms / heavy
    log_s = 0.16 - 1.5 * logp - 0.0062 * mw + 0.066 * rotb + 0.066 * aromatic_proportion
    return {
        "value": round(log_s, 6),
        "unit": "log10(mol/L)",
        "method": "ESOL (Delaney) + RDKit",
        "status": "baseline",
        "descriptors": {
            "logp": round(logp, 6),
            "molecular_weight": round(mw, 6),
            "rotatable_bonds": int(rotb),
            "aromatic_proportion": round(aromatic_proportion, 6),
        },
        "formula": "0.16 − 1.5·LogP − 0.0062·MW + 0.066·RotB + 0.066·AP",
    }


def gse(mol: Any, melting_point_c: Optional[float] = None) -> Dict[str, Any]:
    """Calcula GSE apenas quando o ponto de fusão foi fornecido."""
    logp = float(Crippen.MolLogP(mol))
    mp = _finite(melting_point_c)
    if mp is None:
        return {
            "value": None,
            "unit": "log10(mol/L)",
            "method": "GSE — General Solubility Equation",
            "status": "unavailable",
            "reason": "Ponto de fusão não informado; GSE não foi calculada.",
        }
    log_s = 0.5 - 0.01 * (mp - 25.0) - logp
    return {
        "value": round(log_s, 6),
        "unit": "log10(mol/L)",
        "method": "GSE — General Solubility Equation",
        "status": "baseline",
        "melting_point_c": round(mp, 3),
        "formula": "logS = 0.5 − 0.01·(Tm − 25) − LogP",
    }


def _ionization_factor(ph: float, pka_acidic: Optional[float], pka_basic: Optional[float]) -> float:
    factor = 1.0
    if pka_acidic is not None:
        factor += 10.0 ** (ph - pka_acidic)
    if pka_basic is not None:
        factor += 10.0 ** (pka_basic - ph)
    return factor


def ph_curve(log_s0: Optional[float], pka_acidic: Optional[float], pka_basic: Optional[float], ph_values: Iterable[float] = range(0, 15)) -> Dict[str, Any]:
    """Gera S(pH) aparente a partir de S0 e pKa indicativos.

    Para anfólitos ou múltiplos centros, a soma é uma aproximação de primeira
    ordem e deve permanecer marcada como estimativa até existir um modelo de
    microespécies calibrado.
    """
    points = []
    if log_s0 is not None:
        for ph_value in ph_values:
            ph = float(ph_value)
            factor = _ionization_factor(ph, pka_acidic, pka_basic)
            points.append({
                "pH": round(ph, 3),
                "logS": round(float(log_s0) + math.log10(factor), 6),
                "ionization_factor": round(factor, 6),
            })
    svg = ""
    if points:
        values = [point["logS"] for point in points]
        low, high = min(values), max(values)
        span = max(high - low, 0.001)
        coords = []
        for index, point in enumerate(points):
            x = 42 + (index / max(len(points) - 1, 1)) * 518
            y = 154 - ((point["logS"] - low) / span) * 120
            coords.append((x, y))
        path = " ".join(f"{x:.2f},{y:.2f}" for x, y in coords)
        dots = "".join(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.2" fill="#166534"/>' for x, y in coords)
        svg = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 600 190" role="img" aria-label="Curva de solubilidade aparente em função do pH"><rect x="0" y="0" width="600" height="190" rx="10" fill="#f8fffa"/><line x1="42" y1="154" x2="560" y2="154" stroke="#94a3b8" stroke-width="1"/><line x1="42" y1="28" x2="42" y2="154" stroke="#94a3b8" stroke-width="1"/><polyline points="{path}" fill="none" stroke="#166534" stroke-width="2.5"/>{dots}<text x="42" y="177" font-size="10" fill="#64748b">pH 0</text><text x="292" y="177" font-size="10" fill="#64748b">pH 7</text><text x="535" y="177" font-size="10" fill="#64748b">pH 14</text><text x="8" y="34" font-size="9" fill="#64748b">logS</text></svg>'
    return {
        "points": points,
        "log_s0": log_s0,
        "pka_acidic": pka_acidic,
        "pka_basic": pka_basic,
        "svg": svg,
        "method": "Equação de ionização monoprotônica combinada; estimativa",
        "status": "estimated" if points else "unavailable",
        "warning": "A curva assume pKa indicativos e não substitui um modelo quantitativo de microespécies.",
    }


def _hybrid_features(mol: Any, geometry: Optional[Dict[str, Any]] = None) -> list[float]:
    """Retorna features auditáveis para o baseline e para artefatos híbridos."""
    features = [
        float(Descriptors.MolWt(mol)),
        float(Crippen.MolLogP(mol)),
        float(Descriptors.TPSA(mol)),
        float(Lipinski.NumRotatableBonds(mol)),
    ]
    summary = (geometry or {}).get("distance_summary", {})
    features.extend([
        float((geometry or {}).get("atom_count_3d", 0.0)),
        float(summary.get("min_angstrom") or 0.0),
        float(summary.get("mean_angstrom") or 0.0),
        float(summary.get("max_angstrom") or 0.0),
    ])
    return features


def calibrated_prediction(mol: Any, canonical_smiles: str, artifact_path: Optional[str] = None, geometry: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Carrega um artefato calibrado aprovado, se configurado.

    O artefato deve conter um manifesto JSON ao lado do modelo e expor
    ``predict(features)`` ou ``predict_from_smiles(smiles)``. Sem artefato,
    retorna indisponível e nunca inventa uma previsão calibrada.
    """
    path_value = artifact_path or os.getenv("PHARMASCI_SOLUBILITY_ARTIFACT")
    if not path_value:
        return {"value": None, "unit": "log10(mol/L)", "status": "unavailable", "method": "Modelo calibrado Mol.FQ", "reason": "Nenhum artefato aprovado foi configurado."}
    path = Path(path_value)
    manifest_path = path.with_suffix(".json")
    if not path.exists():
        return {"value": None, "unit": "log10(mol/L)", "status": "unavailable", "method": "Modelo calibrado Mol.FQ", "reason": f"Artefato não encontrado: {path}"}
    try:
        import joblib
        model = joblib.load(path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
        if manifest.get("status") != "approved":
            return {"value": None, "unit": "log10(mol/L)", "status": "unavailable", "method": "Modelo calibrado Mol.FQ", "reason": "O manifesto não está aprovado."}
        feature_schema = manifest.get("feature_schema", "baseline_v1")
        if hasattr(model, "predict_from_smiles"):
            predicted = model.predict_from_smiles(canonical_smiles)
        elif hasattr(model, "predict"):
            if feature_schema == "hybrid_3d_v1" and geometry:
                features = [_hybrid_features(mol, geometry)]
            else:
                features = [_hybrid_features(mol)[:4]]
            predicted = model.predict(features)
        else:
            raise TypeError("Artefato não expõe predict ou predict_from_smiles")
        value = float(predicted[0] if hasattr(predicted, "__len__") else predicted)
        return {"value": round(value, 6), "unit": "log10(mol/L)", "status": "calibrated", "method": manifest.get("model_name", "Modelo calibrado Mol.FQ"), "model_version": manifest.get("version"), "feature_schema": feature_schema, "metrics": manifest.get("metrics"), "domain": manifest.get("domain")}
    except Exception as exc:
        return {"value": None, "unit": "log10(mol/L)", "status": "unavailable", "method": "Modelo calibrado Mol.FQ", "reason": f"Falha ao carregar o artefato: {exc}"}


def build_solubility(mol: Any, canonical_smiles: str, pka_acidic: Any = None, pka_basic: Any = None, melting_point_c: Any = None, geometry: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    acidic = _finite(pka_acidic)
    basic = _finite(pka_basic)
    esol_result = esol(mol)
    gse_result = gse(mol, melting_point_c)
    calibrated = calibrated_prediction(mol, canonical_smiles, geometry=geometry)
    s0 = calibrated.get("value") if calibrated.get("status") == "calibrated" else esol_result.get("value")
    curve = ph_curve(s0, acidic, basic)
    return {
        "status": "estimated" if calibrated.get("status") != "calibrated" else "calibrated",
        "method": "Modelo calibrado Mol.FQ + ESOL/GSE + curva de ionização" if calibrated.get("status") == "calibrated" else "ESOL baseline + GSE opcional + curva de ionização",
        "s0": {"value": s0, "unit": "log10(mol/L)", "method": calibrated.get("method") if calibrated.get("status") == "calibrated" else esol_result["method"], "status": calibrated.get("status") if calibrated.get("status") == "calibrated" else "baseline"},
        "baselines": {"esol": esol_result, "gse": gse_result},
        "calibrated": calibrated,
        "ph_curve": curve,
        "ph7": next((point for point in curve["points"] if point["pH"] == 7.0), None),
        "domain": {"status": "not_evaluated", "message": "Domínio de aplicabilidade será preenchido pelo artefato calibrado aprovado."},
    }
