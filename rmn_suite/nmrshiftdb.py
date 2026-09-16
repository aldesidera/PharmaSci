"""Provedor experimental/HOSE do NMRShiftDB2 para o PharmaSci.

A API pública do NMRShiftDB2 retorna um espectro medido quando disponível e
uma previsão baseada em HOSE quando não há registro medido. O módulo preserva
essa proveniência e acrescenta cache local persistente com fallback offline.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import tempfile
import time
from typing import Any
from urllib.parse import quote
import xml.etree.ElementTree as ET

import requests
from rdkit import Chem

BASE_URL = "https://nmrshiftdb.nmr.uni-koeln.de/NmrshiftdbServlet/nmrshiftdbaction"
NS = {"cml": "http://www.xml-cml.org/schema", "nmr": "http://www.nmrshiftdb.org/dict"}
_DEFAULT_CACHE_TTL = 900.0
_DEFAULT_STALE_MAX_AGE = 7 * 24 * 60 * 60
_CACHE: dict[tuple[str, str], tuple[float, "RmnPrediction"]] = {}


class RmnProviderError(RuntimeError):
    pass


@dataclass(frozen=True)
class RmnPeak:
    ppm: float
    intensity: float | None = None
    multiplicity: str | None = None
    atom_refs: list[str] | None = None


@dataclass(frozen=True)
class RmnPrediction:
    smiles_input: str
    smiles_canonical: str
    nucleus: str
    provenance: str  # measured | hose-predicted | unavailable
    source_url: str
    spectrum_id: str | None
    molecule_id: str | None
    solvent: str | None
    temperature_kelvin: float | None
    field_mhz: str | None
    peaks: list[RmnPeak]
    warnings: list[str]
    cache_status: str = "live"  # live | fresh-cache | stale-cache

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["peaks"] = [asdict(peak) for peak in self.peaks]
        return value


def _cache_path() -> Path:
    configured = os.environ.get("PHARMASCI_NMR_CACHE_PATH", "")
    if configured.strip():
        return Path(configured).expanduser()
    return Path(__file__).resolve().parents[2] / ".cache" / "nmrshiftdb.json"


def _cache_key(canonical: str, nucleus: str) -> str:
    return f"{canonical}|{nucleus}"


def _prediction_from_dict(value: dict[str, Any]) -> RmnPrediction:
    peaks = [RmnPeak(**peak) for peak in value.get("peaks", []) if isinstance(peak, dict)]
    return RmnPrediction(
        smiles_input=value.get("smiles_input", ""),
        smiles_canonical=value.get("smiles_canonical", ""),
        nucleus=value.get("nucleus", ""),
        provenance=value.get("provenance", "unavailable"),
        source_url=value.get("source_url", ""),
        spectrum_id=value.get("spectrum_id"),
        molecule_id=value.get("molecule_id"),
        solvent=value.get("solvent"),
        temperature_kelvin=value.get("temperature_kelvin"),
        field_mhz=value.get("field_mhz"),
        peaks=peaks,
        warnings=list(value.get("warnings", [])),
        cache_status=value.get("cache_status", "fresh-cache"),
    )


def _read_disk_cache() -> dict[str, dict[str, Any]]:
    path = _cache_path()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        entries = payload.get("entries", {}) if isinstance(payload, dict) else {}
        return entries if isinstance(entries, dict) else {}
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        return {}


def _write_disk_cache(entries: dict[str, dict[str, Any]]) -> None:
    path = _cache_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"version": 1, "updated_at": time.time(), "entries": entries}
    fd, temporary = tempfile.mkstemp(prefix="nmrshiftdb-", suffix=".json", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _store_cache(key: str, prediction: RmnPrediction, fetched_at: float) -> None:
    entries = _read_disk_cache()
    entries[key] = {"fetched_at": fetched_at, "prediction": prediction.to_dict()}
    _write_disk_cache(entries)


def _disk_entry(key: str) -> tuple[float, RmnPrediction] | None:
    entry = _read_disk_cache().get(key)
    if not isinstance(entry, dict) or not isinstance(entry.get("prediction"), dict):
        return None
    try:
        fetched_at = float(entry.get("fetched_at"))
    except (TypeError, ValueError):
        return None
    return fetched_at, _prediction_from_dict(entry["prediction"])


def canonicalize_smiles(smiles: str) -> str:
    value = (smiles or "").strip()
    if not value:
        raise RmnProviderError("Informe um SMILES.")
    mol = Chem.MolFromSmiles(value)
    if mol is None:
        raise RmnProviderError("SMILES inválido para o RDKit.")
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _text(element: ET.Element | None) -> str | None:
    if element is None or element.text is None:
        return None
    return element.text.strip() or None


def _parse_float(value: str | None) -> float | None:
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def parse_cml(xml_text: str, smiles_input: str, smiles_canonical: str, nucleus: str, source_url: str) -> RmnPrediction:
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        raise RmnProviderError(f"Resposta CML inválida: {exc}") from exc

    spectrum = root.find(".//cml:spectrum", NS)
    if spectrum is None:
        return RmnPrediction(smiles_input, smiles_canonical, nucleus, "unavailable", source_url, None, None, None, None, None, [], ["NMRShiftDB2 não retornou um espectro."])

    spectrum_id = spectrum.get("id")
    molecule_id = spectrum.get("moleculeRef")
    measured = bool(spectrum_id and spectrum_id.startswith("nmrshiftdb"))
    peaks: list[RmnPeak] = []
    for peak in spectrum.findall(".//cml:peak", NS):
        ppm = _parse_float(peak.get("xValue"))
        if ppm is None:
            continue
        refs = peak.get("atomRefs")
        peaks.append(RmnPeak(ppm, _parse_float(peak.get("peakHeight")), peak.get("peakMultiplicity"), refs.split() if refs else None))

    def scalar(dict_ref: str) -> str | None:
        node = spectrum.find(f".//cml:scalar[@dictRef='{dict_ref}']", NS)
        return _text(node)

    solvent = _text(spectrum.find(".//cml:substance[@role='subst:solvent']", NS))
    if solvent is None:
        substance = spectrum.find(".//cml:substance", NS)
        solvent = substance.get("title") if substance is not None else None
    temp = _parse_float(scalar("cml:temp"))
    field = scalar("cml:field")
    warnings = [] if measured else ["O NMRShiftDB2 não encontrou registro medido; o resultado é previsão HOSE."]
    return RmnPrediction(smiles_input, smiles_canonical, nucleus, "measured" if measured else "hose-predicted", source_url, spectrum_id, molecule_id, solvent, temp, field, peaks, warnings)


def clear_cache() -> None:
    _CACHE.clear()


def cache_size() -> int:
    return len(_CACHE)


def _configured_seconds(name: str, default: float) -> float:
    try:
        return max(0.0, float(os.environ.get(name, default)))
    except (TypeError, ValueError):
        return default


def predict_from_smiles(
    smiles: str,
    nucleus: str = "13C",
    timeout: float = 30.0,
    cache_ttl: float | None = None,
    stale_max_age: float | None = None,
) -> RmnPrediction:
    if cache_ttl is None:
        cache_ttl = _configured_seconds("PHARMASCI_NMR_CACHE_TTL_SECONDS", _DEFAULT_CACHE_TTL)
    if stale_max_age is None:
        stale_max_age = _configured_seconds("PHARMASCI_NMR_STALE_MAX_AGE_SECONDS", _DEFAULT_STALE_MAX_AGE)
    if nucleus not in {"1H", "13C", "15N", "19F", "31P", "11B", "29Si"}:
        raise RmnProviderError("Núcleo não suportado pelo endpoint NMRShiftDB2.")
    canonical = canonicalize_smiles(smiles)
    key = (canonical, nucleus)
    disk_key = _cache_key(canonical, nucleus)
    now = time.time()
    monotonic_now = time.monotonic()
    cached = _CACHE.get(key)
    if cached is None:
        disk = _disk_entry(disk_key)
        if disk is not None:
            fetched_at, prediction = disk
            _CACHE[key] = (monotonic_now + max(0.0, cache_ttl - (now - fetched_at)), prediction)
            cached = _CACHE[key]
    if cached is not None:
        expires_at, prediction = cached
        if cache_ttl > 0 and expires_at > monotonic_now:
            return RmnPrediction(**{**prediction.to_dict(), "peaks": prediction.peaks, "cache_status": "fresh-cache"})

    encoded = quote(canonical, safe="")
    url = f"{BASE_URL}/searchorpredict/smiles/{encoded}/spectrumtype/{quote(nucleus, safe='')}"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        prediction = parse_cml(response.text, smiles.strip(), canonical, nucleus, url)
    except (requests.RequestException, RmnProviderError) as exc:
        disk = _disk_entry(disk_key)
        if disk is not None and now - disk[0] <= stale_max_age:
            stale = disk[1]
            warnings = list(stale.warnings) + [f"NMRShiftDB2 indisponível; resultado reutilizado do cache local ({int(now - disk[0])} s)."]
            return RmnPrediction(**{**stale.to_dict(), "peaks": stale.peaks, "warnings": warnings, "cache_status": "stale-cache"})
        if isinstance(exc, RmnProviderError):
            raise
        raise RmnProviderError(f"Falha ao consultar NMRShiftDB2: {exc}") from exc

    _CACHE[key] = (time.monotonic() + max(0.0, cache_ttl), prediction)
    if cache_ttl > 0:
        _store_cache(disk_key, prediction, now)
    return prediction
