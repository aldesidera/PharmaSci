"""Build the fixed z-score profile used by Nitro.RA's EMA chemical space."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from analysis import get_fingerprint, get_properties
from chemo_suite.core.chemical_space import (
    calculate_multimodal_space,
    descriptor_vector,
    fit_zscore_profile,
    profile_to_json_dict,
)


def build_profile(snapshot_path: Path, output_path: Path) -> int:
    payload = json.loads(snapshot_path.read_text(encoding="utf-8"))
    records = []
    seen = set()
    for record in payload.get("records", []):
        if record.get("sheet") != "N-nitrosamines":
            continue
        canonical = record.get("canonical_smiles")
        if not canonical or canonical in seen:
            continue
        molecule = Chem.MolFromSmiles(canonical)
        vector = descriptor_vector(get_properties(molecule) if molecule else None)
        if molecule is None or vector is None:
            continue
        seen.add(canonical)
        records.append((canonical, vector.tolist()))

    if not records:
        raise RuntimeError("Nenhuma estrutura elegível foi encontrada na folha N-nitrosamines.")

    vectors = [vector for _, vector in records]
    fingerprints = [get_fingerprint(Chem.MolFromSmiles(canonical), "MACCS") for canonical, _ in records]
    profile = profile_to_json_dict(fit_zscore_profile(vectors))
    metrics = calculate_multimodal_space(fingerprints, vectors, "Tanimoto", profile=profile)
    upper = np.triu_indices(len(records), k=1)
    profile.update({
        "dataset": "EMA Appendix 1 chemical-space reference library",
        "reference_number": payload.get("reference_number"),
        "last_updated": payload.get("last_updated"),
        "source_url": payload.get("source_url"),
        "page_url": payload.get("page_url"),
        "sheet": "N-nitrosamines",
        "canonical_smiles": [canonical for canonical, _ in records],
        "deduplication": "canonical_smiles; first occurrence retained",
        "outlier_policy": "No automatic outlier removal; classical population z-score retained for auditability.",
        "fq_normalizer": round(float(np.max(metrics["physicochemical_distances"][upper])), 12) if len(records) > 1 else 1.0,
        "fq_normalizer_method": "maior distância FQ entre pares da biblioteca EMA após z-score fixo",
        "profile_version": "EMA-Appendix1-Rev13-NN-243-v1",
        "snapshot_record_count": len(payload.get("records", [])),
    })
    output_path.write_text(json.dumps(profile, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return len(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", type=Path, default=Path("chemo_suite/apps/nitro_ra/data/ema_appendix1.json"))
    parser.add_argument("--output", type=Path, default=Path("chemo_suite/apps/nitro_ra/data/ema_chemical_space_profile.json"))
    args = parser.parse_args()
    print(f"Perfil EMA gerado com {build_profile(args.snapshot, args.output)} estruturas.")


if __name__ == "__main__":
    main()
