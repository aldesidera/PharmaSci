from pathlib import Path

from chemo_suite.apps.nitro_ra.cpca import calculate_cpca

output = Path("artifacts/cpca_fluxograma_nitrora_palette.svg")
result = calculate_cpca("CCN(CC)N=O")
output.write_text(result["flowchart_svg"], encoding="utf-8")
print(f"status={result.get('status')} category={result.get('potency_category')} score={result.get('potency_score')}")
print(f"svg={output} bytes={output.stat().st_size}")
