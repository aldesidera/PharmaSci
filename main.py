"""Main entrypoint for the ChemoSuite web app and isolated modules."""

from __future__ import annotations

import argparse
import json
from typing import Optional, Sequence

from chemo_suite.apps.nitro_ra.cpca import calculate_cpca
from chemo_suite.main import run_molsim_web


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="ChemoSuite modular scientific application")
    parser.add_argument("--app", choices=["nitro_ra"], help="Application module to run in isolated mode")
    parser.add_argument("--module", choices=["cpca"], help="Scientific module to run")
    parser.add_argument("--smiles", help="Candidate nitrosamine SMILES")
    parser.add_argument("--mdd-mg", type=float, help="Maximum daily dose in mg/day for ppm conversion")
    return parser


def run_isolated_module(args: argparse.Namespace) -> int:
    if not args.app and not args.module and not args.smiles and args.mdd_mg is None:
        run_molsim_web()
        return 0

    if args.app != "nitro_ra" or args.module != "cpca":
        raise SystemExit("O modo isolado disponível nesta etapa é --app nitro_ra --module cpca.")
    if not args.smiles:
        raise SystemExit("--smiles é obrigatório no modo isolado cPCA.")

    result = calculate_cpca(args.smiles, mdd_mg=args.mdd_mg)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("status") not in {"invalid_input", "invalid_smiles"} else 2


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return run_isolated_module(args)


if __name__ == "__main__":
    raise SystemExit(main())
