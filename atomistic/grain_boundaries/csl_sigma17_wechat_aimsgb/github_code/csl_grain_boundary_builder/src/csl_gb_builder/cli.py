"""Command-line interface for configuration-driven CSL GB generation."""

from __future__ import annotations

import argparse
import json
from typing import Sequence

from .builder import build_grain_boundary
from .config import ConfigurationError, load_config


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate periodic CSL grain boundaries using aimsgb.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("build", help="Validate a JSON config and generate a structure.")
    build.add_argument("config", help="Path to a JSON configuration file.")
    build.add_argument("--dry-run", action="store_true", help="Validate and print the resolved configuration only.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        config = load_config(args.config)
        if args.dry_run:
            print(json.dumps({key: str(value) for key, value in config.__dict__.items()}, indent=2))
            return 0
        metadata = build_grain_boundary(config)
    except (ConfigurationError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 2
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
