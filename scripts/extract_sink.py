#!/usr/bin/env python3
"""
Extract a sink CSV from an SBML model (XML or JSON) using biopathopt.

The sink CSV has two columns:
- Name  : Prefer MetaNetX chemical ID (metanetx.chemical); fallback to metabolite name
- InChI : From metabolite annotation ('inchi'), empty string if missing
"""

import argparse
import logging
import pandas as pd
from pathlib import Path
import biopathopt


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate sink.csv from an SBML model using biopathopt."
    )
    p.add_argument("--model", required=True, help="Path to SBML model (.xml or .json).")
    p.add_argument("--out", required=True, help="Path to output sink CSV (e.g., sink.csv).")
    return p


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    model_path = Path(args.model)
    out_path = Path(args.out)

    # Build the model
    bio_model = biopathopt.ModelBuilder(
        path_to_model=str(model_path),
        low_memory_mode=False,
        use_progressbar=False,
    )

    rows = []
    for m in bio_model.model.metabolites:
        # Prefer MetaNetX chemical ID for Name
        name = None
        mnx = m.annotation.get("metanetx.chemical")
        if mnx:
            if isinstance(mnx, list) and len(mnx) > 0:
                name = mnx[0]
                if len(mnx) > 1:
                    logging.warning(
                        f"There are {len(mnx)} MetaNetX IDs for {m.id}; using the first"
                    )
            elif isinstance(mnx, str):
                name = mnx

        if not name:
            logging.warning(f"No MetaNetX ID for {m.id}; using metabolite name")
            name = m.name

        inchi = m.annotation.get("inchi")
        if inchi:
            rows.append({"Name": name, "InChI": inchi})

    df = pd.DataFrame(rows).drop_duplicates()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
