from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the HRS–CHARLS PEF statistical analysis pipeline."
    )
    parser.add_argument("--hrs-data", required=True, type=Path, help="Path to the RAND HRS Stata file.")
    parser.add_argument(
        "--charls-data", required=True, type=Path, help="Path to the Harmonized CHARLS Stata file."
    )
    parser.add_argument("--skip-figures", action="store_true", help="Run analyses without figure scripts.")
    return parser.parse_args()


def run(command: list[str], env: dict[str, str]) -> None:
    print(f"Running: {' '.join(command)}")
    subprocess.run(command, cwd=ROOT, env=env, check=True)


def main() -> None:
    args = parse_args()
    hrs_data = args.hrs_data.expanduser().resolve()
    charls_data = args.charls_data.expanduser().resolve()
    for path, label in ((hrs_data, "HRS"), (charls_data, "CHARLS")):
        if not path.is_file():
            raise FileNotFoundError(f"{label} source file not found: {path}")

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError("Rscript was not found on PATH.")

    env = os.environ.copy()
    env.update(
        {
            "PEF_PROJECT_ROOT": str(ROOT),
            "HRS_DATA_FILE": str(hrs_data),
            "CHARLS_DATA_FILE": str(charls_data),
        }
    )

    python = sys.executable
    steps = [
        [python, "code/01_prepare_hrs.py"],
        [python, "code/02_prepare_charls.py"],
        [rscript, "code/03_primary_analysis.R"],
        [python, "code/04_build_sensitivity_datasets.py"],
        [rscript, "code/05_sensitivity_analysis.R"],
        [rscript, "code/06_absolute_risk_spline_analysis.R"],
    ]
    if not args.skip_figures:
        steps.extend(
            [
                [rscript, "code/07_create_main_forest.R"],
                [rscript, "code/08_create_sensitivity_forest.R"],
                [rscript, "code/09_create_clinical_value_figures.R"],
            ]
        )

    for command in steps:
        run(command, env)

    print("Pipeline completed. Generated files are under outputs/ and are excluded from Git.")


if __name__ == "__main__":
    main()

