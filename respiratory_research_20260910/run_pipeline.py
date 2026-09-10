"""Run the September 2026 revision in a private working directory."""
from pathlib import Path
import argparse
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', required=True, type=Path)
    parser.add_argument('--work-dir', required=True, type=Path)
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--skip-figures', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    scripts = Path(__file__).resolve().parent
    root = args.data_root.resolve()
    work = args.work_dir.resolve()
    if not args.dry_run:
        required = [root/'CHARLS date/Harmonized CHARLS/H_CHARLS_D_Data/H_CHARLS_D_Data.dta']
        if not all(p.is_file() for p in required):
            parser.error('Missing Harmonized CHARLS file; see README input layout.')
        if not list((root/'HRS data/01_rand_longitudinal').rglob('randhrs1992_2022v1.dta')):
            parser.error('Missing RAND HRS longitudinal file; see README input layout.')
        for year in [2006, 2010, 2014, 2018]:
            files = list((root/f'HRS data/02_rand_fat/{year}').rglob('*.dta'))
            if len(files) != 1:
                parser.error(f'Expected exactly one RAND FAT Stata file for {year}; found {len(files)}.')
        work.mkdir(parents=True, exist_ok=True)
    commands = [
        [sys.executable, str(scripts/'build_upgrade.py'), '--data-root', str(root), '--out', str(work/'data')],
        [args.rscript, str(scripts/'analyze_upgrade.R'), str(work/'data'), str(work/'results')],
        [args.rscript, str(scripts/'predict_upgrade.R'), str(work/'results')],
    ]
    if not args.skip_figures:
        commands.append([args.rscript, str(scripts/'make_figures.R'), str(work/'results'), str(work/'figures')])
    for command in commands:
        print(subprocess.list2cmdline(command), flush=True)
        if not args.dry_run:
            subprocess.run(command, check=True)


if __name__ == '__main__':
    main()
