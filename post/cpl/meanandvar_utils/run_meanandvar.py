#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys


MEANANDVAR_DIR = Path(__file__).resolve().parent
DATA_FILE = Path(
    "/home/ws/th5982/research/projects/HST-S1S2/fortran-code/"
    "box_ar-effect/steady_state/data/Re-1000/B-4-3/A-0/"
    "statistics/Runtimedata"
)


def main() -> int:
    executable = MEANANDVAR_DIR / "meanandvar"

    if not executable.is_file():
        print(f"error: executable not found: {executable}", file=sys.stderr)
        return 1
    if not DATA_FILE.is_file():
        print(f"error: input file not found: {DATA_FILE}", file=sys.stderr)
        return 1

    try:
        with DATA_FILE.open("rb") as input_file:
            completed = subprocess.run(
                ["cpl", "run", str(executable)],
                stdin=input_file,
                check=False,
            )
    except FileNotFoundError:
        print("error: cpl command was not found in PATH", file=sys.stderr)
        return 1

    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())