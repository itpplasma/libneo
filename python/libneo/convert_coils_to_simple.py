#!/usr/bin/env python3
"""
Convert STELLOPT coils format to simple biotsavart format.

Usage:
    python -m libneo.convert_coils_to_simple input.coils output.simple
"""

import sys
from pathlib import Path
from libneo.coils import CoilsFile


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    input_file = Path(sys.argv[1])
    output_file = Path(sys.argv[2])

    if not input_file.exists():
        print(f"Error: Input file {input_file} does not exist")
        sys.exit(1)

    print(f"Reading STELLOPT coils from {input_file}")
    coils = CoilsFile.from_file(str(input_file))

    n_filaments = len(coils.filaments)
    n_points = sum(len(fil.coords) for fil in coils.filaments)

    print(f"Found {n_filaments} filament(s) with {n_points} total points")
    if coils.periods is not None:
        print(f"Periods: {coils.periods}")

    print(f"Writing simple format to {output_file}")
    coils.write_simple(str(output_file))

    print("Conversion complete")


if __name__ == "__main__":
    main()
