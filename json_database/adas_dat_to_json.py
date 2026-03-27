#!/usr/bin/env python3
"""
Convert OpenADAS ADF11-style .dat to Hermes-3 / atomic1D JSON (same schema as plt96_c.json).

Grid layout is the same for 89 and 96 year files (plt96_c, plt89_b, …):
  - Header: IZ  N1  N2  …
  - After the first '---' line: N1 + N2 floats.
      • First N1 values: log10(ne / cm^-3), e.g. 7.7…15 (plt96) or 10…15 (10^10…10^15 cm^-3).
      • Next N2 values: log10(Te / eV), e.g. -0.7…4.2 (plt96) or 0…4.7 (1 eV … ~5×10^4 eV).
  - JSON log_density = file + 6 (log10(ne / m^-3)); log_temperature = file (log10 Te).
  - Coefficient table is n_T × n_n with n_T = N2, n_n = N1; file order row-major over (T, n).
  - log_coeff = dat_value - 6 (same as plt96_c.json).
  - Z1= blocks: n_T * n_n floats per charge state.
  - Output JSON rounds all floats to 5 decimal places (avoids 35.562470000000005).
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any

HEADER_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)")
# Fortran-style floats: -.69877, -32.93616, 10.00000
FLOAT_RE = re.compile(r"-?\d*\.\d+")
Z1_RE = re.compile(r"Z1=\s*(\d+)")


def parse_floats(line: str) -> list[float]:
    return [float(x) for x in FLOAT_RE.findall(line)]


def round_floats(obj: Any, ndigits: int = 5) -> Any:
    """Round every float in nested dict/list structures (JSON output)."""
    if isinstance(obj, float):
        return round(obj, ndigits)
    if isinstance(obj, list):
        return [round_floats(x, ndigits) for x in obj]
    if isinstance(obj, dict):
        return {k: round_floats(v, ndigits) for k, v in obj.items()}
    return obj


def read_dat(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    m = HEADER_RE.match(lines[0])
    if not m:
        raise ValueError(f"No IZ N1 N2 header: {path}")
    iz, n1, n2 = int(m.group(1)), int(m.group(2)), int(m.group(3))

    # Skip until after first line of dashes
    i = 1
    while i < len(lines) and "---" not in lines[i]:
        i += 1
    i += 1

    grid: list[float] = []
    while len(grid) < n1 + n2 and i < len(lines):
        if lines[i].strip().startswith("-" * 5):
            break
        grid.extend(parse_floats(lines[i]))
        i += 1
    if len(grid) != n1 + n2:
        raise ValueError(f"{path}: expected {n1 + n2} grid values, got {len(grid)}")

    # N1 = density axis (log10 cm^-3), N2 = temperature axis (log10 eV); table (n_T, n_n) = (N2, N1).
    n_grid = grid[:n1]
    T_grid = grid[n1 : n1 + n2]
    n_n, n_t = n1, n2

    log_temperature = list(T_grid)
    log_density = [x + 6.0 for x in n_grid]

    # Remaining Z1 blocks
    z1_starts: list[tuple[int, int]] = []
    for j, line in enumerate(lines):
        zm = Z1_RE.search(line)
        if zm:
            z1_starts.append((j, int(zm.group(1))))
    z1_starts.sort(key=lambda x: x[0])

    blocks: list[list[list[float]]] = []
    for bi, (start_line, z1) in enumerate(z1_starts):
        end_line = z1_starts[bi + 1][0] if bi + 1 < len(z1_starts) else len(lines)
        chunk = lines[start_line + 1 : end_line]
        nums: list[float] = []
        for ln in chunk:
            if ln.strip().startswith("-" * 5):
                break
            nums.extend(parse_floats(ln))
        need = n_t * n_n
        if len(nums) < need:
            raise ValueError(f"{path}: Z1={z1} expected {need} coeffs, got {len(nums)}")
        nums = nums[:need]
        table: list[list[float]] = []
        for r in range(n_t):
            row = [nums[r * n_n + c] - 6.0 for c in range(n_n)]
            table.append(row)
        blocks.append(table)

    if len(blocks) != iz:
        raise ValueError(f"{path}: expected {iz} Z1 blocks, got {len(blocks)}")

    stem = path.stem
    # e.g. plt96_c -> class plt, element c; acd89_b -> acd, b
    mname = re.match(r"^([a-z]+)\d+_(.+)$", stem, re.I)
    if not mname:
        raise ValueError(f"Cannot parse class/element from filename {stem}")
    cls, elem = mname.group(1).lower(), mname.group(2).lower()

    return {
        "charge": iz,
        "class": cls,
        "element": elem,
        "help": "JSON file corresponding to an OpenADAS data file\n"
        "Created by hermes-3 json_database/adas_dat_to_json.py (schema as TBody/atomic1D)\n"
        "Documentation at https://github.com/TBody/OpenADAS_to_JSON",
        "log_coeff": blocks,
        "log_density": log_density,
        "log_temperature": log_temperature,
        "name": str(path.resolve()),
        "number_of_charge_states": iz,
        "numpy_ndarrays": ["log_density", "log_temperature", "log_coeff"],
    }


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parent
    inputs = argv[1:] or [
        "acd89_b.dat",
        "ccd89_b.dat",
        "plt89_b.dat",
        "prb89_b.dat",
        "scd89_b.dat",
    ]
    for name in inputs:
        p = Path(name)
        if not p.is_file():
            p = root / name
        data = round_floats(read_dat(p))
        out = p.with_suffix(".json")
        out.write_text(json.dumps(data, indent=4) + "\n", encoding="utf-8")
        print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
