"""
Handle STELLOPT/MAKEGRID filament coil files (e.g., coils.c09r00).

Style aligns with python/libneo/mgrid.py: provide a class with a
`from_file` constructor, a `write` method, metadata adapters, and a plot method.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Optional
import numpy as np


@dataclass(frozen=True)
class Filament:
    coords: np.ndarray   # shape (N, 3) -> columns x, y, z
    current: np.ndarray  # shape (N,)


@dataclass(frozen=True)
class Coils:
    periods: Optional[int]
    mirror: Optional[str]
    filaments: List[Filament]


def _parse_coils(path: str) -> Coils:
    """
    Parse a STELLOPT/MAKEGRID coils file (e.g., coils.c09r00) into Coils.

    - Recognizes: 'periods', 'begin filament', optional 'mirror'.
    - Within a filament block: lines with 4 numbers: x y z current.
    - Stores geometry and per-point current for each filament.
    """
    periods: Optional[int] = None
    mirror: Optional[str] = None
    filaments: List[Filament] = []
    in_filament = False
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    Is: List[float] = []
    with open(path, "r", encoding="utf-8") as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith("!"):
                continue
            # Remove inline comments
            for cmt in ("#", "!"):
                cpos = line.find(cmt)
                if cpos >= 0:
                    line = line[:cpos].strip()
            if not line:
                continue

            low = line.lower()
            if low.startswith("periods"):
                # periods N -> ignore for now
                continue
            if low.startswith("mirror"):
                # mirror NIL/ON/OFF -> record token for completeness
                toks = line.split()
                mirror = toks[1] if len(toks) > 1 else None
                continue
            if low.startswith("begin filament"):
                if in_filament:
                    # Close previous filament
                    if len(xs) == 0:
                        raise ValueError(f"Empty filament at {path}:{lineno}")
                    coords = np.column_stack([np.array(xs), np.array(ys), np.array(zs)])
                    cur = np.array(Is)
                    filaments.append(Filament(coords=coords, current=cur))
                    xs.clear(); ys.clear(); zs.clear(); Is.clear()
                in_filament = True
                continue

            # Inside filament: expect 4 numbers x y z I
            if in_filament:
                parts = line.split()
                if len(parts) < 4:
                    # End of filament
                    if len(xs) == 0:
                        raise ValueError(f"Empty filament at {path}:{lineno}")
                    coords = np.column_stack([np.array(xs), np.array(ys), np.array(zs)])
                    cur = np.array(Is)
                    filaments.append(Filament(coords=coords, current=cur))
                    xs.clear(); ys.clear(); zs.clear(); Is.clear()
                    in_filament = False
                    continue
                try:
                    x = float(parts[0]); y = float(parts[1]); z = float(parts[2]); cur = float(parts[3])
                except ValueError:
                    # End of filament on non-numeric line
                    if len(xs) == 0:
                        raise ValueError(f"Empty filament at {path}:{lineno}")
                    coords = np.column_stack([np.array(xs), np.array(ys), np.array(zs)])
                    cur_arr = np.array(Is)
                    filaments.append(Filament(coords=coords, current=cur_arr))
                    xs.clear(); ys.clear(); zs.clear(); Is.clear()
                    in_filament = False
                    continue
                xs.append(x); ys.append(y); zs.append(z); Is.append(cur)
                continue

            # Lines outside known blocks are ignored in this minimal parser
            continue

    # Close last filament if file ended while in one
    if in_filament:
        if len(xs) == 0:
            raise ValueError(f"Empty filament at end of file: {path}")
        coords = np.column_stack([np.array(xs), np.array(ys), np.array(zs)])
        cur = np.array(Is)
        filaments.append(Filament(coords=coords, current=cur))

    # Attempt to parse periods if present
    # Quick pass at top of file
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            low = line.lower()
            if low.startswith("periods"):
                toks = line.split()
                if len(toks) > 1:
                    try:
                        periods = int(toks[1])
                    except ValueError:
                        periods = None
                break

    return Coils(periods=periods, mirror=mirror, filaments=filaments)

class CoilsFile:
    def __init__(self, *, periods: Optional[int], mirror: Optional[str], filaments: List[Filament]):
        self.periods = periods
        self.mirror = mirror
        self.filaments = filaments

    def to_mgrid_metadata(self) -> Tuple[List[str], List[float]]:
        """
        Convert to Mgrid metadata arrays (coil_group, raw_coil_cur).
        Each filament maps to one group label 'filament' with its averaged current.
        """
        groups: List[str] = []
        currents: List[float] = []
        for fil in self.filaments:
            if fil.current.size == 0:
                continue
            nz = fil.current[np.nonzero(fil.current)]
            cur = float(nz.mean()) if nz.size > 0 else float(fil.current.mean())
            groups.append("filament")
            currents.append(cur)
        return groups, currents

    @classmethod
    def from_file(cls, filename: str) -> "CoilsFile":
        coils = _parse_coils(filename)
        return cls(periods=coils.periods, mirror=coils.mirror, filaments=coils.filaments)

    def write(self, filename: str) -> None:
        with open(filename, "w", encoding="utf-8") as f:
            if self.periods is not None:
                f.write(f"periods {self.periods}\n")
            if self.mirror is not None:
                f.write(f"mirror {self.mirror}\n")
            for fil in self.filaments:
                f.write("begin filament\n")
                for (x, y, z), I in zip(fil.coords, fil.current):
                    f.write(f"{x:.14E}   {y:.14E}   {z:.14E}   {I:.14E}\n")

    def plot(self, show_legend: bool = True):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        all_points = []
        for k, fil in enumerate(self.filaments, start=1):
            pts = fil.coords
            # Label with averaged current
            nz = fil.current[np.nonzero(fil.current)]
            cur = float(nz.mean()) if nz.size > 0 else float(fil.current.mean())
            label = f"filament_{k:02d} (I={cur:.3e} A)"
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], label=label)
            all_points.append(pts)
        if all_points:
            stacked = np.vstack(all_points)
            span = stacked.max(axis=0) - stacked.min(axis=0)
            span[span == 0.0] = 1.0
            ax.set_box_aspect(span)
        ax.set_xlabel("X [m]")
        ax.set_ylabel("Y [m]")
        ax.set_zlabel("Z [m]")
        if show_legend:
            ax.legend(fontsize="small")
        plt.tight_layout()
        plt.show()


__all__ = [
    "Filament",
    "Coils",
    "CoilsFile",
]
