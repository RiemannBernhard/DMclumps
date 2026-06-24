"""Stage 1->2 glue, part B: bridge to the two-fluid TOV solver.

Maps a target DM fraction f_DM (from dm_capture.py) onto a run of
TOV-dark-sector-stars and reads back the predicted observables
(M_total, R, Lambda) for comparison with the measured NS mass.

Three independent capabilities, so the parts that do not need the compiled
solver work without a ROOT build:

1. ``write_tov_config`` -- emit a solver JSON config (README schema).
2. ``run_tov``          -- invoke the compiled binary on a config (needs the
                           ROOT/VDT build; raises a clear error if absent).
3. ``parse_phys`` / ``star_at_fdm`` / ``star_at_max_mass`` -- read the physical
                           output table and pick the star at a target DM
                           fraction or the maximum-mass star. Pure text, always
                           runnable.

facdar vs f_DM
--------------
The solver's ``facdar`` is a DM central-pressure knob, *not* f_DM directly; the
realised fraction is ``massDM/massT`` in the output. To hit a target f_DM you
scan facdar once and interpolate (``facdar_for_fdm``); this module exposes the
achieved fraction so that calibration is explicit rather than assumed.

Physical output table columns (README "Output"), 1-based -> our names:
  1 precNM 2 enNM 3 rhocNM_rho0 4 radNM 5 massNM 6 precDM 7 enDM 8 radDM
  9 massDM 10 massT 11 YR 12 C_total 13 C_ord 14 C_dar 15 k2 16 Lambda
  [17.. rotation: Ibar I_gcm2 I_dm_frac q_kerr ellip fK Mrot_Kep Req_Kep]
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess

# 0-based column indices in rad-mass_phys_*.dat.
PHYS = dict(precNM=0, enNM=1, rhocNM_rho0=2, radNM=3, massNM=4,
            precDM=5, enDM=6, radDM=7, massDM=8, massT=9, YR=10,
            C_total=11, C_ord=12, C_dar=13, k2=14, Lambda=15)

DEFAULT_MN = 0.939565379  # neutron mass [GeV]


# ------------------------------------------------------------------- config I/O
def write_tov_config(path, eos_nm, eos_dm, output_rm, output_phys,
                     facord=1.0, facdar=1.0, idar=2, mdar=1, mn=DEFAULT_MN):
    """Write a solver JSON config (matches TOV-dark-sector-stars README schema).

    facdar=0.0 -> pure baryonic (single-fluid) star; >0 activates the DM fluid.
    idar selects the DM self-interaction column of eos_dm; mdar is the DM mass
    in neutron-mass units.
    """
    cfg = {
        "input_files": {
            "eos_dm": eos_dm, "eos_nm": eos_nm,
            "output_rm": output_rm, "output_phys": output_phys,
        },
        "parameters": {
            "facord": float(facord), "facdar": float(facdar),
            "idar": int(idar), "mdar": int(mdar), "mn": float(mn),
        },
    }
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        json.dump(cfg, fh, indent=2)
    return path


def find_tov_binary(tov_bin=None):
    """Locate the solver binary (arg > $TOV_BIN > sibling build dirs > PATH)."""
    cands = []
    if tov_bin:
        cands.append(tov_bin)
    if os.environ.get("TOV_BIN"):
        cands.append(os.environ["TOV_BIN"])
    here = os.path.dirname(os.path.abspath(__file__))
    repo_parent = os.path.dirname(os.path.dirname(here))  # ~/TestArea
    for d in ("build-rel", "build"):
        cands.append(os.path.join(repo_parent, "TOV-dark-sector-stars",
                                  d, "bin", "tov_2k_dbg"))
    cands.append(shutil.which("tov_2k_dbg") or "")
    for c in cands:
        if c and os.path.exists(c) and os.access(c, os.X_OK):
            return c
    return None


def run_tov(config_path, tov_bin=None, cwd=None, env_extra=None, timeout=None):
    """Run the solver on a config; return the output_phys path.

    Raises FileNotFoundError (with build guidance) if the binary is absent --
    this environment has no ROOT build, so callers should catch it and fall
    back to ``parse_phys`` on a precomputed table.
    """
    exe = find_tov_binary(tov_bin)
    if not exe:
        raise FileNotFoundError(
            "TOV solver binary not found. Build TOV-dark-sector-stars "
            "(needs CERN ROOT + VDT) per its README, or set $TOV_BIN. "
            "Until then use parse_phys() on an existing rad-mass_phys_*.dat.")
    with open(config_path) as fh:
        cfg = json.load(fh)
    out_phys = cfg["input_files"]["output_phys"]
    env = dict(os.environ)
    if env_extra:
        env.update({k: str(v) for k, v in env_extra.items()})
    subprocess.run([exe, config_path], check=True, cwd=cwd, env=env,
                   timeout=timeout)
    return os.path.join(cwd or ".", out_phys) if cwd else out_phys


# ----------------------------------------------------------------- output parse
def parse_phys(path):
    """Parse a rad-mass_phys_*.dat into a list of named-dict rows."""
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s[0] in "#":
                continue
            parts = s.split()
            if len(parts) <= PHYS["Lambda"]:
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            row = {name: vals[i] for name, i in PHYS.items()}
            row["_extra"] = vals[PHYS["Lambda"] + 1:]  # rotation cols if present
            rows.append(row)
    return rows


def star_at_max_mass(rows):
    """Return the maximum-(total-)mass star row."""
    return max(rows, key=lambda r: r["massT"]) if rows else None


def star_at_fdm(rows, target_fdm, tol=None):
    """Return the row whose realised DM fraction massDM/massT is closest to
    target_fdm (optionally require within ``tol``)."""
    if not rows:
        return None
    def fdm(r):
        return r["massDM"] / r["massT"] if r["massT"] > 0 else float("inf")
    best = min(rows, key=lambda r: abs(fdm(r) - target_fdm))
    if tol is not None and abs(fdm(best) - target_fdm) > tol:
        return None
    return best


def star_at_mass(rows, target_mass_Msun):
    """Return the row whose total mass is closest to target_mass_Msun."""
    return min(rows, key=lambda r: abs(r["massT"] - target_mass_Msun)) if rows else None


def facdar_for_fdm(facdar_fdm_pairs, target_fdm):
    """Linear-interpolate the facdar that yields target_fdm from calibration
    pairs [(facdar, achieved_fdm), ...] (sorted by fdm)."""
    pts = sorted(facdar_fdm_pairs, key=lambda p: p[1])
    for (f0, y0), (f1, y1) in zip(pts, pts[1:]):
        if y0 <= target_fdm <= y1 and y1 > y0:
            return f0 + (f1 - f0) * (target_fdm - y0) / (y1 - y0)
    return pts[-1][0] if target_fdm > pts[-1][1] else pts[0][0]


def observables(row):
    """Compact (M_total, R, Lambda, f_DM, k2) tuple from a phys row."""
    if row is None:
        return None
    f = row["massDM"] / row["massT"] if row["massT"] > 0 else float("nan")
    return dict(M_total=row["massT"], M_NM=row["massNM"], M_DM=row["massDM"],
                R_km=row["radNM"], Lambda=row["Lambda"], k2=row["k2"], f_DM=f)
