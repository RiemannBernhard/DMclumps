#!/usr/bin/env python3
"""Stage 2, piece 2: project pulsar mass vs position inside the DM clump.

Implements the seed-paper test (MNRAS 527, 4483): treat each associated DM
clump as a miniature Galaxy and ask whether the accreted-DM / pulsar-mass
behaviour vs position inside the clump reproduces the galactocentric
M_NS(distance) trend of Phys. Dark Univ. 30, 100484 (2020).

Mechanism (the same law used for the whole Galaxy, applied inside each clump):
  rho_DM(r)  : the clump's own Einasto profile (rho_s, r_s, alpha from .drawn),
               evaluated at the pulsar's 3-D offset r from the clump centre.
  M_acc(r)   : accreted DM mass over the pulsar lifetime, accreted_mass()
               ported from convert_clumpy_draw_to_root.C (seed paper Eqs.9-11):
               M_acc = 1.3e43 * (rho/0.3) * t * (0.45*sigma_ratio) * GeV2Msun.
  Because rho_DM peaks in the core and falls toward the edge, M_acc is largest
  near the centre and negligible at the periphery -- the clump analogue of the
  galactocentric trend.

Input  : the piece-1 containment table (stage2_clump_radial.py),
         data/stage2/pulsar_in_clump_Rdelta.csv -- it already carries
         rho_DM_at_r, r, R_delta, the observed pulsar mass, and the age.
Outputs (under --out-dir):
  mass_projection_pairs.csv  : per (pulsar,clump) pair, + M_acc and radial coords
  mass_projection_binned.csv : aggregated vs distance-to-edge (the projection)

Usage:
  python3 pipeline/stage2_mass_projection.py \
      --table data/stage2/pulsar_in_clump_Rdelta.csv \
      --sigma-ratio 1 --time-Gyr 10 --primary-only --bins 10
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mtot_vs_mdm_digitized as mtot  # M_total^max(M_DM), Figs 10/11 of PRD 99 063015

# 1 GeV/c^2 in Msun (1 Msun = 1.111172751e57 GeV/c^2), as in the C++ source.
GEV_TO_MSUN = (1.0 / 1.111172751) * 1e-57


def accreted_mass(rho_gevcm3, sigma_ratio, time_Gyr):
    """Accumulated DM mass [Msun] -- port of accreted_mass() (seed paper)."""
    f = 0.45 * sigma_ratio
    return 1.3e43 * (rho_gevcm3 / 0.3) * time_Gyr * f * GEV_TO_MSUN


def _f(row, key, default=0.0):
    try:
        return float(row.get(key, default) or default)
    except (TypeError, ValueError):
        return default


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "data", "stage2", "pulsar_in_clump_Rdelta.csv"),
        help="piece-1 containment table")
    ap.add_argument("--sigma-ratio", type=float, default=1.0,
                    help="sigma_DM / sigma_crit")
    ap.add_argument("--time-Gyr", type=float, default=10.0,
                    help="accretion time [Gyr]; per-pulsar age used when > 0")
    ap.add_argument("--mchi", type=float, default=100.0,
                    help="DM particle mass [GeV] selecting the M_total(M_DM) "
                         "curve (snapped to nearest digitized: 1/5/10/50/100/200/500)")
    ap.add_argument("--y", type=float, default=0.1,
                    help="DM self-interaction strength: 0.1 (Fig.10 weak) or "
                         "1000 (Fig.11 strong)")
    ap.add_argument("--primary-only", action="store_true",
                    help="keep only the densest container per pulsar (is_primary)")
    ap.add_argument("--bins", type=int, default=10,
                    help="number of radial bins for the projection")
    ap.add_argument("--out-dir", default=None,
                    help="output dir (default: alongside --table)")
    args = ap.parse_args(argv)

    with open(args.table) as fh:
        rows = list(csv.DictReader(fh))
    if args.primary_only:
        rows = [r for r in rows if r.get("is_primary") == "1"]
    if not rows:
        ap.error("no rows (after primary filter) in " + args.table)

    out_dir = args.out_dir or os.path.dirname(args.table) or "."
    os.makedirs(out_dir, exist_ok=True)
    pairs_csv = os.path.join(out_dir, "mass_projection_pairs.csv")
    binned_csv = os.path.join(out_dir, "mass_projection_binned.csv")

    # ---- per-pair: accreted DM mass at the pulsar's location -----------------
    pair_fields = ["psr_index", "name", "category", "psr_mass_Msun",
                   "clump_source", "clump_row", "r_kpc", "Rdelta_kpc",
                   "r_over_Rdelta", "dist_to_edge_kpc", "frac_to_edge",
                   "rho_DM_at_r_GeVcm3", "time_Gyr", "M_acc_Msun",
                   "f_acc_over_Mpsr", "M_total_Msun", "dM_total_Msun"]
    data = []
    with open(pairs_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=pair_fields)
        w.writeheader()
        for r in rows:
            rho = _f(r, "rho_DM_at_r_GeVcm3")
            rk = _f(r, "r_kpc")
            Rd = _f(r, "Rdelta_kpc")
            xr = _f(r, "r_over_Rdelta")
            age = _f(r, "psr_age_Gyr")
            t = age if age > 0 else args.time_Gyr
            Mpsr = _f(r, "psr_mass_Msun")
            Macc = accreted_mass(rho, args.sigma_ratio, t)
            Mtot = mtot.mtot_of_mdm(Macc, y=args.y, m_chi_GeV=args.mchi)
            dMtot = Mtot - mtot.BASELINE_MT
            dist_edge = max(Rd - rk, 0.0)
            frac_edge = max(1.0 - xr, 0.0)
            rec = {
                "psr_index": r.get("psr_index", ""), "name": r.get("name", ""),
                "category": r.get("category", ""), "psr_mass_Msun": f"{Mpsr:.4g}",
                "clump_source": r.get("clump_source", ""),
                "clump_row": r.get("clump_row", ""),
                "r_kpc": f"{rk:.6g}", "Rdelta_kpc": f"{Rd:.6g}",
                "r_over_Rdelta": f"{xr:.5g}", "dist_to_edge_kpc": f"{dist_edge:.6g}",
                "frac_to_edge": f"{frac_edge:.5g}",
                "rho_DM_at_r_GeVcm3": f"{rho:.6g}", "time_Gyr": f"{t:.4g}",
                "M_acc_Msun": f"{Macc:.6g}",
                "f_acc_over_Mpsr": f"{(Macc / Mpsr if Mpsr > 0 else 0):.6g}",
                "M_total_Msun": f"{Mtot:.4g}", "dM_total_Msun": f"{dMtot:.4g}"}
            w.writerow(rec)
            data.append((xr, frac_edge, Mpsr, Macc, rho, Mtot))

    # ---- binned projection vs fractional distance-to-edge --------------------
    # x-axis = frac_to_edge = 1 - r/R_delta  (0 at clump edge, 1 at centre);
    # report mean observed pulsar mass and mean accreted DM mass per bin.
    nb = max(args.bins, 1)
    bins = [[] for _ in range(nb)]
    for xr, fe, Mpsr, Macc, rho, Mtot in data:
        k = min(int(fe * nb), nb - 1)
        bins[k].append((Mpsr, Macc, rho, xr, Mtot))
    with open(binned_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["bin", "frac_to_edge_lo", "frac_to_edge_hi", "n",
                    "mean_r_over_Rdelta", "mean_psr_mass_Msun",
                    "median_psr_mass_Msun", "mean_M_acc_Msun",
                    "mean_M_total_Msun", "mean_rho_GeVcm3"])
        for k in range(nb):
            lo, hi = k / nb, (k + 1) / nb
            b = bins[k]
            if not b:
                w.writerow([k, f"{lo:.3f}", f"{hi:.3f}", 0, "", "", "", "", "", ""])
                continue
            masses = sorted(m for m, _, _, _, _ in b)
            med = masses[len(masses) // 2]
            n = len(b)
            w.writerow([
                k, f"{lo:.3f}", f"{hi:.3f}", n,
                f"{sum(x for _, _, _, x, _ in b) / n:.4g}",
                f"{sum(m for m, _, _, _, _ in b) / n:.4g}", f"{med:.4g}",
                f"{sum(a for _, a, _, _, _ in b) / n:.6g}",
                f"{sum(mt for _, _, _, _, mt in b) / n:.4g}",
                f"{sum(r for _, _, r, _, _ in b) / n:.6g}"])

    snapped = mtot._nearest_mchi(args.y, args.mchi)
    print(f"[in] {len(rows)} pulsar-clump pairs"
          f"{' (primary only)' if args.primary_only else ''}; "
          f"sigma_ratio={args.sigma_ratio}, t={args.time_Gyr} Gyr (age where >0); "
          f"M_total curve: y={args.y}, m_chi={snapped} GeV")
    print(f"[ok] M_acc range:   {min(d[3] for d in data):.3g} .. "
          f"{max(d[3] for d in data):.3g} Msun")
    print(f"[ok] M_total range: {min(d[5] for d in data):.4g} .. "
          f"{max(d[5] for d in data):.4g} Msun (baseline {mtot.BASELINE_MT})")
    print(f"[ok] per-pair    -> {pairs_csv}")
    print(f"[ok] projection  -> {binned_csv}")
    # show the projection inline (centre -> edge)
    print("\nprojection (centre -> edge):")
    with open(binned_csv) as fh:
        for line in fh:
            print("  " + line.rstrip())
    return 0


if __name__ == "__main__":
    sys.exit(main())
