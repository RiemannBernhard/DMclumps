#!/usr/bin/env python3
"""Stage 2, piece 1: locate each pulsar inside its associated DM clump(s).

The seed paper (MNRAS 527, 4483) builds the relation "neutron-star mass vs
position inside the clump" by treating each associated DM clump as a miniature
Galaxy (its internal Einasto profile is self-similar to the galactic one) and
asking whether the M_NS(r_clump) trend reproduces the galactocentric
M_NS(distance) prediction of Phys. Dark Univ. 30, 100484 (2020).

This module produces the geometric foundation for that test: for every pulsar,
across every CLUMPY ``.drawn`` realization, it finds the clump(s) that
*contain* the pulsar and records where inside each clump the pulsar sits,
together with that clump's known properties (rho_s, r_s, R_tidal, R_delta,
M_tidal, M_delta, M_equidens, R_equidens) -- everything the Stage-2 mass model
needs. No mass model is applied here; this is the (pulsar, clump, position)
table the projection is built on.

Containment
-----------
A pulsar is "inside" a clump when its 3-D offset from the clump centre is below
the clump's outer radius. R_tidal is the default boundary (the physical extent
of the subhalo); ``--boundary`` can switch to R_delta or R_equidens. Both
populations are placed in heliocentric Cartesian coordinates (kpc):
  x = d cos b cos l, y = d cos b sin l, z = d sin b.
A pulsar may be contained in several (nested) clumps; every containing pair is
emitted (long format). The densest container (largest rho_DM at the pulsar's
position) is flagged ``is_primary=1`` as a convenience for one-clump-per-pulsar
projections.

Output (long format, one row per contained pulsar-clump pair):
  psr_index name category psr_mass_Msun psr_mass_ep psr_mass_em psr_age_Gyr
  psr_l_deg psr_b_deg psr_dist_kpc
  clump_source clump_row r_kpc r_over_Rtidal r_over_Rdelta r_over_rs
  rho_s_GeVcm3 rho_DM_at_r_GeVcm3
  Rtidal_kpc Rdelta_kpc Requidens_kpc Mtidal_Msun Mdelta_Msun Mequidens_Msun
  Dgal_kpc is_primary

Usage:
  python3 pipeline/stage2_clump_radial.py \
      --drawn 'data/raw/100GeV/*.drawn' 'data/raw/500GeV/*.drawn' \
      --boundary Rtidal --out data/stage2/pulsar_in_clump.csv
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "CLUMPY_utilities"))
import crossmatch_pulsars_clumps as xm  # canonical .drawn / NS-catalog readers

MSUNKPC3_TO_GEVCM3 = xm.MSUNKPC3_TO_GEVCM3

BOUNDARY_COL = {"Rtidal": "Rtid", "Rdelta": "Rdelta", "Requidens": "Requidens"}


def contained_pairs(psr_xyz, clump_xyz, boundary_radius):
    """Yield (pulsar_index, clump_index, r_kpc) for every containment.

    boundary_radius is the per-clump outer radius (kpc). Vectorised over clumps,
    looped over pulsars (219 pulsars x N clumps).
    """
    b2 = boundary_radius ** 2
    for i, p in enumerate(psr_xyz):
        d2 = np.sum((clump_xyz - p) ** 2, axis=1)
        inside = np.nonzero(d2 < b2)[0]
        for j in inside:
            yield i, int(j), float(np.sqrt(d2[j]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--drawn", nargs="+", required=True,
                    help="CLUMPY .drawn file(s) or glob(s)")
    ap.add_argument("--pulsars", default=xm._DEFAULT_PSR,
                    help="NS catalogue .h (default: sibling ADM-NS_GWsignal_templ)")
    ap.add_argument("--boundary", choices=list(BOUNDARY_COL), default="Rtidal",
                    help="clump outer radius defining containment (default Rtidal)")
    ap.add_argument("--out", default=os.path.join(REPO, "data", "stage2",
                                                  "pulsar_in_clump.csv"))
    args = ap.parse_args(argv)

    drawn_files = []
    for pat in args.drawn:
        hits = sorted(glob.glob(pat, recursive=True))
        drawn_files.extend(hits if hits else ([pat] if os.path.exists(pat) else []))
    if not drawn_files:
        ap.error("no .drawn files matched: " + " ".join(args.drawn))

    psr = xm.read_ns_catalog_h(args.pulsars)
    psr_xyz = xm.to_cartesian(psr["l"], psr["b"], psr["dist"])
    bcol = BOUNDARY_COL[args.boundary]
    print(f"[in] {psr['n']} pulsars; {len(drawn_files)} .drawn file(s); "
          f"boundary={args.boundary}")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fields = ["psr_index", "name", "category", "psr_mass_Msun", "psr_mass_ep",
              "psr_mass_em", "psr_age_Gyr", "psr_l_deg", "psr_b_deg",
              "psr_dist_kpc", "clump_source", "clump_row", "r_kpc",
              "r_over_Rtidal", "r_over_Rdelta", "r_over_rs",
              "rho_s_GeVcm3", "rho_DM_at_r_GeVcm3", "Rtidal_kpc", "Rdelta_kpc",
              "Requidens_kpc", "Mtidal_Msun", "Mdelta_Msun", "Mequidens_Msun",
              "Dgal_kpc", "is_primary"]

    n_pairs = 0
    pulsars_contained = set()
    # accumulate per-pulsar to mark the densest container as primary
    per_psr = {}  # i -> list of dict rows
    n_skipped_inf = 0
    for fpath in drawn_files:
        cl = xm.read_drawn([fpath])
        clump_xyz = xm.to_cartesian(cl["l"], cl["b"], cl["d"])
        # CLUMPY writes a non-finite (inf/nan) outer radius for a handful of
        # clumps; those would "contain" every pulsar. Zero out their boundary
        # so they can never match (b2 = 0 -> d2 < 0 is never true).
        boundary = np.array(cl[bcol], dtype=float)
        bad = ~(np.isfinite(boundary) & (boundary > 0.0))
        n_skipped_inf += int(bad.sum())
        boundary[bad] = 0.0
        nf = 0
        for i, j, r in contained_pairs(psr_xyz, clump_xyz, boundary):
            rho_r = float(xm.einasto_rho_gevcm3(cl["rhos"][j], cl["rs"][j],
                                                cl["alpha"][j], max(r, 1e-12)))
            rho_s = float(cl["rhos"][j] * MSUNKPC3_TO_GEVCM3)
            row = {
                "psr_index": i, "name": psr["name"][i],
                "category": psr["category"][i],
                "psr_mass_Msun": f"{psr['mass'][i]:.4f}",
                "psr_mass_ep": f"{psr['mass_ep'][i]:.4f}",
                "psr_mass_em": f"{psr['mass_em'][i]:.4f}",
                "psr_age_Gyr": f"{psr['age_Gyr'][i]:.4g}",
                "psr_l_deg": f"{psr['l'][i]:.5f}", "psr_b_deg": f"{psr['b'][i]:.5f}",
                "psr_dist_kpc": f"{psr['dist'][i]:.4f}",
                "clump_source": cl["source"][j], "clump_row": cl["row"][j],
                "r_kpc": f"{r:.6g}",
                "r_over_Rtidal": f"{r / cl['Rtid'][j]:.5g}",
                "r_over_Rdelta": (f"{r / cl['Rdelta'][j]:.5g}"
                                  if cl["Rdelta"][j] > 0 else ""),
                "r_over_rs": f"{r / cl['rs'][j]:.5g}",
                "rho_s_GeVcm3": f"{rho_s:.6g}",
                "rho_DM_at_r_GeVcm3": f"{rho_r:.6g}",
                "Rtidal_kpc": f"{cl['Rtid'][j]:.6g}",
                "Rdelta_kpc": f"{cl['Rdelta'][j]:.6g}",
                "Requidens_kpc": f"{cl['Requidens'][j]:.6g}",
                "Mtidal_Msun": f"{cl['Mtid'][j]:.6g}",
                "Mdelta_Msun": f"{cl['Mdelta'][j]:.6g}",
                "Mequidens_Msun": f"{cl['Mequidens'][j]:.6g}",
                "Dgal_kpc": f"{cl['Dgal'][j]:.5g}",
                "_rho_r": rho_r,  # for primary selection, dropped on write
            }
            per_psr.setdefault(i, []).append(row)
            pulsars_contained.add(i)
            n_pairs += 1
            nf += 1
        print(f"[..] {os.path.basename(fpath):60s} clumps={cl['n']:>7d} "
              f"containments={nf}")

    # mark densest container per pulsar as primary
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for i, rows in sorted(per_psr.items()):
            jmax = max(range(len(rows)), key=lambda k: rows[k]["_rho_r"])
            for k, row in enumerate(rows):
                row = dict(row)
                row["is_primary"] = 1 if k == jmax else 0
                row.pop("_rho_r", None)
                w.writerow(row)

    print(f"[ok] {n_pairs} pulsar-clump containments; "
          f"{len(pulsars_contained)}/{psr['n']} pulsars contained -> {args.out}")
    if n_skipped_inf:
        print(f"[ok] excluded {n_skipped_inf} clump(s) with non-finite/zero "
              f"{args.boundary} (CLUMPY sentinels)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
