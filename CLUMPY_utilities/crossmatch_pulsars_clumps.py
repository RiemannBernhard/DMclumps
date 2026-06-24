#!/usr/bin/env python3
"""3-D cross-matcher: associate CLUMPY dark-matter clumps with pulsars.

This is the step the repo README lists as "not yet implemented". For every
pulsar it finds the dark-matter clump(s) whose centre lies within a tunable
multiple of the clump tidal radius in *three* dimensions (sky position **and**
heliocentric distance), evaluates the local Einasto density at the pulsar's
offset from the clump centre, and estimates the false-association rate by
re-running the identical matcher on Galactic-longitude-scrambled clump
catalogues (the look-elsewhere control the README asks for).

Inputs
------
Pulsars : the authoritative C-header catalogue shipped with the GW-template
  repo, ``NSmass_dist2GC_corr_data_extended_ext.h`` -- a ``[219][43]`` numeric
  array (Galactic l/b, heliocentric & galactocentric distance, NS mass, radius,
  age, spin frequency ...) plus an aligned ``NS_JName_catalog[219][3]`` name
  table. Galactic coordinates are already present, so no frame conversion is
  done here. Override the path with --pulsars or $DMCLUMPS_PSR_CATALOG.

Clumps  : one or more CLUMPY ``.drawn`` catalogues, read **directly** (not via
  the ROOT converter, whose committed output is stale and is documented to drop
  half the clumps and swap l/b). 21 whitespace columns; we use
  l(col2) b(col3) d(col4) Rdelta(col6) rhos(col7) rs(col8) alpha_Einasto(col10)
  Mtid(col16) Rtid(col17) Dgal(col20).  (0-based indices.)

Matching geometry
-----------------
Both populations are placed in heliocentric Cartesian coordinates (kpc):
  x = d cos b cos l,  y = d cos b sin l,  z = d sin b
and the Euclidean 3-D separation d3d is compared with alpha * R_tidal of each
clump. The nearest clump is reported for *every* pulsar regardless of the
threshold, so non-associations are visible too. Local density uses the Einasto
profile the converter assumes:
  rho(r) = rho_s * exp( -(2/alpha) * ((r/rs)^alpha - 1) )   [Msun/kpc^3]
converted to GeV/cm^3 via 37.96e-9. Because rs << R_tidal for these subhaloes,
rho falls steeply outside a few rs; d3d/rs is reported so the user can see how
deep inside the core (if at all) the pulsar actually sits.

Usage
-----
  python3 crossmatch_pulsars_clumps.py --drawn data/raw/100GeV/annihil_seed4448_*.drawn
  python3 crossmatch_pulsars_clumps.py --drawn 'data/raw/**/*.drawn' --alpha 1.0 \
          --scramble 100 --out data/crossmatch

Outputs (under --out, default data/crossmatch/):
  nearest_clump_per_pulsar.csv   one row per pulsar: nearest clump + association flag
  associations.csv               one row per (pulsar, clump) pair with d3d < alpha*Rtid
  crossmatch_summary.txt         counts, scramble-based false-association statistics
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

# Default pulsar catalogue: sibling GW-template repo under the same parent dir.
_DEFAULT_PSR = os.environ.get(
    "DMCLUMPS_PSR_CATALOG",
    os.path.join(os.path.dirname(REPO),
                 "ADM-NS_GWsignal_templ", "data",
                 "NSmass_dist2GC_corr_data_extended_ext.h"))

MSUNKPC3_TO_GEVCM3 = 37.96e-9  # rho: Msun/kpc^3 -> GeV/cm^3 (matches converter)

# Column indices in the .h numeric array (0-based), per the header legend.
PSR = dict(ra=0, dec=1, b=2, l=3, dist=4, dgc=5,
           mass=8, mass_ep=9, mass_em=10, radius_km=14, age_Gyr=17, f0=19)
# Column indices in the CLUMPY .drawn rows (0-based).
CL = dict(l=2, b=3, d=4, Rdelta=6, rhos=7, rs=8, alpha=10, Mtid=16, Rtid=17, Dgal=20)


# --------------------------------------------------------------------------- IO
def read_ns_catalog_h(path):
    """Parse the active NS catalogue from the C header.

    Returns a dict of numpy arrays plus parallel name/category lists. Only the
    *uncommented* ``NS_MNRAS_YiZhong_mass_and_unc[...][...]`` array and the
    aligned ``NS_JName_catalog`` are read; any /* ... */ block is skipped.
    """
    with open(path) as fh:
        text = fh.read()
    # Drop C block comments so the commented-out legacy copy is ignored.
    text_nc = re.sub(r"/\*.*?\*/", "", text, flags=re.S)

    def _rows(after_decl):
        m = re.search(after_decl, text_nc)
        if not m:
            raise KeyError(f"declaration not found: {after_decl!r}")
        body = text_nc[m.end():]
        body = body[:body.index("};")]
        return re.findall(r"\{([^{}]*)\}", body)

    num_rows = _rows(r"NS_MNRAS_YiZhong_mass_and_unc\s*\[\s*\d+\s*\]\s*\[\s*\d+\s*\]\s*=\s*\{")
    name_rows = _rows(r"NS_JName_catalog\s*\[\s*\d+\s*\]\s*\[\s*\d+\s*\]\s*=\s*\{")

    data = np.array([[float(x) for x in r.split(",")] for r in num_rows], dtype=float)
    names, cats, srcs = [], [], []
    for r in name_rows:
        parts = [p.strip().strip('"') for p in r.split('","')]
        parts = [p.strip().strip('"') for p in parts]
        srcs.append(parts[0] if len(parts) > 0 else "")
        cats.append(parts[1] if len(parts) > 1 else "")
        names.append(parts[2] if len(parts) > 2 else (parts[-1] if parts else ""))

    if len(names) != len(data):
        raise ValueError(f"name/data misalignment: {len(names)} names vs {len(data)} rows")

    return dict(
        n=len(data),
        name=names, category=cats, source=srcs,
        l=data[:, PSR["l"]], b=data[:, PSR["b"]], dist=data[:, PSR["dist"]],
        dgc=data[:, PSR["dgc"]], mass=data[:, PSR["mass"]],
        mass_ep=data[:, PSR["mass_ep"]], mass_em=data[:, PSR["mass_em"]],
        radius_km=data[:, PSR["radius_km"]], age_Gyr=data[:, PSR["age_Gyr"]],
        f0=data[:, PSR["f0"]],
    )


def read_drawn(paths):
    """Read one or more CLUMPY .drawn files into a single clump record dict.

    ``source`` tags which file each clump came from (basename); ``row`` is its
    line index within that file (1-based over data rows).
    """
    cols_l, cols_b, cols_d = [], [], []
    cols_Rtid, cols_rhos, cols_rs, cols_alpha = [], [], [], []
    cols_Mtid, cols_Dgal = [], []
    src_tag, row_idx = [], []
    nmax = max(CL.values())
    for p in paths:
        base = os.path.basename(p)
        n = 0
        with open(p) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.split()
                if len(f) <= nmax:
                    continue
                try:
                    l = float(f[CL["l"]]); b = float(f[CL["b"]]); d = float(f[CL["d"]])
                    Rtid = float(f[CL["Rtid"]]); rhos = float(f[CL["rhos"]])
                    rs = float(f[CL["rs"]]); alpha = float(f[CL["alpha"]])
                    Mtid = float(f[CL["Mtid"]]); Dgal = float(f[CL["Dgal"]])
                except ValueError:
                    continue
                n += 1
                cols_l.append(l); cols_b.append(b); cols_d.append(d)
                cols_Rtid.append(Rtid); cols_rhos.append(rhos)
                cols_rs.append(rs); cols_alpha.append(alpha)
                cols_Mtid.append(Mtid); cols_Dgal.append(Dgal)
                src_tag.append(base); row_idx.append(n)
    if not cols_l:
        raise ValueError("no clumps parsed from: " + ", ".join(paths))
    return dict(
        n=len(cols_l),
        l=np.array(cols_l), b=np.array(cols_b), d=np.array(cols_d),
        Rtid=np.array(cols_Rtid), rhos=np.array(cols_rhos),
        rs=np.array(cols_rs), alpha=np.array(cols_alpha),
        Mtid=np.array(cols_Mtid), Dgal=np.array(cols_Dgal),
        source=src_tag, row=row_idx,
    )


# ---------------------------------------------------------------------- geometry
def to_cartesian(l_deg, b_deg, d_kpc):
    l = np.radians(l_deg); b = np.radians(b_deg)
    cb = np.cos(b)
    return np.stack([d_kpc * cb * np.cos(l),
                     d_kpc * cb * np.sin(l),
                     d_kpc * np.sin(b)], axis=-1)


def einasto_rho_gevcm3(rhos, rs, alpha, r):
    """Einasto density [GeV/cm^3] at radius r (kpc) from the clump centre."""
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        rho_msun = rhos * np.exp(-(2.0 / alpha) * (np.power(r / rs, alpha) - 1.0))
    return rho_msun * MSUNKPC3_TO_GEVCM3


def einasto_density_radius(rhos, rs, alpha, rho_thresh_gevcm3):
    """Radius r [kpc] at which a clump's Einasto density drops to rho_thresh.

    Inverting rho(r)=rho_thresh gives, with C the Msun/kpc^3->GeV/cm^3 factor,
        (r/rs)^alpha = 1 - (alpha/2) * ln( rho_thresh / (C*rhos) )
    Inside this radius the clump's local DM density exceeds rho_thresh, so a
    pulsar is "in an overdensity" iff its 3-D offset from some clump centre is
    below that clump's density radius. This turns the (expensive) density test
    into a per-clump geometric radius -- exact, and cheap enough to scramble.
    Clumps whose *central* density never reaches rho_thresh get r=0.
    """
    with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
        bracket = 1.0 - (alpha / 2.0) * np.log(rho_thresh_gevcm3 / (MSUNKPC3_TO_GEVCM3 * rhos))
        r = np.where(bracket > 0.0, rs * np.power(bracket, 1.0 / alpha), 0.0)
    return np.nan_to_num(r, nan=0.0, posinf=0.0, neginf=0.0)


def nearest_clump(psr_xyz, clump_xyz):
    """For each pulsar return (index of nearest clump, 3-D distance kpc)."""
    idx = np.empty(len(psr_xyz), dtype=int)
    dmin = np.empty(len(psr_xyz), dtype=float)
    for i, p in enumerate(psr_xyz):
        d2 = np.sum((clump_xyz - p) ** 2, axis=1)
        j = int(np.argmin(d2))
        idx[i] = j
        dmin[i] = np.sqrt(d2[j])
    return idx, dmin


def count_associations(psr_xyz, clump_xyz, thresh):
    """Number of pulsars with >=1 clump within the per-clump radius ``thresh``."""
    t2 = thresh ** 2
    hits = 0
    for p in psr_xyz:
        d2 = np.sum((clump_xyz - p) ** 2, axis=1)
        if np.any(d2 < t2):
            hits += 1
    return hits


def pulsar_local_density(psr_xyz, clump, clump_xyz):
    """Exact local DM density at each pulsar, summed/maxed over all clumps.

    Returns (rho_sum, rho_max, idx_dominant) where rho_* are GeV/cm^3 and
    idx_dominant is the clump contributing the most density at that pulsar.
    Evaluates the Einasto profile of every clump at the pulsar offset, so this
    is the physically meaningful "ambient DM density at the NS" -- the quantity
    that drives capture/accretion in the seed paper.
    """
    rho_sum = np.empty(len(psr_xyz))
    rho_max = np.empty(len(psr_xyz))
    idom = np.empty(len(psr_xyz), dtype=int)
    rhos, rs, al = clump["rhos"], clump["rs"], clump["alpha"]
    for i, p in enumerate(psr_xyz):
        r = np.sqrt(np.sum((clump_xyz - p) ** 2, axis=1))
        rho = einasto_rho_gevcm3(rhos, rs, al, np.maximum(r, 1e-12))
        rho_sum[i] = rho.sum()
        j = int(np.argmax(rho))
        idom[i] = j
        rho_max[i] = rho[j]
    return rho_sum, rho_max, idom


def all_associations(psr_xyz, clump_xyz, thresh):
    """List of (pulsar_index, clump_index, d3d) with d3d < thresh[clump]."""
    out = []
    t2 = thresh ** 2
    for i, p in enumerate(psr_xyz):
        d2 = np.sum((clump_xyz - p) ** 2, axis=1)
        for j in np.nonzero(d2 < t2)[0]:
            out.append((i, int(j), float(np.sqrt(d2[j]))))
    return out


# -------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--drawn", nargs="+", required=True,
                    help="CLUMPY .drawn file(s) or glob(s)")
    ap.add_argument("--pulsars", default=_DEFAULT_PSR,
                    help="NS catalogue .h (default: sibling ADM-NS_GWsignal_templ)")
    ap.add_argument("--alpha", type=float, default=1.0,
                    help="tidal-membership factor: flag if d3d < alpha * R_tidal "
                         "(secondary, geometric; saturates for big clumps)")
    ap.add_argument("--rho-thresh", type=float, default=0.4, dest="rho_thresh",
                    help="primary association: pulsar counts as 'in an "
                         "overdensity' if local clump density exceeds this "
                         "[GeV/cm^3] (default 0.4 = canonical local halo)")
    ap.add_argument("--scramble", type=int, default=50,
                    help="number of longitude-scrambled catalogues for the "
                         "false-association estimate (0 to skip)")
    ap.add_argument("--seed", type=int, default=12345,
                    help="RNG seed for the scramble control (reproducible)")
    ap.add_argument("--out", default=os.path.join(REPO, "data", "crossmatch"),
                    help="output directory")
    args = ap.parse_args(argv)

    drawn_files = []
    for pat in args.drawn:
        hits = sorted(glob.glob(pat, recursive=True))
        drawn_files.extend(hits if hits else ([pat] if os.path.exists(pat) else []))
    if not drawn_files:
        ap.error("no .drawn files matched: " + " ".join(args.drawn))

    print(f"[in] pulsars: {args.pulsars}")
    psr = read_ns_catalog_h(args.pulsars)
    print(f"[in] {psr['n']} pulsars")
    print(f"[in] clumps : {len(drawn_files)} file(s)")
    cl = read_drawn(drawn_files)
    print(f"[in] {cl['n']} clumps")

    psr_xyz = to_cartesian(psr["l"], psr["b"], psr["dist"])
    clump_xyz = to_cartesian(cl["l"], cl["b"], cl["d"])
    tidal_r = args.alpha * cl["Rtid"]                       # geometric membership radius
    rho_r = einasto_density_radius(cl["rhos"], cl["rs"],    # density-threshold radius
                                   cl["alpha"], args.rho_thresh)

    print("[..] exact local DM density per pulsar (Einasto over all clumps)")
    rho_sum, rho_max, idom = pulsar_local_density(psr_xyz, cl, clump_xyz)

    print("[..] nearest clump + threshold associations")
    nidx, nd3d = nearest_clump(psr_xyz, clump_xyz)

    # Primary (density) and secondary (tidal) associations, both many-to-one.
    assoc_rho = all_associations(psr_xyz, clump_xyz, rho_r)
    assoc_tid = all_associations(psr_xyz, clump_xyz, tidal_r)
    n_rho = len({i for i, _, _ in assoc_rho})
    n_tid = len({i for i, _, _ in assoc_tid})

    os.makedirs(args.out, exist_ok=True)

    # ---- per-pulsar table: dominant (densest) clump + nearest clump ----------
    near_csv = os.path.join(args.out, "pulsar_clump_table.csv")
    with open(near_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["psr_index", "name", "category",
                    "psr_l_deg", "psr_b_deg", "psr_dist_kpc",
                    "psr_mass_Msun", "psr_radius_km", "psr_age_Gyr",
                    "rho_DM_sum_GeVcm3", "rho_DM_max_GeVcm3",
                    "in_overdensity",            # local density > rho_thresh
                    "dom_clump_source", "dom_clump_row", "dom_d3d_kpc",
                    "dom_d3d_over_rs", "dom_clump_Mtidal_Msun",
                    "near_clump_row", "near_d3d_kpc", "near_Rtidal_kpc",
                    "near_d3d_over_Rtidal", "in_tidal_radius"])
        for i in range(psr["n"]):
            jd = idom[i]
            dd = float(np.sqrt(np.sum((clump_xyz[jd] - psr_xyz[i]) ** 2)))
            jn = nidx[i]
            w.writerow([
                i, psr["name"][i], psr["category"][i],
                f"{psr['l'][i]:.5f}", f"{psr['b'][i]:.5f}", f"{psr['dist'][i]:.4f}",
                f"{psr['mass'][i]:.4f}", f"{psr['radius_km'][i]:.3g}",
                f"{psr['age_Gyr'][i]:.4g}",
                f"{rho_sum[i]:.6g}", f"{rho_max[i]:.6g}",
                int(rho_max[i] > args.rho_thresh),
                cl["source"][jd], cl["row"][jd], f"{dd:.6g}",
                f"{dd / cl['rs'][jd]:.4g}", f"{cl['Mtid'][jd]:.6g}",
                cl["row"][jn], f"{nd3d[i]:.6g}", f"{cl['Rtid'][jn]:.6g}",
                f"{nd3d[i] / cl['Rtid'][jn]:.4g}", int(nd3d[i] < tidal_r[jn])])

    # ---- density-association pairs (the scientifically meaningful matches) ----
    assoc_csv = os.path.join(args.out, "associations.csv")
    with open(assoc_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["psr_index", "name", "category", "clump_source", "clump_row",
                    "d3d_kpc", "rs_kpc", "d3d_over_rs",
                    "rho_DM_at_psr_GeVcm3", "clump_Mtidal_Msun"])
        rows = []
        for i, j, d3d in assoc_rho:
            rho = float(einasto_rho_gevcm3(cl["rhos"][j], cl["rs"][j],
                                           cl["alpha"][j], max(d3d, 1e-12)))
            rows.append((rho, i, j, d3d))
        for rho, i, j, d3d in sorted(rows, reverse=True):  # densest first
            w.writerow([i, psr["name"][i], psr["category"][i],
                        cl["source"][j], cl["row"][j],
                        f"{d3d:.6g}", f"{cl['rs'][j]:.6g}", f"{d3d / cl['rs'][j]:.4g}",
                        f"{rho:.6g}", f"{cl['Mtid'][j]:.6g}"])

    # ---- false-association control: scramble clump Galactic longitude --------
    # Run on the *density* criterion (the discriminating one); the tidal
    # criterion saturates and is reported only for context.
    scramble_txt = "scramble disabled (--scramble 0)"
    if args.scramble > 0:
        rng = np.random.default_rng(args.seed)
        counts = np.empty(args.scramble, dtype=int)
        for s in range(args.scramble):
            l_scram = rng.uniform(0.0, 360.0, size=cl["n"])
            xyz_s = to_cartesian(l_scram, cl["b"], cl["d"])
            counts[s] = count_associations(psr_xyz, xyz_s, rho_r)
        mu, sd = float(counts.mean()), float(counts.std())
        excess = n_rho - mu
        sig = excess / sd if sd > 0 else float("inf")
        scramble_txt = (
            f"scrambled trials       : {args.scramble} (Galactic-longitude randomised)\n"
            f"null overdensity assoc : {mu:.2f} +/- {sd:.2f} pulsars\n"
            f"observed overdensity   : {n_rho} pulsars\n"
            f"excess over null       : {excess:+.2f} ({sig:+.2f} sigma vs scramble)\n"
            f"false-association rate  : {(mu / max(psr['n'], 1)):.3f} per pulsar")

    summary = os.path.join(args.out, "crossmatch_summary.txt")
    n_over = int(np.sum(rho_max > args.rho_thresh))
    order = np.argsort(rho_max)[::-1]
    with open(summary, "w") as fh:
        fh.write("DMclumps 3-D pulsar<->clump cross-match summary\n")
        fh.write("=" * 52 + "\n")
        fh.write(f"pulsar catalogue       : {args.pulsars}\n")
        fh.write(f"clump files            : {len(drawn_files)}\n")
        for p in drawn_files:
            fh.write(f"    {p}\n")
        fh.write(f"pulsars                : {psr['n']}\n")
        fh.write(f"clumps                 : {cl['n']}\n")
        fh.write(f"rho threshold          : {args.rho_thresh} GeV/cm^3 (primary criterion)\n")
        fh.write(f"tidal alpha            : {args.alpha} (secondary, geometric)\n")
        fh.write("-" * 52 + "\n")
        fh.write(f"pulsars in overdensity : {n_over}  (rho_max > rho_thresh)\n")
        fh.write(f"  density assoc pairs  : {len(assoc_rho)}\n")
        fh.write(f"pulsars in tidal radius: {n_tid}  (d3d < alpha*R_tidal)\n")
        fh.write(f"median nearest d3d     : {np.median(nd3d):.4g} kpc\n")
        fh.write(f"min nearest d3d        : {nd3d.min():.4g} kpc "
                 f"({psr['name'][int(np.argmin(nd3d))]})\n")
        fh.write("-" * 52 + "\n")
        fh.write("top-5 pulsars by local DM density (rho_max):\n")
        for i in order[:5]:
            fh.write(f"  {psr['name'][i]:<22s} rho_max={rho_max[i]:.4g} GeV/cm^3  "
                     f"rho_sum={rho_sum[i]:.4g}  ({psr['category'][i]})\n")
        fh.write("-" * 52 + "\n")
        fh.write(scramble_txt + "\n")

    print(f"[ok] {n_over}/{psr['n']} pulsars in a DM overdensity (rho>{args.rho_thresh}); "
          f"{n_tid}/{psr['n']} inside a tidal radius")
    print(f"[ok] per-pulsar table -> {near_csv}")
    print(f"[ok] associations     -> {assoc_csv}")
    print(f"[ok] summary          -> {summary}")
    with open(summary) as fh:
        print("\n" + fh.read())
    return 0


if __name__ == "__main__":
    sys.exit(main())
