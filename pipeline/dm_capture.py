#!/usr/bin/env python3
"""Stage 1->2 glue, part A: dark-matter capture onto neutron stars.

Turns the cross-matcher's per-pulsar local DM density into an accumulated DM
mass and DM mass fraction f_DM = M_DM / M_NS, which is the quantity the
two-fluid TOV solver (TOV-dark-sector-stars) consumes (see tov_bridge.py).

The capture physics is a faithful Python port of the formulae already used in
the published analysis (``CLUMPY_utilities/convert_clumpy_draw_to_root.C``,
functions ``accreted_mass`` / ``accreted_rate``), i.e. Deliyergiyev et al.
(2024), MNRAS 527, 4483, Eqs. 9-11.

Two equivalent routes to the accumulated mass are provided:

* ``accreted_mass``        -- the time-integrated closed form (independent of the
                              DM particle mass; fiducial NS geometry baked into
                              the 1.3e43 constant).
* ``accreted_mass_from_rate`` -- capture *rate* (Eq. 9) x particle mass x time,
                              scaled by the actual NS mass & radius. Preferred
                              when per-NS M, R are known. The rate goes as
                              1/m_chi and the per-capture mass as m_chi, so the
                              accumulated mass is m_chi-independent -- consistent
                              with ``accreted_mass``.

Sign of the result for realistic inputs: at the canonical local density
rho=0.4 GeV/cm^3, sigma_ratio=1, t=10 Gyr the accreted mass is ~7e-14 Msun
(f_DM ~ 5e-14) -- far below any measurable mass/Lambda change. An observable
f_DM (>~1e-3) needs either sigma_ratio >~ 1e10 or the thermalised central spike
density (rho(r_min) ~ 7.7e9 GeV/cm^3 in the seed paper), not field/clump
densities. ``required_sigma_ratio`` inverts the linear relation to report, per
pulsar, the cross-section that would be needed to reach a target f_DM.

Usage (apply to a cross-matcher table):
  python3 pipeline/dm_capture.py \
      --table data/crossmatch/pulsar_clump_table.csv \
      --rho-col rho_DM_max_GeVcm3 --sigma-ratio 1.0 --time-Gyr 10 \
      --target-fdm 1e-3 --out data/crossmatch/capture.csv
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

# 1 GeV/c^2 in Msun  (1 Msun = 1.111172751e57 GeV/c^2), as in the C++ source.
GEV_TO_MSUN = (1.0 / 1.111172751) * 1e-57
RSUN_KM = 695700.0           # solar radius [km], for R_NS [km] -> [Rsun]
SEC_PER_GYR = 3.1557e16      # Julian Gyr in seconds


def accreted_mass(rho_gevcm3, sigma_ratio, time_Gyr):
    """Accumulated DM mass [Msun] -- closed form (port of ``accreted_mass``).

    rho_gevcm3  : local DM density [GeV/cm^3]
    sigma_ratio : sigma_DM / sigma_crit  (sigma_crit ~ 6e-46 cm^2)
    time_Gyr    : accumulation time [Gyr] (formation time, 0.49..13.7)
    """
    f = 0.45 * sigma_ratio
    return 1.3e43 * (rho_gevcm3 / 0.3) * time_Gyr * f * GEV_TO_MSUN


def accreted_rate(rho_gevcm3, sigma_ratio, m_chi_GeV, M_NS_Msun, R_NS_Rsun,
                  velocity_kms=10.0):
    """DM capture rate [particles/s] -- port of ``accreted_rate`` (Eq. 9)."""
    f = 0.45 * sigma_ratio
    return (1.1e27 * (rho_gevcm3 / 0.3) * (220.0 / velocity_kms)
            * (1e3 / m_chi_GeV) * M_NS_Msun * R_NS_Rsun * f)


def accreted_mass_from_rate(rho_gevcm3, sigma_ratio, m_chi_GeV,
                            M_NS_Msun, R_NS_Rsun, time_Gyr, velocity_kms=10.0):
    """Accumulated DM mass [Msun] via rate x m_chi x time, scaled by NS M,R.

    rate[#/s] * m_chi[GeV] * GEV_TO_MSUN[Msun/GeV] * time[s].
    The explicit m_chi cancels the 1/m_chi in the rate, so the result does not
    depend on m_chi (kept as an argument only for transparency / Pauli-blocking
    extensions).
    """
    rate = accreted_rate(rho_gevcm3, sigma_ratio, m_chi_GeV,
                          M_NS_Msun, R_NS_Rsun, velocity_kms)
    return rate * m_chi_GeV * GEV_TO_MSUN * (time_Gyr * SEC_PER_GYR)


def dm_fraction(M_DM_Msun, M_NS_Msun):
    """DM mass fraction f_DM = M_DM / M_NS (the TOV input knob)."""
    return M_DM_Msun / M_NS_Msun if M_NS_Msun > 0 else float("nan")


def required_sigma_ratio(target_fdm, rho_gevcm3, M_NS_Msun, time_Gyr):
    """sigma_DM/sigma_crit needed to reach ``target_fdm`` via ``accreted_mass``.

    accreted_mass is linear in sigma_ratio, so invert directly. Returns inf if
    the density is zero.
    """
    per_unit = accreted_mass(rho_gevcm3, 1.0, time_Gyr)  # M_DM at sigma_ratio=1
    if per_unit <= 0:
        return float("inf")
    return target_fdm * M_NS_Msun / per_unit


# --------------------------------------------------------------------------- CLI
def _f(row, key, default=0.0):
    try:
        return float(row.get(key, default) or default)
    except (TypeError, ValueError):
        return default


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True,
                    help="cross-matcher pulsar_clump_table.csv")
    ap.add_argument("--rho-col", default="rho_DM_max_GeVcm3",
                    help="density column to use (rho_DM_max_GeVcm3 or "
                         "rho_DM_sum_GeVcm3); or pass --rho to override per-row")
    ap.add_argument("--rho", type=float, default=None,
                    help="fixed DM density [GeV/cm^3] for every pulsar (e.g. the "
                         "thermalised spike value 7.7e9), overrides --rho-col")
    ap.add_argument("--sigma-ratio", type=float, default=1.0,
                    help="sigma_DM / sigma_crit")
    ap.add_argument("--mchi", type=float, default=100.0,
                    help="DM particle mass [GeV] (rate route only; M_DM is "
                         "m_chi-independent)")
    ap.add_argument("--time-Gyr", type=float, default=10.0,
                    help="accumulation time [Gyr]; per-pulsar psr_age_Gyr is "
                         "used when > 0, else this default")
    ap.add_argument("--default-radius-km", type=float, default=11.0,
                    help="NS radius [km] when the catalogue value is 0")
    ap.add_argument("--target-fdm", type=float, default=1e-3,
                    help="reference f_DM for the required-sigma sensitivity column")
    ap.add_argument("--out", default=None, help="output CSV (default: alongside --table)")
    args = ap.parse_args(argv)

    with open(args.table) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        ap.error("empty table: " + args.table)

    out = args.out or os.path.join(os.path.dirname(args.table) or ".", "capture.csv")
    fields = ["psr_index", "name", "category", "rho_DM_GeVcm3", "M_NS_Msun",
              "R_NS_km", "time_Gyr", "sigma_ratio",
              "M_DM_closed_Msun", "M_DM_rate_Msun", "f_DM_closed", "f_DM_rate",
              f"sigma_ratio_for_fdm_{args.target_fdm:g}"]
    n_obs = 0
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            rho = args.rho if args.rho is not None else _f(r, args.rho_col)
            M = _f(r, "psr_mass_Msun")
            R_km = _f(r, "psr_radius_km")
            if R_km <= 0:
                R_km = args.default_radius_km
            age = _f(r, "psr_age_Gyr")
            t = age if age > 0 else args.time_Gyr
            M_closed = accreted_mass(rho, args.sigma_ratio, t)
            M_rate = accreted_mass_from_rate(rho, args.sigma_ratio, args.mchi,
                                             M if M > 0 else 1.4,
                                             R_km / RSUN_KM, t)
            f_closed = dm_fraction(M_closed, M if M > 0 else 1.4)
            f_rate = dm_fraction(M_rate, M if M > 0 else 1.4)
            sig_need = required_sigma_ratio(args.target_fdm, rho,
                                            M if M > 0 else 1.4, t)
            if f_closed > args.target_fdm:
                n_obs += 1
            w.writerow({
                "psr_index": r.get("psr_index", ""), "name": r.get("name", ""),
                "category": r.get("category", ""),
                "rho_DM_GeVcm3": f"{rho:.6g}", "M_NS_Msun": f"{M:.4g}",
                "R_NS_km": f"{R_km:.4g}", "time_Gyr": f"{t:.4g}",
                "sigma_ratio": f"{args.sigma_ratio:.6g}",
                "M_DM_closed_Msun": f"{M_closed:.6g}",
                "M_DM_rate_Msun": f"{M_rate:.6g}",
                "f_DM_closed": f"{f_closed:.6g}", "f_DM_rate": f"{f_rate:.6g}",
                f"sigma_ratio_for_fdm_{args.target_fdm:g}": f"{sig_need:.6g}"})

    fmax = max((_f(r, args.rho_col) for r in rows), default=0.0)
    print(f"[ok] capture applied to {len(rows)} pulsars -> {out}")
    print(f"[ok] sigma_ratio={args.sigma_ratio}, time={args.time_Gyr} Gyr (per-PSR "
          f"age where available); rho from "
          f"{'--rho '+str(args.rho) if args.rho is not None else args.rho_col}")
    print(f"[ok] pulsars exceeding f_DM>{args.target_fdm:g}: {n_obs}/{len(rows)}  "
          f"(max table rho={fmax:.4g} GeV/cm^3)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
