#!/usr/bin/env python3
"""Digitized M_total^max(M_DM) curves -- NS branch of Figs. 10 & 11 of
Deliyergiyev et al., PRD 99, 063015 (2019) [arXiv:1903.01183].

These are the structural-response curves that 100484 (Fig. 2) feeds the accreted
DM mass through to turn M_DM(galactocentric distance) into M_total(distance).
Stage-2 piece 2 uses the same lookup with the clump's M_acc(r).

!! APPROXIMATE !!  Read by eye off the published NS-branch panels (log x-axis,
M_DM in Msun; y-axis M_total in Msun). Cliff/known-feature positions are good to
~0.3 dex in M_DM and ~0.2 Msun in M_T. Replace POINTS with the exact tabulated
data (the author's Rom_pOMpDM / RomM_sorted files) when available -- the API
(mtot_of_mdm) stays the same.

Conventions:
  baseline M_T0 = 2.0 Msun (the DM-free NS maximum mass in these figures).
  Each curve: list of (M_DM[Msun], M_T[Msun]) with increasing M_DM.
  Below the first point  -> baseline 2.0 (DM content negligible).
  Above the last point   -> the last tabulated M_T (held).
  Interpolation is linear in log10(M_DM).
"""
from __future__ import annotations

import bisect
import math

BASELINE_MT = 2.0  # Msun, NS DM-free maximum mass in Figs 10/11

# (y, m_chi_GeV) -> [(M_DM, M_T), ...]   NS branch, eyeball-digitized.
POINTS = {
    # ---- Fig. 10 : weakly interacting, y = 1e-1  (flat 2.0 -> cliff to ~0) ----
    (0.1, 500): [(1e-11, 2.00), (5e-7, 2.01), (1e-6, 1.95), (1.5e-6, 1.40),
                 (1.8e-6, 0.80), (2.2e-6, 0.15)],
    (0.1, 200): [(1e-10, 2.00), (5e-6, 2.02), (1e-5, 1.90), (2e-5, 1.20),
                 (3e-5, 0.50), (4e-5, 0.12)],
    (0.1, 100): [(1e-9, 2.00), (3e-5, 2.03), (6e-5, 1.70), (9e-5, 1.00),
                 (1.3e-4, 0.40), (1.7e-4, 0.12)],
    (0.1, 50):  [(1e-8, 2.00), (1e-4, 2.03), (2e-4, 1.60), (3e-4, 0.90),
                 (4.5e-4, 0.30), (6e-4, 0.12)],
    (0.1, 10):  [(1e-7, 2.00), (1e-3, 2.04), (2e-3, 1.50), (3e-3, 0.80),
                 (4e-3, 0.20)],
    (0.1, 5):   [(1e-6, 2.00), (8e-3, 2.05), (1.3e-2, 1.40), (2e-2, 0.60),
                 (2.5e-2, 0.15)],
    (0.1, 1):   [(1e-4, 2.00), (2e-1, 2.05), (3e-1, 2.00), (5e-1, 1.50),
                 (6e-1, 1.00), (8e-1, 0.30)],
    # ---- Fig. 11 : strongly interacting, y = 1e3 ----
    # heavy DM (>=50 GeV) collapses to ~0; light DM (<=10 GeV) grows.
    (1000, 500): [(1e-9, 2.00), (8e-4, 2.00), (1.3e-3, 1.20), (1.8e-3, 0.40),
                  (2e-3, 0.10)],
    (1000, 200): [(1e-8, 2.00), (7e-3, 2.00), (1e-2, 1.10), (1.5e-2, 0.30),
                  (2e-2, 0.10)],
    (1000, 100): [(1e-7, 2.00), (5e-2, 2.00), (8e-2, 1.10), (1.2e-1, 0.30),
                  (1.5e-1, 0.10)],
    (1000, 50):  [(1e-6, 2.00), (1.5e-1, 2.00), (2.5e-1, 1.00), (4e-1, 0.30)],
    (1000, 10):  [(1e-7, 2.00), (1e-1, 2.10), (3e-1, 2.60), (6e-1, 3.50),
                  (1.0, 4.50)],
    (1000, 5):   [(1e-5, 2.00), (1.0, 2.50), (3.0, 5.00), (6.0, 9.00),
                  (9.0, 12.0)],
    (1000, 1):   [(1e-6, 2.00), (1e-5, 2.30), (1e-4, 5.00), (1e-3, 9.00),
                  (5e-3, 13.0), (1e-2, 16.0)],
}

AVAILABLE = sorted(POINTS)


def _nearest_mchi(y, m_chi_GeV):
    masses = sorted(m for (yy, m) in POINTS if yy == y)
    if not masses:
        raise KeyError(f"no digitized curves for y={y}; have y in "
                       f"{sorted({yy for yy, _ in POINTS})}")
    return min(masses, key=lambda m: abs(math.log10(m) - math.log10(m_chi_GeV)))


def mtot_of_mdm(m_dm_Msun, y=0.1, m_chi_GeV=100.0, snap_mchi=True):
    """Total max NS mass [Msun] for accreted DM mass m_dm_Msun.

    y         : DM self-interaction strength (0.1 = Fig.10, 1000 = Fig.11).
    m_chi_GeV : DM particle mass; snapped to the nearest digitized curve.
    Below the curve's first point returns the 2.0 Msun baseline; the curves are
    flat there, so negligible accreted mass leaves M_total unchanged.
    """
    key = (y, m_chi_GeV)
    if key not in POINTS:
        if not snap_mchi:
            raise KeyError(f"no curve {key}; available {AVAILABLE}")
        m_chi_GeV = _nearest_mchi(y, m_chi_GeV)
        key = (y, m_chi_GeV)
    pts = POINTS[key]
    xs = [math.log10(mx) for mx, _ in pts]
    ys = [mt for _, mt in pts]
    if m_dm_Msun <= 0:
        return BASELINE_MT
    x = math.log10(m_dm_Msun)
    if x <= xs[0]:
        return BASELINE_MT
    if x >= xs[-1]:
        return ys[-1]
    i = bisect.bisect_right(xs, x)
    x0, x1, y0, y1 = xs[i - 1], xs[i], ys[i - 1], ys[i]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def _selftest():
    # baseline far below any cliff
    assert abs(mtot_of_mdm(1e-12, 0.1, 100) - 2.0) < 1e-9
    # weak: 100 GeV crosses M_T=1 near ~9e-5
    assert 0.6 < mtot_of_mdm(9e-5, 0.1, 100) < 1.4
    # weak: lighter DM cliffs at larger M_DM than heavier DM
    half = lambda m: next(d for d in (10**(e/10) for e in range(-110, 1))
                          if mtot_of_mdm(d, 0.1, m) < 1.0)
    assert half(10) > half(100) > half(500)
    # strong + light DM grows above baseline
    assert mtot_of_mdm(1e-3, 1000, 1) > 5.0
    # nearest-mchi snap
    assert abs(mtot_of_mdm(1e-12, 0.1, 120) - 2.0) < 1e-9
    print("mtot_vs_mdm_digitized self-test: OK")
    print("  weak  y=0.1  m=100: M_T(3e-5)=%.2f M_T(9e-5)=%.2f M_T(2e-4)=%.2f"
          % (mtot_of_mdm(3e-5, 0.1, 100), mtot_of_mdm(9e-5, 0.1, 100),
             mtot_of_mdm(2e-4, 0.1, 100)))
    print("  strong y=1e3 m=1:   M_T(1e-5)=%.2f M_T(1e-3)=%.2f M_T(1e-2)=%.2f"
          % (mtot_of_mdm(1e-5, 1000, 1), mtot_of_mdm(1e-3, 1000, 1),
             mtot_of_mdm(1e-2, 1000, 1)))


if __name__ == "__main__":
    _selftest()
