#!/usr/bin/env python3
"""Self-checks for dm_capture: port fidelity + internal consistency.

Run: python3 pipeline/test_dm_capture.py   (exit 0 = all pass)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dm_capture as c


def approx(a, b, rel=1e-9):
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-300)


def test_accreted_mass_matches_cpp_constant():
    # Reproduce the C++ accreted_mass(0.4, 1.0, 10) by hand.
    GeV2Msun = (1.0 / 1.111172751) * 1e-57
    expect = 1.3e43 * (0.4 / 0.3) * 10.0 * (0.45 * 1.0) * GeV2Msun
    got = c.accreted_mass(0.4, 1.0, 10.0)
    assert approx(got, expect), (got, expect)
    assert approx(got, 7.0196e-14, rel=1e-3), got  # sanity magnitude


def test_mass_is_linear_in_sigma_and_time_and_rho():
    base = c.accreted_mass(0.4, 1.0, 10.0)
    assert approx(c.accreted_mass(0.4, 2.0, 10.0), 2 * base)
    assert approx(c.accreted_mass(0.4, 1.0, 20.0), 2 * base)
    assert approx(c.accreted_mass(0.8, 1.0, 10.0), 2 * base)


def test_rate_route_is_mchi_independent():
    # Accumulated mass via the rate route must not depend on m_chi.
    a = c.accreted_mass_from_rate(0.4, 1.0, 1.0, 1.4, 11.0 / c.RSUN_KM, 10.0)
    b = c.accreted_mass_from_rate(0.4, 1.0, 500.0, 1.4, 11.0 / c.RSUN_KM, 10.0)
    assert approx(a, b), (a, b)
    assert a > 0


def test_required_sigma_inverts_accreted_mass():
    # required_sigma_ratio must reproduce a known f_DM round-trip.
    rho, M, t, target = 0.4, 1.4, 10.0, 1e-3
    sig = c.required_sigma_ratio(target, rho, M, t)
    f = c.dm_fraction(c.accreted_mass(rho, sig, t), M)
    assert approx(f, target, rel=1e-9), (f, target)


def test_canonical_fdm_is_negligible_but_spike_is_not():
    f_field = c.dm_fraction(c.accreted_mass(0.4, 1.0, 10.0), 1.4)
    f_spike = c.dm_fraction(c.accreted_mass(7.7e9, 1.0, 10.0), 1.4)
    assert f_field < 1e-13, f_field
    assert f_spike > f_field * 1e9 * 0.5, (f_field, f_spike)


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nall {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
