# Stage 1 → 2 glue: DM capture → two-fluid TOV

Connects the cross-matcher's per-pulsar local DM density to the two-fluid TOV
solver (`TOV-dark-sector-stars`), so a measured neutron-star mass can be
compared with the mass predicted once captured dark matter is included.

```
cross-matcher                dm_capture.py                tov_bridge.py
pulsar_clump_table.csv  -->   rho_DM, M_NS, age,    -->   f_DM -> facdar -> TOV
(rho_DM at each pulsar)       sigma -> M_DM, f_DM          -> (M_total, R, Lambda)
                                                           -> compare vs observed M_NS
```

## Modules

| File | Role | Needs ROOT build? |
|---|---|---|
| `dm_capture.py` | DM capture physics (Eqs. 9–11 of the seed paper): ρ_DM, σ, age, M, R → accumulated M_DM and f_DM = M_DM/M_NS. Faithful port of `convert_clumpy_draw_to_root.C`. | no |
| `tov_bridge.py` | Generate a solver JSON config, run the binary (if built), parse `rad-mass_phys_*.dat`, pick the star at a target f_DM / mass / M_max. | run: yes · parse/config: no |
| `test_dm_capture.py` | Self-checks: port fidelity vs the C++ constants, linearity, m_χ-independence, required-σ round-trip. | no |

## Capture physics — the key numbers

The accumulated mass is **linear** in (ρ_DM, σ_DM/σ_crit, time) and
**independent of the DM particle mass** (the capture rate ∝ 1/m_χ, the per-capture
mass ∝ m_χ). At the canonical local density:

```
rho = 0.4 GeV/cm^3, sigma_ratio = 1, t = 10 Gyr  ->  M_DM ~ 7e-14 Msun  (f_DM ~ 5e-14)
```

So **clump/field densities produce no measurable mass or Λ change** — confirmed
on the cross-match sweep (0/219 pulsars reach f_DM > 1e-3 at any realization;
max local ρ ≈ 0.016 GeV/cm³). An observable f_DM (≳ 1e-3) requires either
σ_ratio ≳ 1e10 or the **thermalised central spike density** ρ(r_min) ≈
7.7×10⁹ GeV/cm³ (seed paper): at that density 52/219 pulsars cross f_DM > 1e-3
with σ_ratio only ~0.2–0.5. `required_sigma_ratio` reports, per pulsar, the
cross-section needed to reach a target f_DM.

```bash
# clump densities (negligible effect)
python3 pipeline/dm_capture.py --table data/crossmatch/pulsar_clump_table.csv \
    --sigma-ratio 1 --time-Gyr 10 --out data/crossmatch/capture_clump.csv
# thermalised spike density (observable regime)
python3 pipeline/dm_capture.py --table data/crossmatch/pulsar_clump_table.csv \
    --rho 7.7e9 --sigma-ratio 1 --out data/crossmatch/capture_spike.csv
```

## f_DM ↔ facdar calibration

The solver's `facdar` is a DM central-pressure knob, **not** f_DM directly; the
realised fraction is `massDM/massT` in the output. To target a specific f_DM,
scan `facdar` once per EoS and interpolate with `tov_bridge.facdar_for_fdm`. The
bridge always reports the achieved fraction so the calibration stays explicit.

## Running TOV (when built)

This environment has no ROOT/VDT, so `tov_bridge.run_tov` raises with build
guidance; `write_tov_config` and `parse_phys` work regardless. Once
`TOV-dark-sector-stars` is built (see its README), point `$TOV_BIN` at
`build-rel/bin/tov_2k_dbg` (the bridge also auto-detects the sibling build dir).

## Next on `dev/p2-tov-glue`

1. **`predict_pulsar_mass.py` driver** — capture → `facdar` → `run_tov` →
   `star_at_fdm` → predicted (M_total, R, Λ) vs measured M_NS, per pulsar.
2. **facdar calibration grid** — one small `facdar` scan per EoS → f_DM↔facdar table.
3. **MCMC over (m_χ, σ_n)** — wrap the forward model in a likelihood against the
   219 NS-mass measurements (Stage-2 endpoint of the roadmap).
