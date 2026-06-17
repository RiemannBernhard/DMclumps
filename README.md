# DMclumps

[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C.svg?logo=cplusplus&logoColor=white)](https://isocpp.org/)
[![Python 3](https://img.shields.io/badge/Python-3-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-2311.00113-b31b1b.svg)](https://arxiv.org/abs/2311.00113)

Convert [CLUMPY](https://clumpy.gitlab.io/CLUMPY/) dark-matter sub-halo
("clump") catalogs into [CERN ROOT](https://root.cern/) `TTree`s, so they can be
cross-matched against a pulsar database and used to associate dark-matter clumps
with individual pulsars along their lines of sight.

> **Project status.** The CLUMPY → ROOT conversion and the supporting
> coordinate/sky-map tooling are in place. The pulsar ↔ clump **cross-matcher
> itself is not yet implemented** — `makeCLUMPY_skyMap.py` currently only
> visualizes the catalog. See [Status & caveats](#status--caveats).

---

## Publications

These utilities were developed for, and used to produce the results in:

> **Neutron star mass in dark matter clumps**  
> Maksym Deliyergiyev, Antonino Del Popolo, Morgan Le Delliou  
> *Monthly Notices of the Royal Astronomical Society* **527** (2023) 4483–4504;
> erratum *MNRAS* **531** (2024) 4, 4263–4274.  
> [doi:10.1093/mnras/stad3311](https://doi.org/10.1093/mnras/stad3311) ·
> arXiv:[2311.00113](https://arxiv.org/abs/2311.00113) [astro-ph.GA]

If you use this code, please cite the paper above:

```bibtex
@article{Deliyergiyev2023_NSinDMclumps,
  author        = {Deliyergiyev, Maksym and Del Popolo, Antonino and Le Delliou, Morgan},
  title         = {Neutron star mass in dark matter clumps},
  journal       = {Monthly Notices of the Royal Astronomical Society},
  volume        = {527},
  number        = {3},
  pages         = {4483--4504},
  year          = {2023},
  doi           = {10.1093/mnras/stad3311},
  eprint        = {2311.00113},
  archivePrefix = {arXiv},
  primaryClass  = {astro-ph.GA},
  note          = {Erratum: MNRAS 531 (2024) 4, 4263--4274}
}
```

---

## Scientific goal

CLUMPY produces Monte-Carlo realizations of the Galactic dark-matter sub-halo
population: for each drawn clump it gives a sky position (Galactic `l, b`), a
distance, a density profile, masses, radii, and an annihilation `J`-factor. The
aim of this repository is to turn those catalogs into a queryable ROOT format
and combine them with a list of observed pulsars, in order to test whether a
given pulsar is spatially coincident with — and plausibly associated to — one or
more dark-matter clumps. Any such association must be done in **3-D** (sky
position *and* distance) and vetted statistically; see the caveats below.

---

## Repository layout

```
CLUMPY_utilities/
  convert_clumpy_draw_to_root.C            .drawn -> ROOT TTree   (primary, FLAG=1)
  convert_clumpy_draw_annih_rs_to_root.C   halo-centre variant    (FLAG=2, see caveats)
  convert_RADEC_to_Galcoord.py             pulsar RA/Dec -> Galactic l, b
  makeCLUMPY_skyMap.py                      HEALPix/Mollweide sky maps of the catalog
data/
  raw/{100GeV,500GeV}/*.drawn              CLUMPY clump catalogs (input)
  root/{100GeV,500GeV}/*.root              converted ROOT TTrees    (see Status)
Dockerfile.clumpy                          ROOT + GreAT + CLUMPY build/run environment
.dockerignore
```

---

## Environment

You need CERN ROOT (≥ 6.x; tested against 6.38/6.39) plus, to *generate* new
catalogs, the CLUMPY stack. Two options:

### Option A — container (recommended, reproducible)

`Dockerfile.clumpy` builds a self-contained image with **ROOT 6.39 + GreAT +
CLUMPY 3.1.1** on Rocky Linux 10. Versions are pinned via build args
(`ROCKY_TAG`, `ROOT_VERSION`, `CLUMPY_VERSION`); an optional scientific-Python
layer (`WITH_PYSCI=true`, on by default) installs `numpy/scipy/astropy/healpy/
uproot/pandas/matplotlib` for the analysis scripts.

```bash
# Build (compiles ROOT from source; first build is long)
docker build -f Dockerfile.clumpy -t clumpy:v3.1.1 .

# Run, mounting the repo at /work
docker run --rm -it -v "$PWD":/work -w /work clumpy:v3.1.1

# CLUMPY-only image (skip the Python layer)
docker build -f Dockerfile.clumpy --build-arg WITH_PYSCI=false -t clumpy:v3.1.1-core .
```

Inside the container `clumpy`, `root`, and the GreAT libraries are on the path
(`$CLUMPY`, `$ROOTSYS`, `$GREAT` are set).

### Option B — local ROOT

```bash
source /path/to/root/bin/thisroot.sh   # sets ROOTSYS / PATH / LD_LIBRARY_PATH
root-config --version
```

The Python helpers additionally require `astropy`, `healpy`, `uproot`,
`numpy`, `pandas` (and `matplotlib` for plotting).

---

## Pipeline / usage

1. **Generate a clump catalog** with CLUMPY (in the container) → produces a
   `*.drawn` text file. *(Smoke-test only: `clumpy -h`. A full draw can take a
   very long time and is not required just to validate the install.)*

2. **Convert `.drawn` → ROOT.** The macro takes input and (optional) output
   paths; if the output is omitted it is derived from the input basename:

   ```bash
   root -l -b -q 'CLUMPY_utilities/convert_clumpy_draw_to_root.C(
       "data/raw/100GeV/annihil_seed4448_..._nside2048.drawn",
       "data/root/100GeV/ntuple_annihil_seed4448_..._nside2048.root")'
   ```

   This writes a lean science file (`TTree` `CLUMPY_output` + provenance) and a
   separate `*_diagnostics.root` with ~200 QA histograms/profiles.

3. **Pulsar coordinates.** `convert_RADEC_to_Galcoord.py` converts a list of
   pulsar equatorial (RA/Dec, ICRS) positions to Galactic `l, b`.

4. **Sky maps.** `makeCLUMPY_skyMap.py` reads the TTree via `uproot` and draws
   HEALPix/Mollweide maps of the clumps with the pulsar list overlaid.

5. **Cross-match clumps ↔ pulsars** — *not yet implemented* (see caveats).

---

## Data formats

### Input — CLUMPY `.drawn`

Plain text, one clump per non-`#` line, whitespace-separated, **21 columns**
(Galactic coordinates):

| # | Column | Unit | # | Column | Unit |
|--:|--------|------|--:|--------|------|
| 1 | ID | – | 12 | profile param #2 | – |
| 2 | type | – | 13 | profile param #3 | – |
| 3 | longitude `l` | deg | 14 | `J` (0.0299° aperture) | GeV²·cm⁻⁵ |
| 4 | latitude `b` | deg | 15 | `J / J_continuum` | – |
| 5 | distance `d` | kpc | 16 | M<sub>Δ</sub> | M<sub>☉</sub> |
| 6 | redshift `z` | – | 17 | M<sub>tidal</sub> | M<sub>☉</sub> |
| 7 | R<sub>Δ</sub> | kpc | 18 | R<sub>tidal</sub> | kpc |
| 8 | ρ<sub>s</sub> | M<sub>☉</sub>·kpc⁻³ | 19 | M<sub>equidens</sub> | M<sub>☉</sub> |
| 9 | scale radius `rs` | kpc | 20 | R<sub>equidens</sub> | kpc |
| 10 | profile (e.g. kEINASTO) | – | 21 | D<sub>gal</sub> | kpc |
| 11 | profile param #1 | – | | | |

> `J` is integrated in a **fixed 0.0299° aperture**, *not* the total clump
> `J`-factor — don't use it as a total.

### Output — ROOT `TTree` `CLUMPY_output`

One entry per clump. Coordinates are Galactic. The branches most relevant to the
pulsar association are:

| Branch | Type | Unit | Meaning |
|--------|------|------|---------|
| `DMclumpID` | int | – | clump index (1…N) |
| `GalacLongitude_deg` | float | deg | Galactic longitude `l` |
| `GalacLatitude_deg` | float | deg | Galactic latitude `b` |
| `Distance_kpc` | float | kpc | heliocentric distance |
| `DM_Dgal_kpc` | float | kpc | galactocentric distance |
| `AngRtidal_deg` | float | deg | on-sky angular radius `asin(R_tidal/d)` |
| `AngRdelta_deg` | float | deg | on-sky angular radius `asin(R_Δ/d)` |
| `DM_Mtidal_Msol` / `DM_Mdelta_Msol` | float | M<sub>☉</sub> | tidal / Δ mass |
| `DM_Rtidal_kpc` / `Rdelta_kpc` | float | kpc | tidal / Δ radius |
| `DM_concentration` | float | – | R<sub>Δ</sub>/r<sub>s</sub> |
| `DMprof_rhos_Msolkpc3` / `DMprof_ScaleRadius_rs_kpc` | float | – | Einasto ρ<sub>s</sub>, r<sub>s</sub> |
| `DM_Jfactor_GeV2cm5` | float | GeV²·cm⁻⁵ | `J` (0.0299° aperture) |

Use **`AngRtidal_deg` / `AngRdelta_deg` as the on-sky matching tolerance** — not
`alpha_S_deg`, which is the much smaller *scale-radius* angle `asin(rs/d)`.
Many additional derived branches (density-profile samples `DMprof_rhoDM_*`,
accreted-mass/rate columns, `*_over_rS` ratios) are also present. The file also
holds provenance objects `clumpy_source_drawn` and `clumpy_Jfactor_aperture`
(`TNamed`), and the QA histograms live in the sibling `*_diagnostics.root`.

---

## Status & caveats

- **The committed `data/root/*.root` are stale.** They were produced by an
  earlier, buggy version of the converter (it parsed only half the clumps and
  swapped `l`/`b`). **Regenerate them** by re-running
  `convert_clumpy_draw_to_root.C` over `data/raw/**/*.drawn`.
- **`convert_clumpy_draw_annih_rs_to_root.C` (FLAG=2) is documented but not
  fixed.** It expects a different `.drawn` column layout (coordinates w.r.t. the
  halo centre) that is not among the committed inputs and could not be
  validated; its header notes the known issues and points to the primary
  converter as the porting template.
- **Statistical vetting is required before any "pulsar X ↔ clump Y" claim:**
  - match in 3-D (`Δθ` *and* a consistent distance), since a foreground/
    background clump can sit on a sightline without being physically associated;
  - account for the look-elsewhere effect over N<sub>pulsar</sub> × N<sub>clump</sub>
    pairs and report a **false-association rate**, e.g. by scrambling the clump
    catalog in `l` and re-running the identical matcher.

---

## Credits

- **CLUMPY** — Charbonnier, Combet & Maurin, *Comput. Phys. Commun.* (clump
  catalogs and `J`-factors). <https://clumpy.gitlab.io/CLUMPY/>
- **ROOT** — Brun & Rademakers et al., CERN. <https://root.cern/>
- **GreAT** — Grenoble Analysis Toolkit (MCMC/Jeans).

## License

Released under the [MIT License](LICENSE). If you use this code in academic
work, please also cite the paper listed under [Publications](#publications).
