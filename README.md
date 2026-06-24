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
> coordinate/sky-map tooling are in place. The pulsar ↔ clump **3-D
> cross-matcher is now implemented** in
> [`CLUMPY_utilities/crossmatch_pulsars_clumps.py`](CLUMPY_utilities/crossmatch_pulsars_clumps.py)
> (reads the `.drawn` catalogs directly, ranks pulsars by local DM density, and
> vets matches with a longitude-scramble false-association control). See
> [Cross-matching](#cross-matching-clumps--pulsars).

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
catalogs, the CLUMPY stack. Three ways: a Docker container (laptop/desktop), an
Apptainer/Singularity image (HPC clusters), or a local ROOT install.

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

### Option B — HPC with Apptainer / Singularity

HPC clusters don't run Docker (no root, no daemon). `clumpy.def` converts the
**already-built** Docker image into a single portable `.sif` — no recompilation,
runs rootless. Build it once (on a host that has the image and Apptainer):

```bash
apptainer build clumpy.sif clumpy.def
# equivalently, without the def file:
#   apptainer build clumpy.sif docker-daemon://my-clumpy:prod
```

Move it to the cluster (the `.sif` is one self-contained file):

```bash
# simplest: build on your laptop, then copy the single file
scp clumpy.sif user@hpc:~/

# no Docker on the build host? export a tarball and build from it there
docker save my-clumpy:prod -o my-clumpy.tar
apptainer build clumpy.sif docker-archive://my-clumpy.tar
# or pull from a registry
apptainer build clumpy.sif docker://<registry>/my-clumpy:prod
```

Run it like the Docker command — Apptainer **auto-mounts your current directory
(writable)**, so the `-v`/`-w` binds are unnecessary and output lands on the host
filesystem, owned by you:

```bash
apptainer exec clumpy.sif \
  clumpy -g7 -i clumpy_config/clumpy_params_g7.txt \
    --gSIM_SEED=4448 --gPP_DM_MASS_GEV=100 --gSIM_HEALPIX_NSIDE=2048 \
    --gSIM_THETA_ORTH_SIZE_DEG=d --gSIM_THETA_SIZE_DEG=360 \
    --gSIM_OUTPUT_DIR=output
```

| Docker | Apptainer |
|---|---|
| `docker run --rm -v "$PWD":/work -w /work my-clumpy:prod clumpy …` | `apptainer exec clumpy.sif clumpy …` |

> **Full-sky FOV.** For a 360° map use *RING mode* — `--gSIM_THETA_ORTH_SIZE_DEG=d`
> (the literal `d`) with `--gSIM_THETA_SIZE_DEG=360`. Setting both dimensions to a
> number ≥ 180 is rejected by CLUMPY. Use a smaller `--gSIM_HEALPIX_NSIDE` for a
> quick test; full-resolution full-sky J-maps are expensive.

Inside a **Slurm** job, just place the `apptainer exec …` line in your batch
script (prefix with `module load apptainer` if your site requires it). The image
ships OpenMPI; for multi-node MPI, bind the host MPI/PMI at run time.

Ready-made (cluster-neutral) batch scripts live in [`slurm/`](slurm/):

| Script | Purpose |
|---|---|
| [`slurm/clumpy_g7_slurm.sbatch`](slurm/clumpy_g7_slurm.sbatch) | one `-g7` draw as a single SMP job |
| [`slurm/clumpy_g7_array.sbatch`](slurm/clumpy_g7_array.sbatch) | a job **array** sweeping a `seed × {100,500} GeV` grid in parallel |

Both are single-node OpenMP (SMP) jobs. Before submitting, set your SLURM
`--account` (`sacctmgr show assoc -n user=$USER format=Account`), confirm the
Apptainer module name (`module avail apptainer`), and adjust `SIF`/`PROJECT`
paths. Submit with `sbatch slurm/clumpy_g7_slurm.sbatch`; results are written to
`/scratch/$USER/$SLURM_JOB_ID` and staged back to `output/`.

### Option C — local ROOT

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

5. **Cross-match clumps ↔ pulsars** — see [Cross-matching](#cross-matching-clumps--pulsars).

### Cross-matching clumps ↔ pulsars

`CLUMPY_utilities/crossmatch_pulsars_clumps.py` associates each pulsar with the
dark-matter clump(s) it is embedded in, in **3-D** (Galactic position *and*
heliocentric distance). It reads the `.drawn` catalog(s) **directly** — not the
stale ROOT files — and the authoritative NS catalog
(`NSmass_dist2GC_corr_data_extended_ext.h`, Galactic `l,b,d` + NS mass/radius/
age) from the sibling `ADM-NS_GWsignal_templ` repo.

```bash
python3 CLUMPY_utilities/crossmatch_pulsars_clumps.py \
    --drawn data/raw/100GeV/annihil_seed4448_*.drawn \
    --rho-thresh 0.4 --alpha 1.0 --scramble 50 --out data/crossmatch
```

For every pulsar it places both populations in heliocentric Cartesian
coordinates, evaluates the Einasto profile of every clump at the pulsar offset,
and reports the **local DM density** (sum and max over clumps) plus the
density-dominant and geometrically-nearest clumps. Two association criteria:

- **primary (density):** the pulsar sits where some clump's local density
  exceeds `--rho-thresh` [GeV/cm³] — implemented as a per-clump *density radius*
  (the analytic Einasto inverse), so it is exact and cheap to scramble;
- **secondary (tidal):** `d3d < alpha · R_tidal` — a pure geometric membership
  test that *saturates* (big clumps' tidal spheres blanket the inner Galaxy), so
  it is reported only for context.

The **false-association rate** is estimated by re-running the density match on
`--scramble` Galactic-longitude-randomised clump catalogs (the look-elsewhere
control). Outputs in `--out`: `pulsar_clump_table.csv` (one row per pulsar),
`associations.csv` (density-match pairs, densest first), `crossmatch_summary.txt`.

> Flag overrides: `--pulsars PATH` (or `$DMCLUMPS_PSR_CATALOG`) for the NS
> catalog; pass multiple `--drawn` files/globs to combine realizations.

### Configuration (`makeCLUMPY_skyMap.py`)

Input/output locations are resolved from environment variables (repo-relative
defaults), so the script is portable — set any you need before running:

| Variable | Default | Purpose |
|---|---|---|
| `DMCLUMPS_NTUPLE` | `data/root/100GeV/ntuple_seed4448_…root` | ROOT ntuple (converted CLUMPY catalog) |
| `DMCLUMPS_MERGED` | `data/root/mergedALL_…root` | merged NS↔clump analysis ROOT file |
| `DMCLUMPS_SKYMAP_DIR` | `data/skymaps` | dir of auxiliary all-sky maps (WMAP/DIRBE/Fermi) |
| `DMCLUMPS_EXPOSURE_MAP` | `<skymap_dir>/exposure_healpix_full_clean_gt20MeV.fits` | Fermi-LAT exposure map |
| `DMCLUMPS_CLUMPY_FITS_DIR` / `DMCLUMPS_CLUMPY_FITS` | `data/skymaps` / `annihil_…nside2048.fits` | CLUMPY HEALPix `.fits` map |
| `DMCLUMPS_PLOT_DIR` | `plots_v2` (auto-created) | output directory for figures/tables |

```bash
DMCLUMPS_NTUPLE=/data/root/my.root DMCLUMPS_PLOT_DIR=out python CLUMPY_utilities/makeCLUMPY_skyMap.py
```

> Running the full script additionally requires `astrotools`, `gammapy` and
> `NPTFit` (plus the core `numpy/healpy/uproot/astropy/matplotlib` stack), and the
> external FITS inputs above.

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
