# clumpy2root

Converts a CLUMPY `.drawn` clump catalogue (ASCII) into a ROOT `TTree` named `clumps`,
one entry per clump.

## Build

    source /opt/root/root-6.38.04/bin/thisroot.sh   # or your ROOT setup
    make

(or in one line, without a Makefile:)

    g++ -O2 -Wall -o clumpy2root clumpy2root.cxx $(root-config --cflags --libs)

## Run

    ./clumpy2root <input.drawn> [-o output.root] [-t treename] [-n maxlines] \
                  [--no-extras] [--rsun R] [--alpha-int A]

Defaults: output = `<input basename>.root`, tree name = `clumps`, all rows.

## Input format

Verified against the catalogues in `/data/53*/*.drawn`. Data rows carry **21**
whitespace-separated fields; the `|` characters and the units line appear only in the
`#`-prefixed header. `prof` is textual (`kEINASTO`), and `Type` is textual (`DSPH`).
A 19-column variant without the leading `Name`/`Type` columns is auto-detected.

## Branches

| branch | type | source column |
|---|---|---|
| `name` | `Long64_t` | `Name` (clump index) |
| `type` | `std::string` | `Type` (e.g. `DSPH`) |
| `long`, `lat` | `Double_t` | galactic longitude / latitude [deg] |
| `d` | `Double_t` | distance [kpc] |
| `z` | `Double_t` | redshift (`-1` when unused) |
| `Rdelta` | `Double_t` | [kpc] |
| `rhos` | `Double_t` | [Msol/kpc^3] |
| `rs` | `Double_t` | [kpc] |
| `prof` | `std::string` | profile enum *name* |
| `prof_code` | `Int_t` | integer code for `prof` (see `prof_codes` in the file) |
| `par1`, `par2`, `par3` | `Double_t` | `#1`, `#2`, `#3` — profile shape parameters |
| `J`, `J_Jcontinuum` | `Double_t` | [GeV^2 cm^-5], ratio |
| `Mdelta`, `Mtid` | `Double_t` | [Msol] |
| `Rtid` | `Double_t` | [kpc] |
| `Mequdens`, `Requdens` | `Double_t` | [Msol], [kpc] |
| `Dgal` | `Double_t` | [kpc] |

`#1/#2/#3` are profile-dependent: for `kEINASTO` `#1` is the shape index alpha
(0.17 in these files) and `#2`/`#3` are unused (`-1`). For `kZHAO` they are
(alpha, beta, gamma). Cut on `prof`/`prof_code` before interpreting them.

## Derived branches (computed on the fly; `--no-extras` disables)

| branch | type | meaning |
|---|---|---|
| `x_gc`, `y_gc`, `z_gc` | `Double_t` | galactocentric cartesian [kpc]; GC at origin, Sun at `(-R_sun,0,0)` |
| `theta_rs`, `theta_delta`, `theta_tid` | `Double_t` | angular radii subtended by `rs`, `Rdelta`, `Rtid` [deg] |
| `Omega_tid` | `Double_t` | solid angle within `Rtid` [sr] |
| `c_delta`, `c_tid` | `Double_t` | `Rdelta/rs`, `Rtid/rs` |
| `f_strip` | `Double_t` | `Mtid/Mdelta` — tidal stripping fraction |
| `rho_mean_tid` | `Double_t` | `3 Mtid / (4pi Rtid^3)` [Msol/kpc^3] |
| `Lann` | `Double_t` | intrinsic annihilation luminosity `int rho^2 dV` [Msol^2/kpc^3], out to `Rtid` (to `Rdelta` when `!tid_valid`) |
| `Lann_delta` | `Double_t` | same, always out to `Rdelta` |
| `J_tot` | `Double_t` | `Lann/d^2` [GeV^2 cm^-5] — **total** J, vs the file's aperture-limited `J` |
| `f_aperture` | `Double_t` | `J / J_tot` — fraction of the clump's J inside the pixel |
| `alpha_int` | `Double_t` | J integration half-angle [deg], parsed from the file header |
| `is_extended` | `Bool_t` | angular radius exceeds `alpha_int` |
| `tid_valid` | `Bool_t` | `Rtid` is a real tidal radius, not the `rs*1e-3` solver floor |

### Conventions, established from the data rather than assumed

* `rs` = r_-2 and `rhos` = rho_-2; `#1` is the Einasto shape index alpha. Re-integrating
  the profile reproduces the file's own `Mdelta` and `Mtid` to its 3 printed digits.
* **CLUMPY's J is `L/d^2`** — integrated over solid angle, with no `1/4pi`. Verified on
  248502 clumps whose tidal radius fits inside the aperture: `J_file/J_tot = 1.0005 +/- 0.005`.
* `R_sun = 8.0 kpc`, recovered from `(long, lat, d)` against the `Dgal` column to 1e-5.
  Override with `--rsun`.
* The analytic Einasto luminosity agrees with direct numerical integration to 6 decimals.
  Non-Einasto profiles (`kZHAO`, `kNFW`) fall back to numerical integration; anything else
  leaves `Lann = J_tot = 0`.

### Caveat: the `Rtid` floor

In these catalogues **14.2%** of clumps (all with `Mdelta` below ~90 Msol) have `Rtid`
pinned at exactly `rs * 1e-3`, the tidal-radius solver's lower bracket rather than a
physical value; their `Mtid` is then ~1e-7 of `Mdelta`. These rows are flagged
`tid_valid == false`. **Always cut on `tid_valid` before using `Rtid`, `Mtid`,
`c_tid`, `f_strip` or `rho_mean_tid`.** Among valid rows, median `f_strip` = 0.60.

The derived doubles compress far worse than the 3-significant-digit source columns, so
the output grows from 66 MB to 217 MB for a 1.76M-clump catalogue. Use `--no-extras`
if you only want the raw columns.

Aliases `glon`/`glat` are set for `long`/`lat` (`long` is a C++ keyword; it works in
`TTree::Draw` but the aliases are handy in compiled code and `TTreeFormula` edge cases).

## Self-describing metadata

The output file explains itself, so a post-mortem years from now needs nothing but
the `.root` file:

* **every branch carries a title** (short description + unit) — visible in `clumps->Print()`
* **`clumps->GetUserInfo()`** is a `TList` of `TNamed(branch, "[unit] brief -- full definition")`,
  the machine-readable data dictionary
* **`README`** is a `TMacro` holding the whole thing as printable text: provenance,
  the conventions that were verified against the data, the known `Rtid` caveat,
  and the full branch dictionary
* **`provenance`** is a `TList` recording the input file, its size and mtime, the exact
  command line, host, user, converter version, ROOT version, conversion time, row
  counts, and the `R_sun` / `alpha_int` actually used
* **`prof_codes`** maps profile names to `prof_code` integers

Read it back with:

    root -l clumps.root
    root [1] README->Print()            // full provenance + dictionary
    root [2] clumps->Print()            // branches with titles and units

or use the bundled macro, which also does single-branch lookup:

    root -l 'describe.C("clumps.root")'
    root -l 'describe.C("clumps.root","J_tot")'

All definitions live in one `kBranchDoc` table in `clumpy2root.cxx` — the same table
sets the branch titles, the UserInfo dictionary and the README, so the three cannot
drift apart. The converter warns if a branch has no dictionary entry.

## Inspect

    root -l clumps.root
    root [1] clumps->Print()
    root [2] clumps->Scan("long:lat:d", "", "", 10)
    root [3] clumps->Draw("J", "J>0")
    root [4] clumps->Draw("Mtid:J", "J>0 && Mtid>0", "COLZ")
    root [5] clumps->Draw("Dgal", "prof==\"kEINASTO\"")
    root [6] clumps->Draw("f_aperture", "is_extended")        // J lost outside the pixel
    root [7] clumps->Draw("f_strip", "tid_valid")             // tidal stripping
    root [8] clumps->Draw("J_tot:Mtid", "tid_valid", "COLZ")

or run the bundled macro:

    root -l 'inspect.C("clumps.root")'

## Error handling

Comment (`#`) and blank lines are skipped. A data row is rejected — with a warning
naming the line number and the offending token, and the run continuing — when the
field count is wrong, a numeric field does not fully parse as a double (`9.71e-`,
`1.14X09`), or a value overflows. The first 20 bad lines are reported individually,
the rest counted. Exit status: 0 = at least one entry written, 2 = I/O error,
3 = no entries written.
