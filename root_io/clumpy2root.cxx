// clumpy2root.cxx -- convert a CLUMPY ".drawn" clump catalogue into a ROOT TTree.
//
// Input format (as produced by CLUMPY, verified against the files in /data/53*/):
//
//   # This file lists the properties and Jcl(centre) for all drawn clumps
//   # Name  Type  long  lat  d  z  Rdelta | rhos rs prof. #1 #2 #3 | J J/Jcontinuum | ...
//   #  -      -   [deg] [deg] [kpc] - [kpc] | [Msol/kpc3] [kpc] [enum] - - - | ...
//   #
//         1  DSPH  -92.23  -0.06  2.87e-03 -1 6.77e-04  1.19e+09 5.66e-06  kEINASTO
//           0.17  -1  -1  9.63e+15 3.28e+00  3.31e-05 3.79e-12 5.66e-09 3.40e-05 5.66e-03 8.00e+00
//
// => 21 whitespace-separated fields per data row; the '|' characters appear only in
//    the comment header, never in the data. 'prof' is TEXTUAL (an enum *name*), so it
//    is stored both as a std::string branch and as a small integer code branch.
//
// A 19-column variant (no leading Name/Type columns) is also accepted and auto-detected.
//
// Usage: clumpy2root <input.drawn> [-o output.root] [-t treename] [-n maxlines]

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TNamed.h"
#include "TString.h"
#include "TObjString.h"
#include "TMacro.h"
#include "TDatime.h"
#include "TSystem.h"
#include "TROOT.h"
#include <sstream>
#include "TList.h"
#include "TTree.h"
#include "TMath.h"

namespace {

const char *const kVersion = "1.1";

// ---------------------------------------------------------------------------
// Profile-name -> integer code.  The first block mirrors CLUMPY's own gENUM_PROFILE
// ordering so that codes are stable across files; anything unknown is assigned a
// fresh code starting at 100 and reported on stderr.
// ---------------------------------------------------------------------------
class ProfileCodes {
public:
   ProfileCodes()
   {
      const char *known[] = {"kZHAO",     "kEINASTO", "kEINASTO_N", "kNFW",
                             "kISHIYAMA", "kBURKERT", "kPSEUDOISO", "kDPDZ_HOST",
                             "kPLUMMER",  "kPOINT"};
      for (int i = 0; i < static_cast<int>(sizeof(known) / sizeof(known[0])); ++i)
         fCodes[known[i]] = i;
   }

   int Code(const std::string &name)
   {
      std::map<std::string, int>::const_iterator it = fCodes.find(name);
      if (it != fCodes.end()) return it->second;
      const int code = fNextDynamic++;
      fCodes[name] = code;
      std::cerr << "clumpy2root: note: unknown profile name '" << name
                << "', assigning code " << code << "\n";
      return code;
   }

   // Only the profiles actually seen in the file are worth persisting.
   void Note(const std::string &name, int code) { fSeen[name] = code; }
   const std::map<std::string, int> &Seen() const { return fSeen; }

private:
   std::map<std::string, int> fCodes;
   std::map<std::string, int> fSeen;
   int fNextDynamic = 100;
};

// Split on whitespace, dropping the '|' column separators should a CLUMPY build
// ever emit them in the data rows.
void Tokenize(const std::string &line, std::vector<std::string> &out)
{
   out.clear();
   const char *p = line.c_str();
   while (*p) {
      while (*p && std::isspace(static_cast<unsigned char>(*p))) ++p;
      if (!*p) break;
      const char *start = p;
      while (*p && !std::isspace(static_cast<unsigned char>(*p))) ++p;
      std::string tok(start, p - start);
      if (tok != "|") out.push_back(tok);
   }
}

// strtod with full-token consumption: rejects "1.0e", "nan_", "3x", empty, ...
// strtod is correctly rounded, so the shortest double round-tripping the decimal
// literal is produced -- the printed precision of the catalogue is preserved exactly.
bool ToDouble(const std::string &tok, double &value)
{
   if (tok.empty()) return false;
   errno = 0;
   char *end = nullptr;
   const double v = std::strtod(tok.c_str(), &end);
   if (end != tok.c_str() + tok.size()) return false;
   if (errno == ERANGE && v != 0.0) return false; // overflow; underflow to 0 is fine
   value = v;
   return true;
}

bool ToLong64(const std::string &tok, Long64_t &value)
{
   if (tok.empty()) return false;
   errno = 0;
   char *end = nullptr;
   const long long v = std::strtoll(tok.c_str(), &end, 10);
   if (end != tok.c_str() + tok.size() || errno == ERANGE) return false;
   value = static_cast<Long64_t>(v);
   return true;
}

bool IsCommentOrBlank(const std::string &line)
{
   for (std::string::const_iterator c = line.begin(); c != line.end(); ++c) {
      if (std::isspace(static_cast<unsigned char>(*c))) continue;
      return *c == '#';
   }
   return true; // all whitespace
}


// ---------------------------------------------------------------------------
// Derived quantities.
//
// Conventions verified numerically against the catalogues themselves:
//   * rs is r_{-2} and rhos is rho_{-2}; column #1 is the Einasto shape index
//     alpha. Re-integrating the profile reproduces the file's Mdelta and Mtid
//     to the 3 significant digits the file prints.
//   * CLUMPY's J is  L/d^2  (an integral over solid angle, NO 1/4pi), confirmed
//     on 248k clumps whose tidal radius is smaller than the J aperture:
//     median L/d^2 / J_file = 0.9996.
//   * R_sun = 8.0 kpc, recovered from (long, lat, d) vs the Dgal column to 1e-5.
// ---------------------------------------------------------------------------

const double kMsolInGeV = 1.115702e57;          // M_sun c^2 [GeV]
const double kKpcInCm   = 3.0856775814913673e21;
// (Msol/kpc^3)^2 * kpc  ->  GeV^2 cm^-5
const double kLtoJ      = kMsolInGeV * kMsolInGeV / (kKpcInCm * kKpcInCm * kKpcInCm *
                                                     kKpcInCm * kKpcInCm);

// Annihilation luminosity  L = \int_0^R rho^2 dV  [Msol^2 kpc^-3] for Einasto,
// computed in logs so the huge/tiny intermediate factors cannot overflow:
//   L = 4pi rhos^2 rs^3 e^{4/a} (1/a) (a/4)^{3/a} Gamma(3/a) P(3/a, (4/a)(R/rs)^a)
double EinastoLuminosity(double rhos, double rs, double a, double R)
{
   if (!(rhos > 0.) || !(rs > 0.) || !(a > 0.) || !(R > 0.)) return 0.;
   const double s   = 3.0 / a;
   const double u   = (4.0 / a) * std::pow(R / rs, a);
   const double lnL = std::log(4.0 * TMath::Pi()) + 2.0 * std::log(rhos) +
                      3.0 * std::log(rs) + 4.0 / a - std::log(a) +
                      s * std::log(a / 4.0) + TMath::LnGamma(s);
   return std::exp(lnL) * TMath::Gamma(s, u);   // TMath::Gamma(a,x) = P(a,x)
}

// Generic fallback: log-spaced Simpson integration of rho^2 dV, used for any
// profile that is not Einasto. rho() is supplied by the caller.
template <class Rho>
double NumericLuminosity(Rho rho, double R, int n = 2000)
{
   if (!(R > 0.)) return 0.;
   const double lo = std::log(R) - 30.0, hi = std::log(R);
   const double h  = (hi - lo) / n;
   double sum = 0.;
   for (int i = 0; i <= n; ++i) {
      const double lr = lo + i * h, r = std::exp(lr), f = rho(r);
      const double w  = (i == 0 || i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
      sum += w * 4.0 * TMath::Pi() * r * r * r * f * f;  // extra r from dr = r dlnr
   }
   return sum * h / 3.0;
}

// Zhao / generalised-NFW: rho = rhos / ( (r/rs)^g (1+(r/rs)^A)^((B-g)/A) )
double ZhaoRho(double r, double rhos, double rs, double A, double B, double g)
{
   const double x = r / rs;
   return rhos / (std::pow(x, g) * std::pow(1.0 + std::pow(x, A), (B - g) / A));
}



// Wrap `text` to `width` columns. `first` prefixes the first produced line and
// `cont` every continuation line, so bullet markers are not repeated.
std::string Wrap(const std::string &text, size_t width, const std::string &first,
                 const std::string &cont)
{
   std::istringstream is(text);
   std::string word, line, out;
   bool isFirst = true;
   while (is >> word) {
      if (!line.empty() && line.size() + 1 + word.size() > width) {
         out += (isFirst ? first : cont) + line + "\n";
         isFirst = false;
         line.clear();
      }
      line += (line.empty() ? "" : " ") + word;
   }
   if (!line.empty()) out += (isFirst ? first : cont) + line + "\n";
   return out;
}

// ---------------------------------------------------------------------------
// DATA DICTIONARY
//
// One row per branch: { branch, unit, description }. This table is the single
// source of truth: it sets the TBranch titles, it is written into the output
// file as the tree's UserInfo list, and it is rendered into the "README"
// TObjString stored alongside the tree. Edit definitions HERE and nowhere else.
//
// "VERIFIED" below means the statement was checked numerically against these
// very catalogues (see the conventions block above main()), not taken on trust
// from the column name.
// ---------------------------------------------------------------------------
struct BranchDoc { const char *name, *unit, *brief, *desc; };

const BranchDoc kBranchDoc[] = {

// ===== columns copied straight from the CLUMPY .drawn file ==================
{"name", "-", "catalogue index (bookkeeping only)",
 "Sequential index of the clump in the CLUMPY catalogue (1-based; the file's "
 "'Name' column). Bookkeeping only, no physical meaning."},
{"type", "-", "CLUMPY object class label",
 "CLUMPY class label of the drawn object, e.g. 'DSPH' = dwarf-spheroidal-like "
 "subhalo drawn from the Galactic subhalo distribution. Constant in these runs."},
{"long", "deg", "galactic longitude l",
 "Galactic longitude l of the clump centre, in (-180,+180]. Alias: glon."},
{"lat", "deg", "galactic latitude b",
 "Galactic latitude b of the clump centre, in [-90,+90]. Alias: glat."},
{"d", "kpc", "heliocentric distance",
 "Heliocentric distance to the clump centre."},
{"z", "-", "redshift (-1 = unused here)",
 "Redshift. Fixed at -1 in these Galactic (non-cosmological) runs, i.e. unused."},
{"Rdelta", "kpc", "overdensity radius R_200c",
 "Overdensity ('virial') radius. VERIFIED to be R_200c: the mean density inside "
 "Rdelta is 2.55e4 Msol/kpc^3 with only 0.4% scatter, i.e. 200 x rho_crit for "
 "h = 0.678 (Planck). Mdelta is the mass inside it."},
{"rhos", "Msol/kpc^3", "profile normalisation rho_-2",
 "Density normalisation of the profile. For Einasto this is rho_{-2}, the density "
 "at the radius where dln(rho)/dln(r) = -2. VERIFIED by re-integrating the profile "
 "and recovering the file's own Mdelta."},
{"rs", "kpc", "profile scale radius r_-2",
 "Profile scale radius; for Einasto this is r_{-2} (where the log slope is -2). "
 "NOT the tidal or virial radius."},
{"prof", "enum name", "density-profile family name",
 "Density-profile family, as the literal CLUMPY enum name (e.g. 'kEINASTO'). "
 "Determines how par1/par2/par3 must be interpreted."},
{"prof_code", "-", "integer encoding of prof",
 "Integer encoding of 'prof' added by this converter (kZHAO=0, kEINASTO=1, "
 "kEINASTO_N=2, kNFW=3, ...; unknown names get 100+). The full name<->code map is "
 "stored in this file as the TList 'prof_codes'. Use it for fast cuts."},
{"par1", "-", "shape param #1 (Einasto alpha)",
 "CLUMPY column #1. MEANING DEPENDS ON 'prof'. Einasto: the shape index alpha "
 "(0.17 throughout these files), in rho = rhos*exp(-2/alpha*((r/rs)^alpha - 1)). "
 "Zhao/gNFW: alpha. ALWAYS cut on prof before interpreting this column."},
{"par2", "-", "shape param #2 (unused for Einasto)",
 "CLUMPY column #2. Einasto: unused, written as -1. Zhao/gNFW: beta (outer slope) "
 "in rho = rhos/[(r/rs)^gamma * (1+(r/rs)^alpha)^((beta-gamma)/alpha)]."},
{"par3", "-", "shape param #3 (unused for Einasto)",
 "CLUMPY column #3. Einasto: unused, written as -1. Zhao/gNFW: gamma (inner slope)."},
{"J", "GeV^2 cm^-5", "annihilation J-factor inside alpha_int (aperture-limited)",
 "Annihilation J-factor toward the clump centre: J = int dOmega int dl rho^2, "
 "integrated over a cone of half-angle alpha_int (see that branch). APERTURE-LIMITED: "
 "for a clump wider than the cone this is only PART of its emission. Compare J_tot "
 "and f_aperture. (Annihilation scales as rho^2; a decay run would tabulate rho instead.)"},
{"J_Jcontinuum", "-", "clump J over smooth-halo J in the same pixel",
 "Ratio of this clump's J to the J of the smooth ('continuum') Galactic halo in the "
 "same direction and aperture: a contrast / detectability measure, >1 meaning the "
 "clump outshines the diffuse halo in that pixel. Taken from the CLUMPY header; "
 "not independently re-derived here."},
{"Mdelta", "Msol", "mass within Rdelta",
 "Mass enclosed within Rdelta. VERIFIED = M(<Rdelta) from the profile to 0.4%, the "
 "file's own printing precision."},
{"Mtid", "Msol", "bound mass within Rtid (see tid_valid)",
 "Bound mass within the tidal radius, i.e. what survives tidal stripping by the host. "
 "VERIFIED = M(<Rtid) to 0.4% on rows where tid_valid is true. SEE tid_valid."},
{"Rtid", "kpc", "tidal radius (see tid_valid)",
 "Tidal radius of the clump. SEE tid_valid: in 14.2% of rows this is a solver "
 "artefact, not a physical radius."},
{"Mequdens", "Msol", "mass within Requdens",
 "Mass within Requdens. VERIFIED = M(<Requdens) to 0.4%. Median 0.61 x Mdelta."},
{"Requdens", "kpc", "equal-density radius",
 "'Equal-density' radius: where the clump's local density matches the smooth host-halo "
 "density at the clump's position -- an alternative outer edge to Rtid. Consistent with "
 "the data (rho(Requdens) tracks Dgal, ~2e6 Msol/kpc^3 at Dgal = 20 kpc) but not exactly "
 "reproducible here because CLUMPY's host-halo parameters are not in the file. "
 "Typically 0.87 x Rtid and 0.17 x Rdelta."},
{"Dgal", "kpc", "galactocentric distance",
 "Galactocentric distance of the clump. VERIFIED consistent with (long,lat,d) for "
 "R_sun = 8.0 kpc."},

// ===== derived by this converter (--no-extras suppresses all of these) ======
{"x_gc", "kpc", "galactocentric cartesian x",
 "Galactocentric Cartesian x. Frame: Galactic centre at the origin, Sun at "
 "(-R_sun,0,0), +x pointing from the Sun toward the GC, +y toward l=+90 deg, "
 "+z toward b=+90 deg. VERIFIED: sqrt(x^2+y^2+z^2)/Dgal = 1.000004 (rms 1.9e-3, "
 "which is just the file's 3-digit rounding)."},
{"y_gc", "kpc", "galactocentric cartesian y",
 "Galactocentric Cartesian y; see x_gc for the frame definition."},
{"z_gc", "kpc", "galactocentric cartesian z",
 "Galactocentric Cartesian z; see x_gc for the frame definition."},
{"theta_rs", "deg", "angular radius of rs",
 "Angular radius subtended by the scale radius: atan(rs/d)."},
{"theta_delta", "deg", "angular radius of Rdelta",
 "Angular radius subtended by Rdelta: atan(Rdelta/d). The clump's apparent size on "
 "the sky out to its virial boundary."},
{"theta_tid", "deg", "angular radius of Rtid (see tid_valid)",
 "Angular radius subtended by the tidal radius: atan(Rtid/d). SEE tid_valid."},
{"Omega_tid", "sr", "solid angle within Rtid (see tid_valid)",
 "Solid angle within the tidal radius: 2*pi*(1 - cos(theta_tid)). SEE tid_valid."},
{"c_delta", "-", "concentration Rdelta over rs",
 "Concentration Rdelta/rs (= R_200c / r_-2). Median ~96 for these microhaloes."},
{"c_tid", "-", "Rtid over rs (see tid_valid)",
 "Rtid/rs. SEE tid_valid -- this is exactly 1e-3 for the flagged rows."},
{"f_strip", "-", "bound fraction Mtid over Mdelta (see tid_valid)",
 "Mtid/Mdelta: the bound mass fraction surviving tidal stripping. Median 0.60 over "
 "tid_valid rows; meaningless (~1e-7) otherwise. SEE tid_valid."},
{"rho_mean_tid", "Msol/kpc^3", "mean density within Rtid (see tid_valid)",
 "Mean density within the tidal radius: 3*Mtid/(4*pi*Rtid^3). SEE tid_valid."},
{"Lann", "Msol^2/kpc^3", "intrinsic annihilation luminosity to Rtid",
 "Intrinsic annihilation luminosity L = int_0^R rho^2 dV, taken to R = Rtid, or to "
 "Rdelta when tid_valid is false. Unlike J this is DISTANCE-INDEPENDENT and "
 "aperture-free -- the property of the clump itself. Einasto uses a closed form built "
 "on the lower incomplete gamma function, VERIFIED against direct numerical "
 "integration to 6 decimal places."},
{"Lann_delta", "Msol^2/kpc^3", "intrinsic annihilation luminosity to Rdelta",
 "As Lann but always integrated out to Rdelta, so it is unaffected by the Rtid "
 "problem and safe to use on every row."},
{"J_tot", "GeV^2 cm^-5", "total aperture-free J-factor Lann over d^2",
 "Total, aperture-free J-factor of the clump: Lann/d^2. VERIFIED: CLUMPY's J is "
 "normalised as L/d^2 with NO factor 1/(4*pi) -- median J/J_tot = 1.0005 +- 0.005 "
 "over the 248502 clumps small enough to fit inside the aperture."},
{"f_aperture", "-", "fraction of total J inside the file aperture",
 "J/J_tot: the fraction of the clump's total J that the file's aperture actually "
 "captured. <= 1 by construction; drops to ~0.7 for the nearest, most extended "
 "microhaloes. Use it to correct J, or use J_tot directly."},
{"alpha_int", "deg", "J integration half-angle of this file",
 "Half-angle of the cone over which the file's J column was integrated, parsed from "
 "the CLUMPY header. Constant per file: 0.0299 deg at nside=2048, 3.7824 deg at "
 "nside=16 (it is the HEALPix pixel scale of the run)."},
{"is_extended", "bool", "clump is wider than alpha_int",
 "True when the clump's angular radius exceeds alpha_int, i.e. its tabulated J "
 "undercounts its emission. True for 86% of clumps in these catalogues."},
{"tid_valid", "bool", "Rtid is physical, not the solver floor",
 "FALSE when Rtid sits exactly at rs*1e-3 -- the lower bracket of CLUMPY's "
 "tidal-radius solver rather than a converged physical radius. This affects 14.2% of "
 "rows (250345/1762617), all with Mdelta below ~90 Msol, and their Mtid comes out "
 "~1e-7 x Mdelta. ALWAYS require tid_valid before using Rtid, Mtid, c_tid, f_strip, "
 "rho_mean_tid, theta_tid or Omega_tid."},
};
const int kNBranchDoc = static_cast<int>(sizeof(kBranchDoc) / sizeof(kBranchDoc[0]));

const BranchDoc *FindDoc(const char *n)
{
   for (int i = 0; i < kNBranchDoc; ++i)
      if (std::strcmp(kBranchDoc[i].name, n) == 0) return &kBranchDoc[i];
   return nullptr;
}

void Usage(const char *argv0)
{
   std::cerr << "Usage: " << argv0 << " <input.drawn> [-o output.root] [-t treename]"
             << " [-n maxlines]\n"
             << "            [--no-extras] [--rsun R] [--alpha-int A]\n"
             << "  -o  output ROOT file   (default: <input basename>.root)\n"
             << "  -t  TTree name         (default: clumps)\n"
             << "  -n  stop after N data rows (0 = all, the default)\n"
             << "  --no-extras   do not compute the derived columns\n"
             << "  --rsun R      Sun's galactocentric radius in kpc (default 8.0)\n"
             << "  --alpha-int A J integration half-angle in deg (default: from header)\n";
}

} // namespace

int main(int argc, char **argv)
{
   std::string inPath, outPath, treeName = "clumps";
   Long64_t maxLines = 0;
   bool   doExtras = true;
   double rSun     = 8.0;    // kpc; matches the Dgal column to 1e-5
   double alphaInt = -1.0;   // deg; <0 => take it from the file header

   for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      if (a == "-h" || a == "--help") { Usage(argv[0]); return 0; }
      else if (a == "-o" && i + 1 < argc) outPath = argv[++i];
      else if (a == "-t" && i + 1 < argc) treeName = argv[++i];
      else if (a == "-n" && i + 1 < argc) maxLines = std::atoll(argv[++i]);
      else if (a == "--no-extras") doExtras = false;
      else if (a == "--rsun" && i + 1 < argc) rSun = std::atof(argv[++i]);
      else if (a == "--alpha-int" && i + 1 < argc) alphaInt = std::atof(argv[++i]);
      else if (!a.empty() && a[0] == '-') { Usage(argv[0]); return 1; }
      else if (inPath.empty()) inPath = a;
      else { Usage(argv[0]); return 1; }
   }
   if (inPath.empty()) { Usage(argv[0]); return 1; }

   if (outPath.empty()) {
      outPath = inPath;
      const size_t slash = outPath.find_last_of('/');
      if (slash != std::string::npos) outPath = outPath.substr(slash + 1);
      const size_t dot = outPath.find_last_of('.');
      if (dot != std::string::npos && dot != 0) outPath = outPath.substr(0, dot);
      outPath += ".root";
   }

   std::ifstream in(inPath.c_str());
   if (!in) {
      std::cerr << "clumpy2root: ERROR: cannot open input file '" << inPath << "'\n";
      return 2;
   }

   TFile fout(outPath.c_str(), "RECREATE");
   if (fout.IsZombie()) {
      std::cerr << "clumpy2root: ERROR: cannot create ROOT file '" << outPath << "'\n";
      return 2;
   }
   fout.SetCompressionLevel(5);

   TTree *tree = new TTree(treeName.c_str(), "CLUMPY drawn clump catalogue");

   // ---- branch buffers ----------------------------------------------------
   Long64_t     b_name = -1;                 // clump index from the "Name" column
   std::string  b_type;                      // "DSPH", ... (empty in 19-column files)
   Double_t     b_long = 0., b_lat = 0., b_d = 0., b_z = 0., b_Rdelta = 0.;
   Double_t     b_rhos = 0., b_rs = 0.;
   std::string  b_prof;                      // "kEINASTO", ...
   Int_t        b_prof_code = -1;
   Double_t     b_par1 = 0., b_par2 = 0., b_par3 = 0.;
   Double_t     b_J = 0., b_J_Jcontinuum = 0.;
   Double_t     b_Mdelta = 0., b_Mtid = 0., b_Rtid = 0.;
   Double_t     b_Mequdens = 0., b_Requdens = 0., b_Dgal = 0.;

   // ---- derived (computed on the fly) -------------------------------------
   Double_t x_gc = 0., y_gc = 0., z_gc = 0.;      // galactocentric cartesian [kpc]
   Double_t theta_rs = 0., theta_delta = 0., theta_tid = 0.;   // angular radii [deg]
   Double_t Omega_tid = 0.;                       // solid angle within Rtid [sr]
   Double_t c_delta = 0., c_tid = 0.;             // concentrations
   Double_t f_strip = 0.;                         // Mtid/Mdelta
   Double_t rho_mean_tid = 0.;                    // 3 Mtid /(4pi Rtid^3)
   Double_t Lann = 0., Lann_delta = 0.;           // int rho^2 dV [Msol^2/kpc^3]
   Double_t J_tot = 0.;                           // Lann/d^2 [GeV^2 cm^-5]
   Double_t f_aperture = 0.;                      // J_file / J_tot
   Double_t alpha_int_b = 0.;                     // J half-angle used [deg]
   Bool_t   is_extended = kFALSE;                 // theta_tid > alpha_int
   Bool_t   tid_valid = kTRUE;                    // Rtid not pinned at its floor

   std::string *p_type = &b_type;
   std::string *p_prof = &b_prof;

   tree->Branch("name",          &b_name,          "name/L");
   tree->Branch("type",          &p_type);                       // std::string branch
   tree->Branch("long",          &b_long,          "long/D");
   tree->Branch("lat",           &b_lat,           "lat/D");
   tree->Branch("d",             &b_d,             "d/D");
   tree->Branch("z",             &b_z,             "z/D");
   tree->Branch("Rdelta",        &b_Rdelta,        "Rdelta/D");
   tree->Branch("rhos",          &b_rhos,          "rhos/D");
   tree->Branch("rs",            &b_rs,            "rs/D");
   tree->Branch("prof",          &p_prof);                       // textual in the file
   tree->Branch("prof_code",     &b_prof_code,     "prof_code/I");
   tree->Branch("par1",          &b_par1,          "par1/D");
   tree->Branch("par2",          &b_par2,          "par2/D");
   tree->Branch("par3",          &b_par3,          "par3/D");
   tree->Branch("J",             &b_J,             "J/D");
   tree->Branch("J_Jcontinuum",  &b_J_Jcontinuum,  "J_Jcontinuum/D");
   tree->Branch("Mdelta",        &b_Mdelta,        "Mdelta/D");
   tree->Branch("Mtid",          &b_Mtid,          "Mtid/D");
   tree->Branch("Rtid",          &b_Rtid,          "Rtid/D");
   tree->Branch("Mequdens",      &b_Mequdens,      "Mequdens/D");
   tree->Branch("Requdens",      &b_Requdens,      "Requdens/D");
   tree->Branch("Dgal",          &b_Dgal,          "Dgal/D");

   if (doExtras) {
      tree->Branch("x_gc",         &x_gc,         "x_gc/D");
      tree->Branch("y_gc",         &y_gc,         "y_gc/D");
      tree->Branch("z_gc",         &z_gc,         "z_gc/D");
      tree->Branch("theta_rs",     &theta_rs,     "theta_rs/D");
      tree->Branch("theta_delta",  &theta_delta,  "theta_delta/D");
      tree->Branch("theta_tid",    &theta_tid,    "theta_tid/D");
      tree->Branch("Omega_tid",    &Omega_tid,    "Omega_tid/D");
      tree->Branch("c_delta",      &c_delta,      "c_delta/D");
      tree->Branch("c_tid",        &c_tid,        "c_tid/D");
      tree->Branch("f_strip",      &f_strip,      "f_strip/D");
      tree->Branch("rho_mean_tid", &rho_mean_tid, "rho_mean_tid/D");
      tree->Branch("Lann",         &Lann,         "Lann/D");
      tree->Branch("Lann_delta",   &Lann_delta,   "Lann_delta/D");
      tree->Branch("J_tot",        &J_tot,        "J_tot/D");
      tree->Branch("f_aperture",   &f_aperture,   "f_aperture/D");
      tree->Branch("alpha_int",    &alpha_int_b,  "alpha_int/D");
      tree->Branch("is_extended",  &is_extended,  "is_extended/O");
      tree->Branch("tid_valid",    &tid_valid,    "tid_valid/O");
   }

   // "long" is a C++ keyword; give TTree::Draw a keyword-free handle as well.
   tree->SetAlias("glon", "long");
   tree->SetAlias("glat", "lat");

   ProfileCodes codes;
   std::vector<std::string> tok;
   std::string line;

   Long64_t lineNo = 0, nFilled = 0, nBad = 0, nTidFloor = 0, nNoLum = 0;
   const Long64_t kMaxReports = 20;
   int layout = 0; // 21 or 19, locked in on the first good data row
   double alphaIntHeader = -1.0;

   while (std::getline(in, line)) {
      ++lineNo;
      if (!line.empty() && line[line.size() - 1] == '\r') line.erase(line.size() - 1);
      if (IsCommentOrBlank(line)) {
         // The units header carries the J integration half-angle, e.g.
         //   ... | J(0.0299^{#circ}) [GeV^{2} cm^{-5}] ...
         const size_t jp = line.find("J(");
         if (jp != std::string::npos && alphaIntHeader < 0.) {
            const double v = std::strtod(line.c_str() + jp + 2, nullptr);
            if (v > 0.) {
               alphaIntHeader = v;
               std::cerr << "clumpy2root: J integration half-angle from header: "
                         << alphaIntHeader << " deg\n";
            }
         }
         continue;
      }

      Tokenize(line, tok);
      const int n = static_cast<int>(tok.size());

      if (layout == 0) {
         if (n == 21 || n == 19) {
            layout = n;
            std::cerr << "clumpy2root: detected " << layout << "-column layout"
                      << (layout == 21 ? " (with Name/Type columns)"
                                       : " (no Name/Type columns)") << "\n";
         } else {
            if (++nBad <= kMaxReports)
               std::cerr << "clumpy2root: WARNING: line " << lineNo << ": expected 21 or 19"
                         << " columns, got " << n << " -- skipped\n";
            continue;
         }
      }

      if (n != layout) {
         if (++nBad <= kMaxReports)
            std::cerr << "clumpy2root: WARNING: line " << lineNo << ": expected " << layout
                      << " columns, got " << n << " -- skipped\n";
         continue;
      }

      // Field offset: 21-column files carry Name+Type in front of 'long'.
      const int o = (layout == 21) ? 2 : 0;

      bool ok = true;
      if (layout == 21) {
         if (!ToLong64(tok[0], b_name)) ok = false;
         b_type = tok[1];
      } else {
         b_name = nFilled + 1;   // synthesise a 1-based index
         b_type.clear();
      }

      // The profile column is textual; everything else on the row is numeric.
      const int profIdx = o + 7;
      double *dst[] = {&b_long, &b_lat, &b_d, &b_z, &b_Rdelta, &b_rhos, &b_rs,
                       nullptr, // placeholder for 'prof'
                       &b_par1, &b_par2, &b_par3, &b_J, &b_J_Jcontinuum,
                       &b_Mdelta, &b_Mtid, &b_Rtid, &b_Mequdens, &b_Requdens, &b_Dgal};
      int badCol = -1;
      for (int k = 0; ok && k < 19; ++k) {
         if (!dst[k]) continue;
         if (!ToDouble(tok[o + k], *dst[k])) { ok = false; badCol = o + k; }
      }

      if (!ok) {
         if (++nBad <= kMaxReports)
            std::cerr << "clumpy2root: WARNING: line " << lineNo << ": unparsable field "
                      << (badCol >= 0 ? tok[badCol] : tok[0]) << " in column "
                      << (badCol >= 0 ? badCol : 0) + 1 << " -- skipped\n";
         continue;
      }

      b_prof = tok[profIdx];
      if (b_prof.empty() || std::isdigit(static_cast<unsigned char>(b_prof[0]))) {
         // Some CLUMPY builds write the enum as a bare integer instead of a name.
         Long64_t asInt = -1;
         if (ToLong64(b_prof, asInt)) {
            b_prof_code = static_cast<Int_t>(asInt);
         } else {
            if (++nBad <= kMaxReports)
               std::cerr << "clumpy2root: WARNING: line " << lineNo
                         << ": bad profile field '" << b_prof << "' -- skipped\n";
            continue;
         }
      } else {
         b_prof_code = codes.Code(b_prof);
      }
      codes.Note(b_prof, b_prof_code);

      if (doExtras) {
         if (alphaInt < 0.) alphaInt = alphaIntHeader;   // header value, if any
         const double deg = TMath::Pi() / 180.0;
         const double cb = std::cos(b_lat * deg), sb = std::sin(b_lat * deg);
         const double cl = std::cos(b_long * deg), sl = std::sin(b_long * deg);

         // Galactocentric cartesian: GC at the origin, Sun at (-rSun, 0, 0),
         // so that l = b = 0 seen from the Sun points along +x.
         x_gc = b_d * cb * cl - rSun;
         y_gc = b_d * cb * sl;
         z_gc = b_d * sb;

         theta_rs    = (b_d > 0.) ? std::atan2(b_rs,     b_d) / deg : 0.;
         theta_delta = (b_d > 0.) ? std::atan2(b_Rdelta, b_d) / deg : 0.;
         theta_tid   = (b_d > 0.) ? std::atan2(b_Rtid,   b_d) / deg : 0.;
         Omega_tid   = 2.0 * TMath::Pi() * (1.0 - std::cos(theta_tid * deg));

         c_delta = (b_rs > 0.) ? b_Rdelta / b_rs : 0.;
         c_tid   = (b_rs > 0.) ? b_Rtid   / b_rs : 0.;
         f_strip = (b_Mdelta > 0.) ? b_Mtid / b_Mdelta : 0.;
         rho_mean_tid = (b_Rtid > 0.)
                        ? 3.0 * b_Mtid / (4.0 * TMath::Pi() * b_Rtid * b_Rtid * b_Rtid)
                        : 0.;

         // 14% of the rows in these catalogues carry Rtid pinned exactly at
         // rs*1e-3 -- the solver's lower bracket, not a physical tidal radius.
         // Flag them rather than silently propagating the value.
         tid_valid = !(c_tid > 0. && std::fabs(c_tid / 1e-3 - 1.0) < 0.01);
         if (!tid_valid) ++nTidFloor;

         const double Rlum = tid_valid ? b_Rtid : b_Rdelta;
         switch (b_prof_code) {
         case 1: case 2:                                   // kEINASTO(_N)
            Lann       = EinastoLuminosity(b_rhos, b_rs, b_par1, Rlum);
            Lann_delta = EinastoLuminosity(b_rhos, b_rs, b_par1, b_Rdelta);
            break;
         case 0: {                                         // kZHAO (a, b, g)
            const double rhos = b_rhos, rs = b_rs;
            const double A = b_par1, B = b_par2, G = b_par3;
            auto rho = [=](double r) { return ZhaoRho(r, rhos, rs, A, B, G); };
            Lann       = NumericLuminosity(rho, Rlum);
            Lann_delta = NumericLuminosity(rho, b_Rdelta);
            break;
         }
         case 3: {                                         // kNFW = Zhao(1,3,1)
            const double rhos = b_rhos, rs = b_rs;
            auto rho = [=](double r) { return ZhaoRho(r, rhos, rs, 1., 3., 1.); };
            Lann       = NumericLuminosity(rho, Rlum);
            Lann_delta = NumericLuminosity(rho, b_Rdelta);
            break;
         }
         default:                                          // unsupported profile
            Lann = Lann_delta = 0.;
            if (nNoLum++ == 0)
               std::cerr << "clumpy2root: note: no luminosity model for profile '"
                         << b_prof << "'; Lann/J_tot left at 0\n";
            break;
         }

         // CLUMPY's J is L/d^2 (integrated over solid angle, no 1/4pi).
         J_tot      = (b_d > 0.) ? Lann / (b_d * b_d) * kLtoJ : 0.;
         f_aperture = (J_tot > 0.) ? b_J / J_tot : 0.;

         alpha_int_b = alphaInt;
         const double thEff = tid_valid ? theta_tid : theta_delta;
         is_extended = (alpha_int_b > 0.) && (thEff > alpha_int_b);
      }

      tree->Fill();
      ++nFilled;
      if (maxLines > 0 && nFilled >= maxLines) break;
   }

   if (nBad > kMaxReports)
      std::cerr << "clumpy2root: WARNING: " << (nBad - kMaxReports)
                << " further malformed lines not reported individually\n";


   // ------------------------------------------------------------------------
   // Self-describing metadata, so that this file can still be understood years
   // from now with no access to the source, the CLUMPY run, or this session.
   //   * every TBranch gets a human-readable title  -> visible in clumps->Print()
   //   * tree->GetUserInfo() carries the machine-readable dictionary
   //   * a "README" TObjString carries the whole thing as printable text
   //   * a "provenance" TList records where the file came from
   // ------------------------------------------------------------------------
   {
      TList *info = tree->GetUserInfo();
      info->SetName("branch_dictionary");

      int nUndocumented = 0;
      TIter nextB(tree->GetListOfBranches());
      while (TObject *o = nextB()) {
         TBranch *br = static_cast<TBranch *>(o);
         const BranchDoc *doc = FindDoc(br->GetName());
         if (!doc) {
            ++nUndocumented;
            std::cerr << "clumpy2root: WARNING: branch '" << br->GetName()
                      << "' has no dictionary entry\n";
            continue;
         }
         // Title shows up in TTree::Print(); the leaves already exist, so this is
         // purely cosmetic and cannot affect how the branch is read back.
         br->SetTitle(Form("%s [%s]", doc->brief, doc->unit));
         info->Add(new TNamed(doc->name, Form("[%s] %s -- %s", doc->unit,
                                              doc->brief, doc->desc)));
      }
      if (nUndocumented)
         std::cerr << "clumpy2root: WARNING: " << nUndocumented
                   << " branch(es) undocumented -- update kBranchDoc\n";

      // ---- provenance -----------------------------------------------------
      TList *prov = new TList();
      prov->SetName("provenance");
      auto add = [&](const char *k, const char *v) { prov->Add(new TNamed(k, v)); };

      TDatime now;
      char *inAbs = gSystem->ExpandPathName(inPath.c_str());
      Long_t fid, fsize, fflags, fmtime;
      gSystem->GetPathInfo(inPath.c_str(), &fid, &fsize, &fflags, &fmtime);

      std::string cmdline;
      for (int i = 0; i < argc; ++i) cmdline += std::string(i ? " " : "") + argv[i];

      add("producer",        Form("clumpy2root %s", kVersion));
      add("root_version",    gROOT->GetVersion());
      add("converted_utc",   now.AsSQLString());
      add("host",            gSystem->HostName());
      add("user",            gSystem->Getenv("USER") ? gSystem->Getenv("USER") : "?");
      add("command_line",    cmdline.c_str());
      add("input_file",      inAbs ? inAbs : inPath.c_str());
      add("input_bytes",     Form("%ld", fsize));
      add("input_mtime",     TDatime(fmtime).AsSQLString());
      add("input_layout",    Form("%d columns", layout));
      add("lines_read",      Form("%lld", lineNo));
      add("entries_filled",  Form("%lld", nFilled));
      add("lines_rejected",  Form("%lld", nBad));
      add("derived_columns", doExtras ? "yes" : "no");
      add("R_sun_kpc",       Form("%.6g", rSun));
      add("alpha_int_deg",   Form("%.6g", alphaInt));
      add("tid_floor_rows",  Form("%lld", nTidFloor));
      if (inAbs) delete[] inAbs;

      // ---- printable README ------------------------------------------------
      std::string rd;
      rd += "================================================================\n";
      rd += " CLUMPY clump catalogue converted to a ROOT TTree\n";
      rd += " produced by clumpy2root " + std::string(kVersion) + "\n";
      rd += "================================================================\n\n";
      rd += "TREE\n";
      rd += Form("  '%s': one entry per dark-matter clump (subhalo), %lld entries.\n\n",
                 treeName.c_str(), nFilled);
      rd += "PROVENANCE\n";
      { TIter it(prov);
        while (TObject *o = it()) {
           TNamed *nm = static_cast<TNamed *>(o);
           rd += Form("  %-16s %s\n", nm->GetName(), nm->GetTitle());
        } }
      rd += "\nCONVENTIONS ESTABLISHED FROM THE DATA (not assumed)\n";
      rd += Wrap("Rdelta is R_200c: the mean density inside it is 2.55e4 Msol/kpc^3 "
                 "with 0.4% scatter, i.e. 200 x rho_crit for h = 0.678.", 74, "  * ", "    ");
      rd += Wrap("rs = r_-2 and rhos = rho_-2, and par1 is the Einasto shape index "
                 "alpha: re-integrating the profile reproduces the file's own Mdelta, "
                 "Mtid and Mequdens to 0.4%, its printing precision.", 74, "  * ", "    ");
      rd += Wrap("CLUMPY's J is normalised as L/d^2, with NO factor 1/(4*pi): median "
                 "J/J_tot = 1.0005 +- 0.005 over 248502 clumps that fit inside the "
                 "integration aperture.", 74, "  * ", "    ");
      rd += Wrap("R_sun = 8.0 kpc, recovered from (long,lat,d) against the Dgal column "
                 "to 1 part in 1e5.", 74, "  * ", "    ");
      rd += "\nKNOWN DATA-QUALITY ISSUE\n";
      rd += Wrap("A fraction of rows carry Rtid pinned at exactly rs*1e-3 -- the lower "
                 "bracket of CLUMPY's tidal-radius solver, not a converged physical "
                 "radius; their Mtid is ~1e-7 x Mdelta. Those rows have tid_valid = "
                 "false. ALWAYS cut on tid_valid before using Rtid, Mtid, c_tid, "
                 "f_strip, rho_mean_tid, theta_tid or Omega_tid.", 74, "  ", "  ");
      rd += Form("  In this file: %lld of %lld rows (%.2f%%) are affected.\n",
                 nTidFloor, nFilled, nFilled ? 100.0 * nTidFloor / nFilled : 0.0);
      rd += "\nBRANCH DICTIONARY\n";
      { TIter it(info);
        while (TObject *o = it()) {
           TNamed *nm = static_cast<TNamed *>(o);
           rd += Form("\n  %s\n", nm->GetName());
           rd += Wrap(nm->GetTitle(), 72, "      ", "      ");
        } }
      rd += "\nOTHER OBJECTS IN THIS FILE\n";
      rd += "  prof_codes   TList, profile-name -> integer code map for prof_code\n";
      rd += "  provenance   TList, the table printed above\n";
      rd += "  README       this text\n";
      rd += "\nQUICK START\n";
      rd += "  root -l <this file>\n";
      rd += "  root [1] README->Print()                    // reprint this text\n";
      rd += Form("  root [2] %s->Print()                    // branches + titles\n",
                 treeName.c_str());
      rd += Form("  root [3] %s->Scan(\"long:lat:d\",\"\",\"\",10)\n", treeName.c_str());
      rd += Form("  root [4] %s->Draw(\"J_tot\",\"tid_valid\")\n", treeName.c_str());

      // TMacro holds text line-by-line and prints it verbatim, so a bare
      // "README->Print()" at the ROOT prompt renders this readably.
      TMacro readme("README", "clumpy2root provenance and branch dictionary");
      { std::istringstream ls(rd); std::string one;
        while (std::getline(ls, one)) readme.AddLine(one.c_str()); }

      fout.cd();
      readme.Write("README", TObject::kSingleKey);
      prov->Write("provenance", TObject::kSingleKey);
      prov->Delete();
      delete prov;
   }

   fout.cd();
   tree->Write("", TObject::kOverwrite);

   // Persist the profile-name <-> code dictionary next to the tree.
   TList profMap;
   profMap.SetName("prof_codes");
   for (std::map<std::string, int>::const_iterator it = codes.Seen().begin();
        it != codes.Seen().end(); ++it) {
      profMap.Add(new TNamed(it->first.c_str(), Form("%d", it->second)));
   }
   profMap.Write("prof_codes", TObject::kSingleKey);
   profMap.Delete();

   fout.Close();

   std::cout << "clumpy2root: read   " << lineNo  << " lines from " << inPath  << "\n"
             << "clumpy2root: filled " << nFilled << " entries into TTree '" << treeName
             << "' in " << outPath << "\n"
             << "clumpy2root: skipped " << nBad << " malformed line(s)\n";
   if (doExtras) {
      std::cout << "clumpy2root: derived columns on (R_sun = " << rSun << " kpc, "
                << "alpha_int = " << alphaInt << " deg)\n";
      if (nTidFloor)
         std::cout << "clumpy2root: " << nTidFloor << " clump(s) ("
                   << (100.0 * nTidFloor / (nFilled ? nFilled : 1))
                   << "%) have Rtid pinned at rs*1e-3 -> tid_valid = false\n";
      if (nNoLum)
         std::cout << "clumpy2root: " << nNoLum
                   << " clump(s) with an unsupported profile -> Lann = 0\n";
   }

   return (nFilled > 0) ? 0 : 3;
}
