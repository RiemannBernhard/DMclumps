#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <sstream>

#include <cstdio>
#include <string>
#include <math.h>

#include <algorithm>
#include <iterator>



#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"

#include "/home/blessed/Documents/_HistogramScripts/DAMPE/draw_plots/include/MyUtility_func.h"

char *PLOT_PATH = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/plots";

char *histoMain_Xaxis_title = (char*)"redshift, z";
char *histoMain_Yaxis_title = (char*)"H(z) [km sec^{-1} Mpc^{-1}]";

Float_t Latex1ScaleFactor(0.032);
Float_t Latex2ScaleFactor(0.025);
const Int_t rebinFactor = 1;
const Int_t YesLogScale = 0;
const Int_t XAxis_SetNdivisions   = 10+100*12;
const Int_t YAxis_SetNdivisions   = 10+100*6;

void  draw_fits(const char *file, 
                bool Xaxis_LOGSCALE, 
                bool Yaxis_LOGSCALE, 
                Float_t xMin_l,Float_t xMax_l,
                Float_t yMin_l,Float_t yMax_l,
                char *fileNameSuff);



int FLAG=2;

void convert_clumpy_draw_annih_rs_to_root(){
  
    //draw_fits("/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/annihil_gal2D_LOS180_0_FOVdiameter360.0deg_nside1024.drawn",
    draw_fits("/opt/CLUMPY/output/annihil_rs01_gamma052D_FOVdiameter8.0deg_nside512.drawn",
              false, false,
              0, 2.5,
              0, 300,
              //(char*)"clumpy_annihil_gal2D_LOS180_0_FOVdiameter360.0deg_nside1024");
              (char*)"clumpy_annihil_rs01_seed666_gamma052D_FOV8deg_nside512");

}

void draw_fits(const char *file, 
               bool Xaxis_LOGSCALE, 
               bool Yaxis_LOGSCALE, 
               Float_t xMin_l,Float_t xMax_l,
               Float_t yMin_l,Float_t yMax_l,
               char *fileNameSuff)
{


  std::ifstream ifstr(file, std::ifstream::in);
  if (ifstr.fail()) {
    std::cerr << "\n Could not open " << file << "\n\n";
    return;
  }
  else {
     cout << "File is successfully opened.\n" << endl;
     //cout << "ifstream.rdbuf  = " << ifstr.rdbuf()   << endl;
     cout << "ifstream.get    = " << ifstr.get()   << endl;
  }
  
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_gal2D_LOS180_0_FOVdiameter360.0deg_nside1024.root";
    const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_ZHAO_rs01_seed666_gamma052D_FOV8deg_nside512.root";

    TFile *file_output = new TFile (output_fname,"RECREATE");

   float mGlong;
   float mGlat;
   float mDist;
   float mGlong_wrtHaloCenter;
   float mGlat_wrtHaloCenter;

   float mZ;
   float mRdelta;
   float mRhos;
   float mRhos2;
   float mRs;

   float mParam1;
   float mParam2;
   float mParam3;
   float mJfactor;
   float mJJcont;

   float mMdelta;
   float mMtidal;
   float mRtidal;
   float mMequdens;
   float mRequdens;
   float mDgal;

  auto mTree = new TTree("CLUMPY_output", "CLUMPY_output");
  mTree->Branch("GalacLongitude_deg",       &mGlong,          "mGlong/F");
  mTree->Branch("GalacLatitude_deg",        &mGlat,           "mGlat/F");

  mTree->Branch("GalacLongitude_wrtHaloCenter_deg",       &mGlong_wrtHaloCenter,          "mGlong_wrtHaloCenter/F");
  mTree->Branch("GalacLatitude_wrtHaloCenter_deg",        &mGlat_wrtHaloCenter,           "mGlat_wrtHaloCenter/F");


  mTree->Branch("Distance_kpc",             &mDist,           "mDist/F");
  mTree->Branch("z",                        &mZ,              "mZ/F");
  mTree->Branch("Rdelta_kpc",               &mRdelta,         "mRdelta/F");
  mTree->Branch("DMprof_rhos_Msolkpc3",     &mRhos,           "mRhos/F");
  mTree->Branch("DMprof_rhos2_Msolkpc3_2",  &mRhos2,          "mRhos2/F");

  mTree->Branch("DMprof_ScaleRadius_rs_kpc",  &mRs,             "mRs/F");
  mTree->Branch("DMprof_p1",                  &mParam1,         "mParam1/F");
  mTree->Branch("DMprof_p2",                  &mParam2,         "mParam2/F");
  mTree->Branch("DMprof_p3",                  &mParam3,         "mParam3/F");
  mTree->Branch("DM_Jfactor_GeV2cm5",         &mJfactor,        "mJfactor/F");
  mTree->Branch("DM_JJcontinuum_GeV2cm5",     &mJJcont,         "mJJcont/F");
  mTree->Branch("DM_Mdelta_Msol",             &mMdelta,         "mMdelta/F");
  mTree->Branch("DM_Mtidal_Msol",             &mMtidal,         "mMtidal/F");
  mTree->Branch("DM_Rtidal_kpc",              &mRtidal,         "mRtidal/F");
  mTree->Branch("DM_Mequdens_Msol",           &mMequdens,       "mMequdens/F");
  mTree->Branch("DM_Requdens_kpc",            &mRequdens,       "mRequdens/F");
  mTree->Branch("DM_Dgal_kpc",                &mDgal,           "mDgal/F");

  TProfile *prf_rhos_Mdelta     = new TProfile ("rhos_Mdelta",     ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT M_{#Delta} #GT [M_{sun}]",             10000, 0, 1.5e9, 0, 2000);
  TProfile *prf_rhos_Mtid       = new TProfile ("rhos_Mtid",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT M_{tidal} #GT [M_{sun}]",              10000, 0, 1.5e9, 0, 2000);
  TProfile *prf_rhos_Rtid       = new TProfile ("rhos_Rtid",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT R_{tidal} #GT [kpc]",                  10000, 0, 1.5e9, 0, 1);
  TProfile *prf_rhos_Mequdens   = new TProfile ("rhos_Mequdens",   ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT M_{equdens} #GT [M_{sun}]",            10000, 0, 1.5e9, 0, 2000);
  TProfile *prf_rhos_Requdens   = new TProfile ("rhos_Requdens",   ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT R_{equdens} #GT [kpc]",                10000, 0, 1.5e9, 0, 1);
  TProfile *prf_rhos_J          = new TProfile ("rhos_J",          ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        10000, 0, 1.5e9, 0, 1e8);
  TProfile *prf_rhos_JJcont     = new TProfile ("rhos_JJcont",     ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 10000, 0, 1.5e9, 0, 100);
  TProfile *prf_rhos_dist       = new TProfile ("rhos_dist",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT dist #GT [kpc]",                       10000, 0, 1.5e9, 0, 120);

  TProfile *prf_rs_Mdelta     = new TProfile ("rs_Mdelta",     ";scaled radius r_{s} [kpc];#LT M_{#Delta} #GT [M_{sun}]",                   40000, 0, 4e-3, 0, 2000);
  TProfile *prf_rs_Mtid       = new TProfile ("rs_Mtid",       ";scaled radius r_{s} [kpc];#LT M_{tidal} #GT [M_{sun}]",                    40000, 0, 4e-3, 0, 2000);
  TProfile *prf_rs_Rtid       = new TProfile ("rs_Rtid",       ";scaled radius r_{s} [kpc];#LT R_{tidal} #GT [kpc]",                        40000, 0, 4e-3, 0, 1);
  TProfile *prf_rs_Mequdens   = new TProfile ("rs_Mequdens",   ";scaled radius r_{s} [kpc];#LT M_{equdens} #GT [M_{sun}]",                  40000, 0, 4e-3, 0, 2000);
  TProfile *prf_rs_Requdens   = new TProfile ("rs_Requdens",   ";scaled radius r_{s} [kpc];#LT R_{equdens} #GT [kpc]",                      40000, 0, 4e-3, 0, 1);
  TProfile *prf_rs_J          = new TProfile ("rs_J",          ";scaled radius r_{s} [kpc];#LT J-factor #GT [GeV^{2}cm^{-5}]",              40000, 0, 4e-3, 0, 1e8);
  TProfile *prf_rs_JJcont     = new TProfile ("rs_JJcont",     ";scaled radius r_{s} [kpc];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]",       40000, 0, 4e-3, 0, 100);
  TProfile *prf_rs_rhos       = new TProfile ("rs_rhos",       ";scaled radius r_{s} [kpc];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",             40000, 0, 4e-3, 0, 2e9);
  TProfile *prf_rs_rhos2      = new TProfile ("rs_rhos2",      ";scaled radius r_{s} [kpc];#LT #rho_{s}^{2} #GT [(M_{sun}/kpc^{3})^{2}]",   40000, 0, 4e-3, 0, 1e18);
  TProfile *prf_rs_dist       = new TProfile ("rs_dist",       ";scaled radius r_{s} [kpc];#LT dist #GT [kpc]",                             40000, 0, 4e-3, 0, 120);

  TProfile *prf_Glong_Mdelta     = new TProfile ("Glong_Mdelta",    ";Galactic longitude [deg];#LT M_{#Delta} #GT [M_{sun}]",             10000, -200, 200, 0, 2000);
  TProfile *prf_Glong_Mtid       = new TProfile ("Glong_Mtid",      ";Galactic longitude [deg];#LT M_{tidal} #GT [M_{sun}]",              10000, -200, 200, 0, 2000);
  TProfile *prf_Glong_Rtid       = new TProfile ("Glong_Rtid",      ";Galactic longitude [deg];#LT R_{tidal} #GT [kpc]",                  10000, -200, 200, 0, 1);
  TProfile *prf_Glong_Mequdens   = new TProfile ("Glong_Mequdens",  ";Galactic longitude [deg];#LT M_{equdens} #GT [M_{sun}]",            10000, -200, 200, 0, 2000);
  TProfile *prf_Glong_Requdens   = new TProfile ("Glong_Requdens",  ";Galactic longitude [deg];#LT R_{equdens} #GT [kpc]",                10000, -200, 200, 0, 1);
  TProfile *prf_Glong_J          = new TProfile ("Glong_J",         ";Galactic longitude [deg];#LT J-factor #GT [GeV^{2}cm^{-5}]",        10000, -200, 200, 0, 1e8);
  TProfile *prf_Glong_JJcont     = new TProfile ("Glong_JJcont",    ";Galactic longitude [deg];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 10000, -200, 200, 0, 100);

  TProfile *prf_Glat_Mdelta     = new TProfile ("Glat_Mdelta",    ";Galactic latitude [deg];#LT M_{#Delta} #GT [M_{sun}]",             10000, -100, 100, 0, 2000);
  TProfile *prf_Glat_Mtid       = new TProfile ("Glat_Mtid",      ";Galactic latitude [deg];#LT M_{tidal} #GT [M_{sun}]",              10000, -100, 100, 0, 2000);
  TProfile *prf_Glat_Rtid       = new TProfile ("Glat_Rtid",      ";Galactic latitude [deg];#LT R_{tidal} #GT [kpc]",                  10000, -100, 100, 0, 1);
  TProfile *prf_Glat_Mequdens   = new TProfile ("Glat_Mequdens",  ";Galactic latitude [deg];#LT M_{equdens} #GT [M_{sun}]",            10000, -100, 100, 0, 2000);
  TProfile *prf_Glat_Requdens   = new TProfile ("Glat_Requdens",  ";Galactic latitude [deg];#LT R_{equdens} #GT [kpc]",                10000, -100, 100, 0, 1);
  TProfile *prf_Glat_J          = new TProfile ("Glat_J",         ";Galactic latitude [deg];#LT J-factor #GT [GeV^{2}cm^{-5}]",        10000, -100, 100, 0, 1e8);
  TProfile *prf_Glat_JJcont     = new TProfile ("Glat_JJcont",    ";Galactic latitude [deg];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 10000, -100, 100, 0, 100);

  TProfile *prf_dist_Mdelta     = new TProfile ("dist_Mdelta",    ";Distance [kpc];#LT M_{#Delta} #GT [M_{sun}]",             20000, 0, 120, 0, 2000);
  TProfile *prf_dist_Mtid       = new TProfile ("dist_Mtid",      ";Distance [kpc];#LT M_{tidal} #GT [M_{sun}]",              20000, 0, 120, 0, 2000);
  TProfile *prf_dist_Rtid       = new TProfile ("dist_Rtid",      ";Distance [kpc];#LT R_{tidal} #GT [kpc]",                  20000, 0, 120, 0, 1);
  TProfile *prf_dist_Mequdens   = new TProfile ("dist_Mequdens",  ";Distance [kpc];#LT M_{equdens} #GT [M_{sun}]",            20000, 0, 120, 0, 2000);
  TProfile *prf_dist_Requdens   = new TProfile ("dist_Requdens",  ";Distance [kpc];#LT R_{equdens} #GT [kpc]",                20000, 0, 120, 0, 1);
  TProfile *prf_dist_J          = new TProfile ("dist_J",         ";Distance [kpc];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 120, 0, 1e8);
  TProfile *prf_dist_JJcont     = new TProfile ("dist_JJcont",    ";Distance [kpc];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 120, 0, 100);

  TProfile *prf_Mdelta_rhos       = new TProfile ("Mdelta_rhos",       ";M_{#Delta} [M_{sun}];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       20000, 0, 200, 0, 2e9);
  TProfile *prf_Mdelta_rs         = new TProfile ("Mdelta_rs",         ";M_{#Delta} [M_{sun}];#LT r_{s} #GT [kpc]",                      20000, 0, 200, 0, 1e-3);
  TProfile *prf_Mdelta_Mtid       = new TProfile ("Mdelta_Mtid",       ";M_{#Delta} [M_{sun}];#LT M_{tidal} #GT [M_{sun}]",              20000, 0, 200, 0, 2000);
  TProfile *prf_Mdelta_Rtid       = new TProfile ("Mdelta_Rtid",       ";M_{#Delta} [M_{sun}];#LT R_{tidal} #GT [kpc]",                  20000, 0, 200, 0, 1);
  TProfile *prf_Mdelta_Mequdens   = new TProfile ("Mdelta_Mequdens",   ";M_{#Delta} [M_{sun}];#LT M_{equdens} #GT [M_{sun}]",            20000, 0, 200, 0, 2000);
  TProfile *prf_Mdelta_Requdens   = new TProfile ("Mdelta_Requdens",   ";M_{#Delta} [M_{sun}];#LT R_{equdens} #GT [kpc]",                20000, 0, 200, 0, 1);
  TProfile *prf_Mdelta_J          = new TProfile ("Mdelta_J",          ";M_{#Delta} [M_{sun}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 200, 0, 1e8);
  TProfile *prf_Mdelta_JJcont     = new TProfile ("Mdelta_JJcont",     ";M_{#Delta} [M_{sun}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 200, 0, 100);

  TProfile *prf_Mtid_rhos       = new TProfile ("Mtid_rhos",      ";M_{tidal} [M_{sun}];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       20000, 0, 800, 0, 2e9);
  TProfile *prf_Mtid_rs         = new TProfile ("Mtid_rs",        ";M_{tidal} [M_{sun}];#LT r_{s} #GT [kpc]",                      20000, 0, 800, 0, 1e-3);
  TProfile *prf_Mtid_Mdelta     = new TProfile ("Mtid_Mdelta",    ";M_{tidal} [M_{sun}];#LT M_{#Delta} #GT [M_{sun}]",             20000, 0, 800, 0, 2000);
  TProfile *prf_Mtid_Rtid       = new TProfile ("Mtid_Rtid",      ";M_{tidal} [M_{sun}];#LT R_{tidal} #GT [kpc]",                  20000, 0, 800, 0, 1);
  TProfile *prf_Mtid_Mequdens   = new TProfile ("Mtid_Mequdens",  ";M_{tidal} [M_{sun}];#LT M_{equdens} #GT [M_{sun}]",            20000, 0, 800, 0, 2000);
  TProfile *prf_Mtid_Requdens   = new TProfile ("Mtid_Requdens",  ";M_{tidal} [M_{sun}];#LT R_{equdens} #GT [kpc]",                20000, 0, 800, 0, 1);
  TProfile *prf_Mtid_J          = new TProfile ("Mtid_J",         ";M_{tidal} [M_{sun}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 800, 0, 1e8);
  TProfile *prf_Mtid_JJcont     = new TProfile ("Mtid_JJcont",    ";M_{tidal} [M_{sun}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 800, 0, 100);

  TProfile *prf_Rtid_rhos       = new TProfile ("Rtid_rhos",      ";R_{tidal} [kpc];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",        40000, 0.002, 0.05, 0, 2e9);
  TProfile *prf_Rtid_rs         = new TProfile ("Rtid_rs",        ";R_{tidal} [kpc];#LT r_{s} #GT [kpc]",                       40000, 0.002, 0.05, 0, 1e-3);
  TProfile *prf_Rtid_Mdelta     = new TProfile ("Rtid_Mdelta",    ";R_{tidal} [kpc];#LT M_{#Delta} #GT [M_{sun}]",              40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Rtid_Mtid       = new TProfile ("Rtid_Mtid",      ";R_{tidal} [kpc];#LT M_{tidal} #GT [M_{sun}]",               40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Rtid_Mequdens   = new TProfile ("Rtid_Mequdens",  ";R_{tidal} [kpc];#LT M_{equdens} #GT [M_{sun}]",             40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Rtid_Requdens   = new TProfile ("Rtid_Requdens",  ";R_{tidal} [kpc];#LT R_{equdens} #GT [kpc]",                 40000, 0.002, 0.05, 0, 1);
  TProfile *prf_Rtid_J          = new TProfile ("Rtid_J",         ";R_{tidal} [kpc];#LT J-factor #GT [GeV^{2}cm^{-5}]",         40000, 0.002, 0.05, 0, 1e8);
  TProfile *prf_Rtid_JJcont     = new TProfile ("Rtid_JJcont",    ";R_{tidal} [kpc];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]",  40000, 0.002, 0.05, 0, 100);

  TProfile *prf_Mequdens_rhos       = new TProfile ("Mequdens_rhos",      ";M_{equdens} [M_{sun}];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       40000, 0, 500, 0, 2e9);
  TProfile *prf_Mequdens_rs         = new TProfile ("Mequdens_rs",        ";M_{equdens} [M_{sun}];#LT r_{s} #GT [kpc]",                      40000, 0, 500, 0, 1e-3);
  TProfile *prf_Mequdens_Mdelta     = new TProfile ("Mequdens_Mdelta",    ";M_{equdens} [M_{sun}];#LT M_{#Delta} #GT [M_{sun}]",             40000, 0, 500, 0, 2000);
  TProfile *prf_Mequdens_Requdens   = new TProfile ("Mequdens_Requdens",  ";M_{equdens} [M_{sun}];#LT R_{equdens} #GT [kpc]",                40000, 0, 500, 0, 1);
  TProfile *prf_Mequdens_Mtid       = new TProfile ("Mequdens_Mtid",      ";M_{equdens} [M_{sun}];#LT M_{tidal} #GT [M_{sun}]",              40000, 0, 500, 0, 2000);
  TProfile *prf_Mequdens_Rtid       = new TProfile ("Mequdens_Rtid",      ";M_{equdens} [M_{sun}];#LT R_{tidal} #GT [kpc]",                  40000, 0, 500, 0, 2000);
  TProfile *prf_Mequdens_J          = new TProfile ("Mequdens_J",         ";M_{equdens} [M_{sun}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        40000, 0, 500, 0, 1e8);
  TProfile *prf_Mequdens_JJcont     = new TProfile ("Mequdens_JJcont",    ";M_{equdens} [M_{sun}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 40000, 0, 500, 0, 100);

  TProfile *prf_Requdens_rhos       = new TProfile ("Requdens_rhos",      ";R_{equdens} [kpc];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       40000, 0.002, 0.05, 0, 2e9);
  TProfile *prf_Requdens_rs         = new TProfile ("Requdens_rs",        ";R_{equdens} [kpc];#LT r_{s} #GT [kpc]",                      40000, 0.002, 0.05, 0, 1e-3);
  TProfile *prf_Requdens_Mdelta     = new TProfile ("Requdens_Mdelta",    ";R_{equdens} [kpc];#LT M_{#Delta} #GT [M_{sun}]",             40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Requdens_Mequdens   = new TProfile ("Requdens_Mequdens",  ";R_{equdens} [kpc];#LT M_{equdens} #GT [M_{sun}]",            40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Requdens_Rtid       = new TProfile ("Requdens_Rtid",      ";R_{equdens} [kpc];#LT R_{tidal} #GT [kpc]",                  40000, 0.002, 0.05, 0, 1);
  TProfile *prf_Requdens_Mtid       = new TProfile ("Requdens_Mtid",      ";R_{equdens} [kpc];#LT M_{tidal} #GT [M_{sun}]",              40000, 0.002, 0.05, 0, 1000);
  TProfile *prf_Requdens_J          = new TProfile ("Requdens_J",         ";R_{equdens} [kpc];#LT J-factor #GT [GeV^{2}cm^{-5}]",        40000, 0.002, 0.05, 0, 1e8);
  TProfile *prf_Requdens_JJcont     = new TProfile ("Requdens_JJcont",    ";R_{equdens} [kpc];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 40000, 0.002, 0.05, 0, 100);





  TH2F *h2d_rhos_Mdelta     = new TH2F ("h2d_rhos_Mdelta",     ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];M_{#Delta} [M_{sun}]",             4000, 0, 1.0e9, 1000, 0, 1000);
  TH2F *h2d_rhos_Mtid       = new TH2F ("h2d_rhos_Mtid",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];M_{tidal} [M_{sun}]",              4000, 0, 1.0e9, 1000, 0, 800);
  TH2F *h2d_rhos_Rtid       = new TH2F ("h2d_rhos_Rtid",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];R_{tidal} [kpc]",                  4000, 0, 1.0e9, 1000, 1e-8, 5e-2);
  TH2F *h2d_rhos_Mequdens   = new TH2F ("h2d_rhos_Mequdens",   ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];M_{equdens} [M_{sun}]",            4000, 0, 1.0e9, 1000, 0, 1000);
  TH2F *h2d_rhos_Requdens   = new TH2F ("h2d_rhos_Requdens",   ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];R_{equdens} [kpc]",                4000, 0, 1.0e9, 1000, 1e-8, 5e-2);
  TH2F *h2d_rhos_J          = new TH2F ("h2d_rhos_J",          ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];J-factor [GeV^{2}cm^{-5}]",        4000, 0, 1.0e9, 1000, 0, 1e8);
  TH2F *h2d_rhos_JJcont     = new TH2F ("h2d_rhos_JJcont",     ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];J/J_{continuum} [GeV^{2}cm^{-5}]", 4000, 0, 1.0e9, 1000, 0, 40);

  TH2F *h2d_rs_Mdelta     = new TH2F ("h2d_rs_Mdelta",     ";scaled radius r_{s} [kpc];M_{#Delta} [M_{sun}]",             4000, 0, 4e-3, 1000, 0, 1000);
  TH2F *h2d_rs_Mtid       = new TH2F ("h2d_rs_Mtid",       ";scaled radius r_{s} [kpc];M_{tidal} [M_{sun}]",              4000, 0, 4e-3, 1000, 0, 800);
  TH2F *h2d_rs_Rtid       = new TH2F ("h2d_rs_Rtid",       ";scaled radius r_{s} [kpc];R_{tidal} [kpc]",                  4000, 0, 4e-3, 1000, 1e-8, 5e-2);
  TH2F *h2d_rs_Mequdens   = new TH2F ("h2d_rs_Mequdens",   ";scaled radius r_{s} [kpc];M_{equdens} [M_{sun}]",            4000, 0, 4e-3, 1000, 0, 1000);
  TH2F *h2d_rs_Requdens   = new TH2F ("h2d_rs_Requdens",   ";scaled radius r_{s} [kpc];R_{equdens} [kpc]",                4000, 0, 4e-3, 1000, 1e-8, 5e-2);
  TH2F *h2d_rs_J          = new TH2F ("h2d_rs_J",          ";scaled radius r_{s} [kpc];J-factor [GeV^{2}cm^{-5}]",        4000, 0, 4e-3, 1000, 0, 1e8);
  TH2F *h2d_rs_JJcont     = new TH2F ("h2d_rs_JJcont",     ";scaled radius r_{s} [kpc];J/J_{continuum} [GeV^{2}cm^{-5}]", 4000, 0, 4e-3, 1000, 0, 40);


  TH2F *h2d_Glong_Mdelta     = new TH2F ("h2d_Glong_Mdelta",     ";Galactic longitude [deg];M_{#Delta} [M_{sun}]",              5000, -240, 240, 1000, 0, 1000);
  TH2F *h2d_Glong_Mtid       = new TH2F ("h2d_Glong_Mtid",       ";Galactic longitude [deg];M_{tidal} [M_{sun}]",               5000, -240, 240, 1000, 0, 800);
  TH2F *h2d_Glong_Rtid       = new TH2F ("h2d_Glong_Rtid",       ";Galactic longitude [deg];R_{tidal} [kpc]",                   5000, -240, 240, 1000, 1e-8, 5e-2);
  TH2F *h2d_Glong_Mequdens   = new TH2F ("h2d_Glong_Mequdens",   ";Galactic longitude [deg];M_{equdens} [M_{sun}]",             5000, -240, 240, 1000, 0, 1000);
  TH2F *h2d_Glong_Requdens   = new TH2F ("h2d_Glong_Requdens",   ";Galactic longitude [deg];R_{equdens} [kpc]",                 5000, -240, 240, 1000, 1e-8, 5e-2);
  TH2F *h2d_Glong_J          = new TH2F ("h2d_Glong_J",          ";Galactic longitude [deg];J-factor [GeV^{2}cm^{-5}]",         5000, -240, 240, 1000, 0, 1e8);
  TH2F *h2d_Glong_JJcont     = new TH2F ("h2d_Glong_JJcont",     ";Galactic longitude [deg];J/J_{continuum} [GeV^{2}cm^{-5}]",  5000, -240, 240, 1000, 0, 40);

  TH2F *h2d_Glat_Mdelta     = new TH2F ("h2d_Glat_Mdelta",     ";Galactic latitude [deg];M_{#Delta} [M_{sun}]",              5000, -120, 120, 1000, 0, 1000);
  TH2F *h2d_Glat_Mtid       = new TH2F ("h2d_Glat_Mtid",       ";Galactic latitude [deg];M_{tidal} [M_{sun}]",               5000, -120, 120, 1000, 0, 800);
  TH2F *h2d_Glat_Rtid       = new TH2F ("h2d_Glat_Rtid",       ";Galactic latitude [deg];R_{tidal} [kpc]",                   5000, -120, 120, 1000, 1e-8, 5e-2);
  TH2F *h2d_Glat_Mequdens   = new TH2F ("h2d_Glat_Mequdens",   ";Galactic latitude [deg];M_{equdens} [M_{sun}]",             5000, -120, 120, 1000, 0, 1000);
  TH2F *h2d_Glat_Requdens   = new TH2F ("h2d_Glat_Requdens",   ";Galactic latitude [deg];R_{equdens} [kpc]",                 5000, -120, 120, 1000, 1e-8, 5e-2);
  TH2F *h2d_Glat_J          = new TH2F ("h2d_Glat_J",          ";Galactic latitude [deg];J-factor [GeV^{2}cm^{-5}]",         5000, -120, 120, 1000, 0, 1e8);
  TH2F *h2d_Glat_JJcont     = new TH2F ("h2d_Glat_JJcont",     ";Galactic latitude [deg];J/J_{continuum} [GeV^{2}cm^{-5}]",  5000, -120, 120, 1000, 0, 40);

  TH2F *h2d_dist_Mdelta     = new TH2F ("h2d_dist_Mdelta",     ";Distance [kpc];M_{#Delta} [M_{sun}]",               10000, 0, 120, 1000, 0, 1000);
  TH2F *h2d_dist_Mtid       = new TH2F ("h2d_dist_Mtid",       ";Distance [kpc];M_{tidal} [M_{sun}]",                10000, 0, 120, 1000, 0, 800);
  TH2F *h2d_dist_Rtid       = new TH2F ("h2d_dist_Rtid",       ";Distance [kpc];R_{tidal} [kpc]",                    10000, 0, 120, 1000, 1e-8, 5e-2);
  TH2F *h2d_dist_Mequdens   = new TH2F ("h2d_dist_Mequdens",   ";Distance [kpc];M_{equdens} [M_{sun}]",              10000, 0, 120, 1000, 0, 1000);
  TH2F *h2d_dist_Requdens   = new TH2F ("h2d_dist_Requdens",   ";Distance [kpc];R_{equdens} [kpc]",                  10000, 0, 120, 1000, 1e-8, 5e-2);
  TH2F *h2d_dist_J          = new TH2F ("h2d_dist_J",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}]",          10000, 0, 120, 1000, 0, 1e8);
  TH2F *h2d_dist_JJcont     = new TH2F ("h2d_dist_JJcont",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}]",   10000, 0, 120, 1000, 0, 40);

  TH2F *h2d_dist_J_cut0          = new TH2F ("h2d_dist_J_cut0",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}:  0-50]",            10000, 0, 120, 1000, 0, 1e8);
  TH2F *h2d_dist_J_cut1          = new TH2F ("h2d_dist_J_cut1",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}:  50-100]",          10000, 0, 120, 1000, 0, 1e8);
  TH2F *h2d_dist_J_cut2          = new TH2F ("h2d_dist_J_cut2",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 100-200]",          10000, 0, 120, 1000, 0, 1e8);
  TH2F *h2d_dist_J_cut3          = new TH2F ("h2d_dist_J_cut3",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 200-300]",          10000, 0, 120, 1000, 0, 1e8);
  TH2F *h2d_dist_J_cut4          = new TH2F ("h2d_dist_J_cut4",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 300-400]",          10000, 0, 120, 1000, 0, 1e8);
  TH2F *h2d_dist_JJcont_cut0     = new TH2F ("h2d_dist_JJcont_cut0",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 0-100]",     10000, 0, 120, 1000, 0, 40);
  TH2F *h2d_dist_JJcont_cut1     = new TH2F ("h2d_dist_JJcont_cut1",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 0-100]",     10000, 0, 120, 1000, 0, 40);
  TH2F *h2d_dist_JJcont_cut2     = new TH2F ("h2d_dist_JJcont_cut2",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 100-200]",   10000, 0, 120, 1000, 0, 40);
  TH2F *h2d_dist_JJcont_cut3     = new TH2F ("h2d_dist_JJcont_cut3",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 200-300]",   10000, 0, 120, 1000, 0, 40);
  TH2F *h2d_dist_JJcont_cut4     = new TH2F ("h2d_dist_JJcont_cut4",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 300-400]",   10000, 0, 120, 1000, 0, 40);




  TH2F *h2d_Mdelta_rhos       = new TH2F ("h2d_Mdelta_rhos",       ";M_{#Delta} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",        20000, 0, 200, 1000, 1e8, 2e9);
  TH2F *h2d_Mdelta_rs         = new TH2F ("h2d_Mdelta_rs",         ";M_{#Delta} [M_{sun}];r_{s} [kpc]",                       20000, 0, 200, 1000, 0, 1e-3);
  TH2F *h2d_Mdelta_Mtid       = new TH2F ("h2d_Mdelta_Mtid",       ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}]",               20000, 0, 200, 1000, 0, 800);
  TH2F *h2d_Mdelta_Rtid       = new TH2F ("h2d_Mdelta_Rtid",       ";M_{#Delta} [M_{sun}];R_{tidal} [kpc]",                   20000, 0, 200, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mdelta_Mequdens   = new TH2F ("h2d_Mdelta_Mequdens",   ";M_{#Delta} [M_{sun}];M_{equdens} [M_{sun}]",             20000, 0, 200, 1000, 0, 1000);
  TH2F *h2d_Mdelta_Requdens   = new TH2F ("h2d_Mdelta_Requdens",   ";M_{#Delta} [M_{sun}];R_{equdens} [kpc]",                 20000, 0, 200, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mdelta_J          = new TH2F ("h2d_Mdelta_J",          ";M_{#Delta} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",         20000, 0, 200, 1000, 0, 1e8);
  TH2F *h2d_Mdelta_JJcont     = new TH2F ("h2d_Mdelta_JJcont",     ";M_{#Delta} [M_{sun}];J/J_{continuum} [GeV^{2}cm^{-5}]",  20000, 0, 200, 1000, 0, 20);

  TH2F *h2d_Mtid_rhos       = new TH2F ("h2d_Mtid_rhos",       ";M_{tidal} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",       20000, 0, 1000, 1000, 1e8, 2e9);
  TH2F *h2d_Mtid_rs         = new TH2F ("h2d_Mtid_rs",         ";M_{tidal} [M_{sun}];r_{s} [kpc]",                      20000, 0, 1000, 1000, 0, 1e-3);
  TH2F *h2d_Mtid_Mdelta     = new TH2F ("h2d_Mtid_Mdelta",     ";M_{tidal} [M_{sun}];M_{#Delta} [M_{sun}]",             20000, 0, 1000, 1000, 0, 800);
  TH2F *h2d_Mtid_Rtid       = new TH2F ("h2d_Mtid_Rtid",       ";M_{tidal} [M_{sun}];R_{tidal} [kpc]",                  20000, 0, 1000, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mtid_Mequdens   = new TH2F ("h2d_Mtid_Mequdens",   ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}]",            20000, 0, 1000, 1000, 0, 1000);
  TH2F *h2d_Mtid_Requdens   = new TH2F ("h2d_Mtid_Requdens",   ";M_{tidal} [M_{sun}];R_{equdens} [kpc]",                20000, 0, 1000, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mtid_J          = new TH2F ("h2d_Mtid_J",          ";M_{tidal} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",        20000, 0, 1000, 1000, 0, 1e8);
  TH2F *h2d_Mtid_JJcont     = new TH2F ("h2d_Mtid_JJcont",     ";M_{tidal} [M_{sun}];J/J_{continuum} [GeV^{2}cm^{-5}]", 20000, 0, 1000, 1000, 0, 20);

  TH2F *h2d_Rtid_rhos       = new TH2F ("h2d_Rtid_rhos",       ";R_{tidal} [kpc];#rho_{s} [M_{sun}/kpc^{3}]",        40000, 0.002, 0.05, 1000, 1e8, 2e9);
  TH2F *h2d_Rtid_rs         = new TH2F ("h2d_Rtid_rs",         ";R_{tidal} [kpc];r_{s} [kpc]",                       40000, 0.002, 0.05, 1000, 0, 1e-3);
  TH2F *h2d_Rtid_Mdelta     = new TH2F ("h2d_Rtid_Mdelta",     ";R_{tidal} [kpc];M_{#Delta} [M_{sun}]",              40000, 0.002, 0.05, 1000, 0, 800);
  TH2F *h2d_Rtid_Mtid       = new TH2F ("h2d_Rtid_Mtid",       ";R_{tidal} [kpc];M_{tidal} [M_{sun}]",               40000, 0.002, 0.05, 1000, 0, 1000);
  TH2F *h2d_Rtid_Mequdens   = new TH2F ("h2d_Rtid_Mequdens",   ";R_{tidal} [kpc];M_{equdens} [M_{sun}]",             40000, 0.002, 0.05, 1000, 0, 1000);
  TH2F *h2d_Rtid_Requdens   = new TH2F ("h2d_Rtid_Requdens",   ";R_{tidal} [kpc];R_{equdens} [kpc]",                 40000, 0.002, 0.05, 1000, 1e-8, 5e-2);
  TH2F *h2d_Rtid_J          = new TH2F ("h2d_Rtid_J",          ";R_{tidal} [kpc];J-factor [GeV^{2}cm^{-5}]",         40000, 0.002, 0.05, 1000, 0, 1e8);
  TH2F *h2d_Rtid_JJcont     = new TH2F ("h2d_Rtid_JJcont",     ";R_{tidal} [kpc];J/J_{continuum} [GeV^{2}cm^{-5}]",  40000, 0.002, 0.05, 1000, 0, 20);

  TH2F *h2d_Requdens_rhos       = new TH2F ("h2d_Requdens_rhos",       ";R_{equdensal} [kpc];#rho_{s} [M_{sun}/kpc^{3}]",       40000, 0.002, 0.05, 1000, 1e8, 2e9);
  TH2F *h2d_Requdens_rs         = new TH2F ("h2d_Requdens_rs",         ";R_{equdensal} [kpc];r_{s} [kpc]",                      40000, 0.002, 0.05, 1000, 0, 1e-3);
  TH2F *h2d_Requdens_Mdelta     = new TH2F ("h2d_Requdens_Mdelta",     ";R_{equdensal} [kpc];M_{#Delta} [M_{sun}]",             40000, 0.002, 0.05, 1000, 0, 1000);
  TH2F *h2d_Requdens_Mequdens   = new TH2F ("h2d_Requdens_Mequdens",   ";R_{equdensal} [kpc];M_{equdens} [M_{sun}]",            40000, 0.002, 0.05, 1000, 0, 1000);
  TH2F *h2d_Requdens_Rtid       = new TH2F ("h2d_Requdens_Rtid",       ";R_{equdensal} [kpc];R_{tidal} [kpc]",                  40000, 0.002, 0.05, 1000, 1e-8, 5e-2);
  TH2F *h2d_Requdens_Mtid       = new TH2F ("h2d_Requdens_Mtid",       ";R_{equdensal} [kpc];M_{tidal} [M_{sun}]",              40000, 0.002, 0.05, 1000, 0, 800);
  TH2F *h2d_Requdens_J          = new TH2F ("h2d_Requdens_J",          ";R_{equdensal} [kpc];J-factor [GeV^{2}cm^{-5}]",        40000, 0.002, 0.05, 1000, 0, 1e8);
  TH2F *h2d_Requdens_JJcont     = new TH2F ("h2d_Requdens_JJcont",     ";R_{equdensal} [kpc];J/J_{continuum} [GeV^{2}cm^{-5}]", 40000, 0.002, 0.05, 1000, 0, 20);

  TH2F *h2d_Mequdens_rhos       = new TH2F ("h2d_Mequdens_rhos",       ";M_{equdensal} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",       40000, 0, 500, 1000, 1e8, 2e9);
  TH2F *h2d_Mequdens_rs         = new TH2F ("h2d_Mequdens_rs",         ";M_{equdensal} [M_{sun}];r_{s} [kpc]",                      40000, 0, 500, 1000, 0, 1e-3);
  TH2F *h2d_Mequdens_Mdelta     = new TH2F ("h2d_Mequdens_Mdelta",     ";M_{equdensal} [M_{sun}];M_{#Delta} [M_{sun}]",             40000, 0, 500, 1000, 0, 1000);
  TH2F *h2d_Mequdens_Requdens   = new TH2F ("h2d_Mequdens_Requdens",   ";M_{equdensal} [M_{sun}];R_{equdens} [kpc]",                40000, 0, 500, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mequdens_Mtid       = new TH2F ("h2d_Mequdens_Mtid",       ";M_{equdensal} [M_{sun}];M_{tidal} [M_{sun}]",              40000, 0, 500, 1000, 0, 800);
  TH2F *h2d_Mequdens_Rtid       = new TH2F ("h2d_Mequdens_Rtid",       ";M_{equdensal} [M_{sun}];R_{tidal} [kpc]",                  40000, 0, 500, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mequdens_J          = new TH2F ("h2d_Mequdens_J",          ";M_{equdensal} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",        40000, 0, 500, 1000, 0, 1e8);
  TH2F *h2d_Mequdens_JJcont     = new TH2F ("h2d_Mequdens_JJcont",     ";M_{equdensal} [M_{sun}];J/J_{continuum} [GeV^{2}cm^{-5}]", 40000, 0, 500, 1000, 0, 20);


  TProfile2D *prf2d_dist_J_Mdelta        = new TProfile2D ("prf2d_dist_J_Mdelta",         ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{#Delta} [M_{sun}]",               2000, 0, 120, 1000, 0, 1e8, 0, 1000);
  TProfile2D *prf2d_dist_J_Mtidal        = new TProfile2D ("prf2d_dist_J_Mtidal",         ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}]",                 2000, 0, 120, 1000, 0, 1e8, 0, 800);
  TProfile2D *prf2d_dist_J_Mequdens      = new TProfile2D ("prf2d_dist_J_Mequdens",       ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{equdens}[M_{sun}]",               2000, 0, 120, 1000, 0, 1e8, 0, 800);
  TProfile2D *prf2d_dist_J_Rtidal        = new TProfile2D ("prf2d_dist_J_Rtidal",         ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];R_{tidal}[kpc]",                     2000, 0, 120, 1000, 0, 1e8, 1e-8, 5e-2);
  TProfile2D *prf2d_dist_J_Requdens      = new TProfile2D ("prf2d_dist_J_Requdens",       ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];R_{equdens}[kpc]",                   2000, 0, 120, 1000, 0, 1e8, 1e-8, 5e-2);
  TProfile2D *prf2d_dist_J_rhos          = new TProfile2D ("prf2d_dist_J_rhos",           ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];#rho_{s} [M_{sun}/kpc^{3}]",         2000, 0, 120, 1000, 0, 1e8, 0, 1.0e9);
  TProfile2D *prf2d_dist_J_Glong         = new TProfile2D ("prf2d_dist_J_Glong",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];Galactic longitude [deg]",           2000, 0, 120, 1000, 0, 1e8, -240, 240);
  TProfile2D *prf2d_dist_J_Glat          = new TProfile2D ("prf2d_dist_J_Glat",           ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];Galactic latitude [deg]",            2000, 0, 120, 1000, 0, 1e8, -120, 120);
  TProfile2D *prf2d_dist_J_GlongHC       = new TProfile2D ("prf2d_dist_J_GlongHC",        ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];Galactic longitude H.C. [deg]",      2000, 0, 120, 1000, 0, 1e8, -240, 240);
  TProfile2D *prf2d_dist_J_GlatHC        = new TProfile2D ("prf2d_dist_J_GlatHC",         ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];Galactic latitude H.C. [deg]",       2000, 0, 120, 1000, 0, 1e8, -120, 120);

  TProfile2D *prf2d_dist_JJcont_Mtidal   = new TProfile2D ("prf2d_dist_JJcont_Mtidal",    ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}]",          2000, 0, 120, 1000, 0, 40, 0, 800);
  TProfile2D *prf2d_dist_JJcont_Mequdens = new TProfile2D ("prf2d_dist_JJcont_Mequdens",  ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{equdens}[M_{sun}]",        2000, 0, 120, 1000, 0, 40, 0, 800);
  TProfile2D *prf2d_dist_JJcont_Mdelta   = new TProfile2D ("prf2d_dist_JJcont_Mdelta",    ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{#Delta} [M_{sun}]",        2000, 0, 120, 1000, 0, 40, 0, 1000);
  TProfile2D *prf2d_dist_JJcont_Rtidal   = new TProfile2D ("prf2d_dist_JJcont_Rtidal",    ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];R_{tidal}[kpc]",              2000, 0, 120, 1000, 0, 40, 1e-8, 5e-2);
  TProfile2D *prf2d_dist_JJcont_Requdens = new TProfile2D ("prf2d_dist_JJcont_Requdens",  ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];R_{equdens}[kpc]",            2000, 0, 120, 1000, 0, 40, 1e-8, 5e-2);
  TProfile2D *prf2d_dist_JJcont_rhos     = new TProfile2D ("prf2d_dist_JJcont_rhos",      ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];#rho_{s} [M_{sun}/kpc^{3}]",  2000, 0, 120, 1000, 0, 40, 0, 1.0e9);
  TProfile2D *prf2d_dist_JJcont_Glong    = new TProfile2D ("prf2d_dist_JJcont_Glong",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];Galactic longitude [deg]",    2000, 0, 120, 1000, 0, 40, -240, 240);
  TProfile2D *prf2d_dist_JJcont_Glat     = new TProfile2D ("prf2d_dist_JJcont_Glat",      ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];Galactic latitude [deg]",     2000, 0, 120, 1000, 0, 40, -120, 120);
  TProfile2D *prf2d_dist_JJcont_GlongHC  = new TProfile2D ("prf2d_dist_JJcont_GlongHC",   ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];Galactic longitude H.C. [deg]",2000, 0, 120, 1000, 0, 40, -240, 240);
  TProfile2D *prf2d_dist_JJcont_GlatHC   = new TProfile2D ("prf2d_dist_JJcont_GlatHC",    ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];Galactic latitude H.C. [deg]", 2000, 0, 120, 1000, 0, 40, -120, 120);

  TProfile2D *prf2d_Glon_Glat_Mdelta     = new TProfile2D ("prf2d_Glon_Glat_Mdelta",      ";Galactic longitude [deg];Galactic latitude [deg];M_{#Delta} [M_{sun}]",       1000, -200, 200, 1000, -100, 100, 0, 1000);
  TProfile2D *prf2d_Glon_Glat_Mtidal     = new TProfile2D ("prf2d_Glon_Glat_Mtidal",      ";Galactic longitude [deg];Galactic latitude [deg];M_{tidal} [M_{sun}]",        1000, -200, 200, 1000, -100, 100, 0, 1000);
  TProfile2D *prf2d_Glon_Glat_Mequdens   = new TProfile2D ("prf2d_Glon_Glat_Mequdens",    ";Galactic longitude [deg];Galactic latitude [deg];M_{equdens} [M_{sun}]",      1000, -200, 200, 1000, -100, 100, 0, 1000);
  TProfile2D *prf2d_Glon_Glat_Rtidal     = new TProfile2D ("prf2d_Glon_Glat_Rtidal",      ";Galactic longitude [deg];Galactic latitude [deg];R_{tidal} [kpc]",            1000, -200, 200, 1000, -100, 100, 1e-8, 5e-2);
  TProfile2D *prf2d_Glon_Glat_Requdens   = new TProfile2D ("prf2d_Glon_Glat_Requdens",    ";Galactic longitude [deg];Galactic latitude [deg];R_{equdens} [kpc]",          1000, -200, 200, 1000, -100, 100, 1e-8, 5e-2);
  TProfile2D *prf2d_Glon_Glat_J          = new TProfile2D ("prf2d_Glon_Glat_J",           ";Galactic longitude [deg];Galactic latitude [deg];J-factor [GeV^{2}cm^{-5}]",  1000, -200, 200, 1000, -100, 100, 0, 1e8);
  TProfile2D *prf2d_Glon_Glat_JJcont     = new TProfile2D ("prf2d_Glon_Glat_JJcont",      ";Galactic longitude [deg];Galactic latitude [deg];J/J_{cont} [GeV^{2}cm^{-5}]",1000, -200, 200, 1000, -100, 100, 0, 40);
  TProfile2D *prf2d_Glon_Glat_rhos       = new TProfile2D ("prf2d_Glon_Glat_rhos",        ";Galactic longitude [deg];Galactic latitude [deg];#rho_{s} [M_{sun}/kpc^{3}]", 1000, -200, 200, 1000, -100, 100, 0, 1.5e9);
  TProfile2D *prf2d_Glon_Glat_rs         = new TProfile2D ("prf2d_Glon_Glat_rs",          ";Galactic longitude [deg];Galactic latitude [deg];r_{s} [kpc]",                1000, -200, 200, 1000, -100, 100, 0, 4e-3);
  TProfile2D *prf2d_GlonHC_GlatHC_J      = new TProfile2D ("prf2d_GlonHC_GlatHC_J",       ";Galactic longitude w.r.t. to HC [deg];Galactic latitude w.r.t. to HC [deg];J-factor [GeV^{2}cm^{-5}]",  2000, -200, 200, 2000, -100, 100, 0, 1e8);
  TProfile2D *prf2d_GlonHC_GlatHC_JJcont = new TProfile2D ("prf2d_GlonHC_GlatHC_JJcont",  ";Galactic longitude w.r.t. to HC [deg];Galactic latitude w.r.t. to HC [deg];J/J_{cont} [GeV^{2}cm^{-5}]",2000, -200, 200, 2000, -100, 100, 0, 40);
  TProfile2D *prf2d_GlonHC_GlatHC_rhos   = new TProfile2D ("prf2d_GlonHC_GlatHC_rhos",    ";Galactic longitude w.r.t. to HC [deg];Galactic latitude w.r.t. to HC [deg];#rho_{s} [M_{sun}/kpc^{3}]", 2000, -200, 200, 2000, -100, 100, 0, 1.5e9);
  TProfile2D *prf2d_GlonHC_GlatHC_rs     = new TProfile2D ("prf2d_Glon_Glat_rs",          ";Galactic longitude w.r.t. to HC [deg];Galactic latitude w.r.t. to HC [deg];r_{s} [kpc]",                2000, -200, 200, 2000, -100, 100, 0, 4e-3);


  TProfile2D *prf2d_Mdelta_Mtidal_rs         = new TProfile2D ("prf2d_Mdelta_Mtidal_rs",          ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 800, 0, 4e-3);
  TProfile2D *prf2d_Mdelta_Mtidal_rhos       = new TProfile2D ("prf2d_Mdelta_Mtidal_rhos",        ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 800, 0, 1.5e9);
  TProfile2D *prf2d_Mdelta_Mtidal_Rtidal     = new TProfile2D ("prf2d_Mdelta_Mtidal_Rtidal",      ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 800, 1e-8, 5e-2);
  TProfile2D *prf2d_Mdelta_Mtidal_Mequdens   = new TProfile2D ("prf2d_Mdelta_Mtidal_Mequdens",    ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];M_{equdens} [M_{sun}]",         1000,  0, 800, 1000, 0, 800, 0, 800);
  TProfile2D *prf2d_Mdelta_Mtidal_Requdens   = new TProfile2D ("prf2d_Mdelta_Mtidal_Requdens",    ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 800, 1e-8, 5e-2);
  TProfile2D *prf2d_Mdelta_Mtidal_Jfact      = new TProfile2D ("prf2d_Mdelta_Mtidal_Jfact",       ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 800, 0, 1e8);
  TProfile2D *prf2d_Mdelta_Mtidal_JJ         = new TProfile2D ("prf2d_Mdelta_Mtidal_JJ",          ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 800, 0, 40);
  TProfile2D *prf2d_Mdelta_Mtidal_dist       = new TProfile2D ("prf2d_Mdelta_Mtidal_dist",        ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];distance [kpc]",                1000,  0, 800, 1000, 0, 800, 0, 120);
  TProfile2D *prf2d_Mdelta_Mtidal_Glong      = new TProfile2D ("prf2d_Mdelta_Mtidal_Glong",       ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 800, -200, 200);
  TProfile2D *prf2d_Mdelta_Mtidal_Glat       = new TProfile2D ("prf2d_Mdelta_Mtidal_Glat",        ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 800, -100, 100);

  TProfile2D *prf2d_Mtidal_Mequdens_rs         = new TProfile2D ("prf2d_Mtidal_Mequdens_rs",          ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 800, 0, 4e-3);
  TProfile2D *prf2d_Mtidal_Mequdens_rhos       = new TProfile2D ("prf2d_Mtidal_Mequdens_rhos",        ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 800, 0, 1.5e9);
  TProfile2D *prf2d_Mtidal_Mequdens_Mdelta     = new TProfile2D ("prf2d_Mtidal_Mequdens_Mdelta",      ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];M_{#Delta} [M_{sun}]",          1000,  0, 800, 1000, 0, 800, 0, 800);
  TProfile2D *prf2d_Mtidal_Mequdens_Rtidal     = new TProfile2D ("prf2d_Mtidal_Mequdens_Rtidal",      ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 800, 1e-8, 5e-2);
  TProfile2D *prf2d_Mtidal_Mequdens_Requdens   = new TProfile2D ("prf2d_Mtidal_Mequdens_Requdens",    ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 800, 1e-8, 5e-2);
  TProfile2D *prf2d_Mtidal_Mequdens_Jfact      = new TProfile2D ("prf2d_Mtidal_Mequdens_Jfact",       ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 800, 0, 1e8);
  TProfile2D *prf2d_Mtidal_Mequdens_JJ         = new TProfile2D ("prf2d_Mtidal_Mequdens_JJ",          ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 800, 0, 40);
  TProfile2D *prf2d_Mtidal_Mequdens_dist       = new TProfile2D ("prf2d_Mtidal_Mequdens_dist",        ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];distance [kpc]",                1000,  0, 800, 1000, 0, 800, 0, 120);
  TProfile2D *prf2d_Mtidal_Mequdens_Glong      = new TProfile2D ("prf2d_Mtidal_Mequdens_Glong",       ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 800, -200, 200);
  TProfile2D *prf2d_Mtidal_Mequdens_Glat       = new TProfile2D ("prf2d_Mtidal_Mequdens_Glat",        ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 800, -100, 100);




  string DMclumID;
  string DMtype;
  string Glong, Glat, DM_clump_dist;
  string z;
  string Rdelta, rhos, rs;
  string DMprof;

  string DMprof_p1, DMprof_p2, DMprof_p3;
  string DM_J, DM_JJcont;
  string Mdelta, Mtid, Rtid, Mequdens, Requdens, Dgal;

  string Glong_wrtHaloCenter;
  string Glat_wrtHaloCenter;

  string dummy;



  // Import data from a file.
  Int_t counter = 0;

  std::string line;
  ifstream energy_dat;

  energy_dat.open(file);

std::string parsed;
std::stringstream input_stringstream;
/**/
  while (getline(energy_dat, line)) {
        getline(energy_dat, line);
        input_stringstream << line;


        std::stringstream ss(line);
        std::cout << "line# " << counter    << std::endl;
        std::cout << "ss: " << line << std::endl;
        std::cout << "length: = " << line.length()    << std::endl;

        //if (counter>3 && line.length() < 184){
        if (counter>1 && line.length() < 160){
          //break;
          if (counter>=49840 && counter<478138){
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                //std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,DMclumID,' ');   //std::cout << "DMclumID: " << DMclumID << " --> length:" << DMclumID.length()  << std::endl;
                if (DMclumID.length() == 0){
                  std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                  if (DMclumID.length() == 0){
                    std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                    if (DMclumID.length() == 0){
                       std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;

                        if (DMclumID.length() == 0){
                           std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                    }
                  }
                 }
                }
            }
          else if (counter>478137){
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                //std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                //std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,DMclumID,' ');   //std::cout << "DMclumID: " << DMclumID << " --> length:" << DMclumID.length()  << std::endl;
                if (DMclumID.length() == 0){
                  std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                  if (DMclumID.length() == 0){
                    std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                    if (DMclumID.length() == 0){
                       std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;

                        if (DMclumID.length() == 0){
                           std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                    }
                  }
                 }
                }
          }
            else{
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;

                std::getline(ss,DMclumID,' ');   //std::cout << "DMclumID: " << DMclumID << " --> length:" << DMclumID.length()  << std::endl;
                if (DMclumID.length() == 0){
                  std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                  if (DMclumID.length() == 0){
                    std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                    if (DMclumID.length() == 0){
                       std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;

                        if (DMclumID.length() == 0){
                           std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                    }
                  }
                 }
                }

            }
                std::cout << "DMclumID: " << DMclumID << std::endl;


                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,DMtype,' ');     std::cout << "DMtype: " << DMtype << " --> length:" << DMtype.length()  << std::endl;

                //std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << " --> length:" << dummy.length()  << std::endl;
                std::getline(ss,Glong,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                if (Glong.length() == 0){
                  std::getline(ss,Glong,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                  if (Glong.length() == 0){
                    std::getline(ss,Glong,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                    if (Glong.length() == 0)
                       std::getline(ss,Glong,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;
                  }

                }
                std::cout << "Glong: " << Glong << std::endl;



                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << " --> length:" << dummy.length()  << std::endl;

                std::getline(ss,Glat,' ');      //std::cout << "Glat: " << Glat << " --> length:" << Glat.length()  << std::endl;
                if (Glat.length() == 0){
                  std::getline(ss,Glat,' ');      //std::cout << "Glat: " << Glat << " --> length:" << Glat.length()  << std::endl;

                  if (Glat.length() == 0)
                    std::getline(ss,Glat,' ');      //std::cout << "Glat: " << Glat << " --> length:" << Glat.length()  << std::endl;
                }
                std::cout << "Glat: " << Glat  << std::endl;

                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,DM_clump_dist,' ');     //std::cout << "DM_clump_dist: " << DM_clump_dist  << std::endl;
                //std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,z,' ');     //std::cout << "z: " << z  << std::endl;
                //std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,Rdelta,' ');     //std::cout << "Rdelta: " << Rdelta  << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,rhos,' ');       //std::cout << "rhos: " << rhos  << std::endl;
                std::getline(ss,rs,' ');         //std::cout << "rs: " << rs  << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,DMprof,' ');     //std::cout << "DMprof: " << DMprof  << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,DMprof_p1,' ');   //std::cout << "p1: " << DMprof_p1  << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,DMprof_p2,' ');   //std::cout << "p2: " << DMprof_p2  << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,DMprof_p3,' ');   //std::cout << "p3: " << DMprof_p3  << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,DM_J,' ');        //std::cout << "J: " << DM_J  << std::endl;
                std::getline(ss,DM_JJcont,' ');   //std::cout << "JJcont: " << DM_JJcont  << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,Mdelta,' ');     //std::cout << "Mdelta: " << Mdelta  << std::endl;

                if (FLAG==1){
                std::getline(ss,Mtid,' ');       //std::cout << "Mtid: " << Mtid  << std::endl;
                std::getline(ss,Rtid,' ');       //std::cout << "Rtid: " << Rtid  << std::endl;
                std::getline(ss,Mequdens,' ');   //std::cout << "Mequdens: " << Mequdens  << std::endl;
                std::getline(ss,Requdens,' ');   //std::cout << "Requdens: " << Requdens  << std::endl;
                std::getline(ss,Dgal,' ');       //std::cout << "Dgal: " << Dgal  << std::endl;
                }
                if (FLAG==2){
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,Glong_wrtHaloCenter,' ');       //std::cout << "Mtid: " << Mtid  << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                std::getline(ss,Glat_wrtHaloCenter,' ');       //std::cout << "Rtid: " << Rtid  << std::endl;
                }
                std::cout << "Glong (w.r.t. haloCenter): " << Glong_wrtHaloCenter  << std::endl;
                std::cout << "Glat (w.r.t. haloCenter): " << Glat_wrtHaloCenter  << std::endl;



                int dm_ID      = std::atoi(DMclumID.c_str());
                float dm_glong = std::stof(Glong);
                float dm_glat  = std::stof(Glat);
                float dm_dist  = std::stof(DM_clump_dist);


                float dm_z      = std::stof(z);
                float dm_Rdelta = std::stof(Rdelta);
                float dm_rhos   = std::stof(rhos);
                float dm_rs     = std::stof(rs);

                float dm_p1     = std::stof(DMprof_p1);
                float dm_p2     = std::stof(DMprof_p2);
                float dm_p3     = std::stof(DMprof_p3);

                float dm_Jfactor = std::stof(DM_J);
                float dm_JJcont  = std::stof(DM_JJcont);


                float dm_Mdelta = 0.;//std::stof(Mdelta);
                float dm_Mtid   = 0.;//std::stof(Mtid);
                float dm_Rtid   = 0.;//std::stof(Rtid);

                float dm_Mequdens = 0.;//std::stof(Mequdens);
                float dm_Requdens = 0.;//std::stof(Requdens);
                float dm_Dgal     = 0.;//std::stof(Dgal);

                if (FLAG==1){
                  dm_Mdelta = std::stof(Mdelta);
                  dm_Mtid   = std::stof(Mtid);
                  dm_Rtid   = std::stof(Rtid);

                  dm_Mequdens = std::stof(Mequdens);
                  dm_Requdens = std::stof(Requdens);
                  dm_Dgal     = std::stof(Dgal);
                }

                float dm_Glong_wrtHaloCenter   = 0.0;
                float dm_Glat_wrtHaloCenter    = 0.0;
                if (FLAG==2){
                   dm_Glong_wrtHaloCenter   = std::stof(Glong_wrtHaloCenter);
                   dm_Glat_wrtHaloCenter    = std::stof(Glat_wrtHaloCenter);
                }


                 std::cout << "dm_glong: " << dm_glong  << std::endl;
                 std::cout << "dm_glat: " << dm_glat  << std::endl;
                 std::cout << "dm_dist: " << dm_dist  << std::endl;
                 std::cout << "dm_Rdelta: " << dm_Rdelta  << std::endl;
                 std::cout << "dm_rhos: " << dm_rhos  << std::endl;
                 std::cout << "dm_rs: " << dm_rs  << std::endl;


                mGlong = dm_glong;
                mGlat  = dm_glat;
                mDist  = dm_dist;

                mZ      = dm_z;
                mRdelta = dm_Rdelta;
                mRhos   = dm_rhos;
                mRhos2  = mRhos*mRhos;
                mRs     = dm_rs;

                mParam1    = dm_p1;
                mParam2    = dm_p2;
                mParam3    = dm_p3;
                mJfactor   = dm_Jfactor;
                mJJcont    = dm_JJcont;

                mMdelta    = dm_Mdelta;
                mMtidal    = dm_Mtid;
                mRtidal    = dm_Rtid;
                mMequdens  = dm_Mequdens;
                mRequdens  = dm_Requdens;
                mDgal      = dm_Dgal;

                mGlong_wrtHaloCenter = dm_Glong_wrtHaloCenter;
                mGlat_wrtHaloCenter  = dm_Glat_wrtHaloCenter;


                prf_rhos_Mdelta   -> Fill (mRhos, mMdelta);
                prf_rhos_Mtid     -> Fill (mRhos, mMtidal);
                prf_rhos_Rtid     -> Fill (mRhos, mRtidal);
                prf_rhos_Mequdens -> Fill (mRhos, mMequdens);
                prf_rhos_Requdens -> Fill (mRhos, mRequdens);
                prf_rhos_J        -> Fill (mRhos, mJfactor);
                prf_rhos_JJcont   -> Fill (mRhos, mJJcont);
                prf_rhos_dist     -> Fill (mRhos, mDist);

                prf_rs_Mdelta   -> Fill (mRs, mMdelta);
                prf_rs_Mtid     -> Fill (mRs, mMtidal);
                prf_rs_Rtid     -> Fill (mRs, mRtidal);
                prf_rs_Mequdens -> Fill (mRs, mMequdens);
                prf_rs_Requdens -> Fill (mRs, mRequdens);
                prf_rs_J        -> Fill (mRs, mJfactor);
                prf_rs_JJcont   -> Fill (mRs, mJJcont);
                prf_rs_rhos     -> Fill (mRs, mRhos);
                prf_rs_rhos2    -> Fill (mRs, mRhos2);
                prf_rs_dist     -> Fill (mRs, mDist);

                prf_Glong_Mdelta   -> Fill (mGlong, mMdelta);
                prf_Glong_Mtid     -> Fill (mGlong, mMtidal);
                prf_Glong_Rtid     -> Fill (mGlong, mRtidal);
                prf_Glong_Mequdens -> Fill (mGlong, mMequdens);
                prf_Glong_Requdens -> Fill (mGlong, mRequdens);
                prf_Glong_J        -> Fill (mGlong, mJfactor);
                prf_Glong_JJcont   -> Fill (mGlong, mJJcont);

                prf_Glat_Mdelta   -> Fill (mGlat, mMdelta);
                prf_Glat_Mtid     -> Fill (mGlat, mMtidal);
                prf_Glat_Rtid     -> Fill (mGlat, mRtidal);
                prf_Glat_Mequdens -> Fill (mGlat, mMequdens);
                prf_Glat_Requdens -> Fill (mGlat, mRequdens);
                prf_Glat_J        -> Fill (mGlat, mJfactor);
                prf_Glat_JJcont   -> Fill (mGlat, mJJcont);

                prf_dist_Mdelta   -> Fill (mDist, mMdelta);
                prf_dist_Mtid     -> Fill (mDist, mMtidal);
                prf_dist_Rtid     -> Fill (mDist, mRtidal);
                prf_dist_Mequdens -> Fill (mDist, mMequdens);
                prf_dist_Requdens -> Fill (mDist, mRequdens);
                prf_dist_J        -> Fill (mDist, mJfactor);
                prf_dist_JJcont   -> Fill (mDist, mJJcont);

                prf_Mdelta_rhos     -> Fill (mMdelta, mRhos);
                prf_Mdelta_rs       -> Fill (mMdelta, mRs);
                prf_Mdelta_Mtid     -> Fill (mMdelta, mMtidal);
                prf_Mdelta_Rtid     -> Fill (mMdelta, mRtidal);
                prf_Mdelta_Mequdens -> Fill (mMdelta, mMequdens);
                prf_Mdelta_Requdens -> Fill (mMdelta, mRequdens);
                prf_Mdelta_J        -> Fill (mMdelta, mJfactor);
                prf_Mdelta_JJcont   -> Fill (mMdelta, mJJcont);

                prf_Mtid_rhos     -> Fill (mMtidal, mRhos);
                prf_Mtid_rs       -> Fill (mMtidal, mRs);
                prf_Mtid_Mdelta   -> Fill (mMtidal, mMdelta);
                prf_Mtid_Rtid     -> Fill (mMtidal, mRtidal);
                prf_Mtid_Mequdens -> Fill (mMtidal, mMequdens);
                prf_Mtid_Requdens -> Fill (mMtidal, mRequdens);
                prf_Mtid_J        -> Fill (mMtidal, mJfactor);
                prf_Mtid_JJcont   -> Fill (mMtidal, mJJcont);

                prf_Rtid_rhos     -> Fill (mRtidal, mRhos);
                prf_Rtid_rs       -> Fill (mRtidal, mRs);
                prf_Rtid_Mtid     -> Fill (mRtidal, mMtidal);
                prf_Rtid_Mdelta   -> Fill (mRtidal, mMdelta);
                prf_Rtid_Mequdens -> Fill (mRtidal, mMequdens);
                prf_Rtid_Requdens -> Fill (mRtidal, mRequdens);
                prf_Rtid_J        -> Fill (mRtidal, mJfactor);
                prf_Rtid_JJcont   -> Fill (mRtidal, mJJcont);


                prf_Mequdens_rhos     -> Fill (mMequdens, mRhos);
                prf_Mequdens_rs       -> Fill (mMequdens, mRs);
                prf_Mequdens_Mdelta   -> Fill (mMequdens, mMdelta);
                prf_Mequdens_Requdens -> Fill (mMequdens, mRequdens);
                prf_Mequdens_Rtid     -> Fill (mMequdens, mRtidal);
                prf_Mequdens_Mtid     -> Fill (mMequdens, mMtidal);
                prf_Mequdens_J        -> Fill (mMequdens, mJfactor);
                prf_Mequdens_JJcont   -> Fill (mMequdens, mJJcont);

                prf_Requdens_rhos     -> Fill (mRequdens, mRhos);
                prf_Requdens_rs       -> Fill (mRequdens, mRs);
                prf_Requdens_Mtid     -> Fill (mRequdens, mMtidal);
                prf_Requdens_Rtid     -> Fill (mRequdens, mRtidal);
                prf_Requdens_Mdelta   -> Fill (mRequdens, mMdelta);
                prf_Requdens_Mequdens -> Fill (mRequdens, mMequdens);
                prf_Requdens_J        -> Fill (mRequdens, mJfactor);
                prf_Requdens_JJcont   -> Fill (mRequdens, mJJcont);





                h2d_rhos_Mdelta   -> Fill (mRhos, mMdelta);
                h2d_rhos_Mtid     -> Fill (mRhos, mMtidal);
                h2d_rhos_Rtid     -> Fill (mRhos, mRtidal);
                h2d_rhos_Mequdens -> Fill (mRhos, mMequdens);
                h2d_rhos_Requdens -> Fill (mRhos, mRequdens);
                h2d_rhos_J        -> Fill (mRhos, mJfactor);
                h2d_rhos_JJcont   -> Fill (mRhos, mJJcont);

                h2d_rs_Mdelta   -> Fill (mRs, mMdelta);
                h2d_rs_Mtid     -> Fill (mRs, mMtidal);
                h2d_rs_Rtid     -> Fill (mRs, mRtidal);
                h2d_rs_Mequdens -> Fill (mRs, mMequdens);
                h2d_rs_Requdens -> Fill (mRs, mRequdens);
                h2d_rs_J        -> Fill (mRs, mJfactor);
                h2d_rs_JJcont   -> Fill (mRs, mJJcont);


                h2d_Glong_Mdelta   -> Fill (mGlong, mMdelta);
                h2d_Glong_Mtid     -> Fill (mGlong, mMtidal);
                h2d_Glong_Rtid     -> Fill (mGlong, mRtidal);
                h2d_Glong_Mequdens -> Fill (mGlong, mMequdens);
                h2d_Glong_Requdens -> Fill (mGlong, mRequdens);
                h2d_Glong_J        -> Fill (mGlong, mJfactor);
                h2d_Glong_JJcont   -> Fill (mGlong, mJJcont);

                h2d_Glat_Mdelta   -> Fill (mGlat, mMdelta);
                h2d_Glat_Mtid     -> Fill (mGlat, mMtidal);
                h2d_Glat_Rtid     -> Fill (mGlat, mRtidal);
                h2d_Glat_Mequdens -> Fill (mGlat, mMequdens);
                h2d_Glat_Requdens -> Fill (mGlat, mRequdens);
                h2d_Glat_J        -> Fill (mGlat, mJfactor);
                h2d_Glat_JJcont   -> Fill (mGlat, mJJcont);

                h2d_dist_Mdelta   -> Fill (mDist, mMdelta);
                h2d_dist_Mtid     -> Fill (mDist, mMtidal);
                h2d_dist_Rtid     -> Fill (mDist, mRtidal);
                h2d_dist_Mequdens -> Fill (mDist, mMequdens);
                h2d_dist_Requdens -> Fill (mDist, mRequdens);
                h2d_dist_J        -> Fill (mDist, mJfactor);
                h2d_dist_JJcont   -> Fill (mDist, mJJcont);

                if (mMtidal>0 && mMtidal<50){
                  h2d_dist_J_cut0        -> Fill (mDist, mJfactor);
                  h2d_dist_JJcont_cut0   -> Fill (mDist, mJJcont);
                }

                if (mMtidal>50 && mMtidal<100){
                  h2d_dist_J_cut1        -> Fill (mDist, mJfactor);
                  h2d_dist_JJcont_cut1   -> Fill (mDist, mJJcont);
                }
                if (mMtidal>100 && mMtidal<200){
                  h2d_dist_J_cut2        -> Fill (mDist, mJfactor);
                  h2d_dist_JJcont_cut2   -> Fill (mDist, mJJcont);
                }
                if (mMtidal>200 && mMtidal<300){
                  h2d_dist_J_cut3        -> Fill (mDist, mJfactor);
                  h2d_dist_JJcont_cut3   -> Fill (mDist, mJJcont);
                }
                if (mMtidal>300 && mMtidal<400){
                  h2d_dist_J_cut4        -> Fill (mDist, mJfactor);
                  h2d_dist_JJcont_cut4   -> Fill (mDist, mJJcont);
                }

                prf2d_dist_J_Mtidal        -> Fill (mDist, mJfactor, mMtidal);
                prf2d_dist_J_Mequdens      -> Fill (mDist, mJfactor, mMequdens);
                prf2d_dist_J_Mdelta        -> Fill (mDist, mJfactor, mMdelta);
                prf2d_dist_J_Rtidal        -> Fill (mDist, mJfactor, mRtidal);
                prf2d_dist_J_Requdens      -> Fill (mDist, mJfactor, mRequdens);
                prf2d_dist_J_rhos          -> Fill (mDist, mJfactor, mRhos);
                prf2d_dist_J_Glong         -> Fill (mDist, mJfactor, mGlong);
                prf2d_dist_J_Glat          -> Fill (mDist, mJfactor, mGlat);
                prf2d_dist_J_GlongHC       -> Fill (mDist, mJfactor, mGlong_wrtHaloCenter);
                prf2d_dist_J_GlatHC        -> Fill (mDist, mJfactor, mGlat_wrtHaloCenter);

                prf2d_dist_JJcont_Mtidal   -> Fill (mDist, mJJcont,  mMtidal);
                prf2d_dist_JJcont_Mequdens -> Fill (mDist, mJJcont,  mMequdens);
                prf2d_dist_JJcont_Mdelta   -> Fill (mDist, mJJcont,  mMdelta);
                prf2d_dist_JJcont_Rtidal   -> Fill (mDist, mJJcont,  mRtidal);
                prf2d_dist_JJcont_Requdens -> Fill (mDist, mJJcont,  mRequdens);
                prf2d_dist_JJcont_rhos     -> Fill (mDist, mJJcont,  mRhos);
                prf2d_dist_JJcont_Glong    -> Fill (mDist, mJJcont,  mGlong);
                prf2d_dist_JJcont_Glat     -> Fill (mDist, mJJcont,  mGlat);
                prf2d_dist_JJcont_GlongHC  -> Fill (mDist, mJJcont, mGlong_wrtHaloCenter);
                prf2d_dist_JJcont_GlatHC   -> Fill (mDist, mJJcont, mGlat_wrtHaloCenter);

                prf2d_Glon_Glat_Mdelta    -> Fill (mGlong, mGlat, mMdelta);
                prf2d_Glon_Glat_Mtidal    -> Fill (mGlong, mGlat, mMtidal);
                prf2d_Glon_Glat_Rtidal    -> Fill (mGlong, mGlat, mRtidal);
                prf2d_Glon_Glat_Mequdens  -> Fill (mGlong, mGlat, mMequdens);
                prf2d_Glon_Glat_Requdens  -> Fill (mGlong, mGlat, mRequdens);
                prf2d_Glon_Glat_J         -> Fill (mGlong, mGlat, mJfactor);
                prf2d_Glon_Glat_JJcont    -> Fill (mGlong, mGlat, mJJcont);
                prf2d_Glon_Glat_rhos      -> Fill (mGlong, mGlat, mRhos);
                prf2d_Glon_Glat_rs        -> Fill (mGlong, mGlat, mRs);
                prf2d_GlonHC_GlatHC_J         -> Fill (mGlong_wrtHaloCenter, mGlat_wrtHaloCenter, mJfactor);
                prf2d_GlonHC_GlatHC_JJcont    -> Fill (mGlong_wrtHaloCenter, mGlat_wrtHaloCenter, mJJcont);
                prf2d_GlonHC_GlatHC_rhos      -> Fill (mGlong_wrtHaloCenter, mGlat_wrtHaloCenter, mRhos);
                prf2d_GlonHC_GlatHC_rs        -> Fill (mGlong_wrtHaloCenter, mGlat_wrtHaloCenter, mRs);



                prf2d_Mdelta_Mtidal_rs       -> Fill (mMdelta, mMtidal, mRs);
                prf2d_Mdelta_Mtidal_rhos     -> Fill (mMdelta, mMtidal, mRhos);
                prf2d_Mdelta_Mtidal_Rtidal   -> Fill (mMdelta, mMtidal, mRtidal);
                prf2d_Mdelta_Mtidal_Mequdens -> Fill (mMdelta, mMtidal, mMequdens);
                prf2d_Mdelta_Mtidal_Requdens -> Fill (mMdelta, mMtidal, mRequdens);
                prf2d_Mdelta_Mtidal_Jfact    -> Fill (mMdelta, mMtidal, mJfactor);
                prf2d_Mdelta_Mtidal_JJ       -> Fill (mMdelta, mMtidal, mJJcont);
                prf2d_Mdelta_Mtidal_dist     -> Fill (mMdelta, mMtidal, mDist);
                prf2d_Mdelta_Mtidal_Glong    -> Fill (mMdelta, mMtidal, mGlong);
                prf2d_Mdelta_Mtidal_Glat     -> Fill (mMdelta, mMtidal, mGlat);

                prf2d_Mtidal_Mequdens_rs       -> Fill (mMtidal, mMequdens, mRs);
                prf2d_Mtidal_Mequdens_rhos     -> Fill (mMtidal, mMequdens, mRhos);
                prf2d_Mtidal_Mequdens_Mdelta   -> Fill (mMtidal, mMequdens, mMdelta);
                prf2d_Mtidal_Mequdens_Rtidal   -> Fill (mMtidal, mMequdens, mRtidal);
                prf2d_Mtidal_Mequdens_Requdens -> Fill (mMtidal, mMequdens, mRequdens);
                prf2d_Mtidal_Mequdens_Jfact    -> Fill (mMtidal, mMequdens, mJfactor);
                prf2d_Mtidal_Mequdens_JJ       -> Fill (mMtidal, mMequdens, mJJcont);
                prf2d_Mtidal_Mequdens_dist     -> Fill (mMtidal, mMequdens, mDist);
                prf2d_Mtidal_Mequdens_Glong    -> Fill (mMtidal, mMequdens, mGlong);
                prf2d_Mtidal_Mequdens_Glat     -> Fill (mMtidal, mMequdens, mGlat);



                h2d_Mdelta_rhos     -> Fill (mMdelta, mRhos);
                h2d_Mdelta_rs       -> Fill (mMdelta, mRs);
                h2d_Mdelta_Mtid     -> Fill (mMdelta, mMtidal);
                h2d_Mdelta_Rtid     -> Fill (mMdelta, mRtidal);
                h2d_Mdelta_Mequdens -> Fill (mMdelta, mMequdens);
                h2d_Mdelta_Requdens -> Fill (mMdelta, mRequdens);
                h2d_Mdelta_J        -> Fill (mMdelta, mJfactor);
                h2d_Mdelta_JJcont   -> Fill (mMdelta, mJJcont);

                h2d_Mtid_rhos     -> Fill (mMtidal, mRhos);
                h2d_Mtid_rs       -> Fill (mMtidal, mRs);
                h2d_Mtid_Mdelta   -> Fill (mMtidal, mMdelta);
                h2d_Mtid_Rtid     -> Fill (mMtidal, mRtidal);
                h2d_Mtid_Mequdens -> Fill (mMtidal, mMequdens);
                h2d_Mtid_Requdens -> Fill (mMtidal, mRequdens);
                h2d_Mtid_J        -> Fill (mMtidal, mJfactor);
                h2d_Mtid_JJcont   -> Fill (mMtidal, mJJcont);

                h2d_Rtid_rhos     -> Fill (mRtidal, mRhos);
                h2d_Rtid_rs       -> Fill (mRtidal, mRs);
                h2d_Rtid_Mtid     -> Fill (mRtidal, mMtidal);
                h2d_Rtid_Mdelta   -> Fill (mRtidal, mMdelta);
                h2d_Rtid_Mequdens -> Fill (mRtidal, mMequdens);
                h2d_Rtid_Requdens -> Fill (mRtidal, mRequdens);
                h2d_Rtid_J        -> Fill (mRtidal, mJfactor);
                h2d_Rtid_JJcont   -> Fill (mRtidal, mJJcont);

                h2d_Mequdens_rhos     -> Fill (mMequdens, mRhos);
                h2d_Mequdens_rs       -> Fill (mMequdens, mRs);
                h2d_Mequdens_Mdelta   -> Fill (mMequdens, mMdelta);
                h2d_Mequdens_Requdens -> Fill (mMequdens, mRequdens);
                h2d_Mequdens_Rtid     -> Fill (mMequdens, mRtidal);
                h2d_Mequdens_Mtid     -> Fill (mMequdens, mMtidal);
                h2d_Mequdens_J        -> Fill (mMequdens, mJfactor);
                h2d_Mequdens_JJcont   -> Fill (mMequdens, mJJcont);

                h2d_Requdens_rhos     -> Fill (mRequdens, mRhos);
                h2d_Requdens_rs       -> Fill (mRequdens, mRs);
                h2d_Requdens_Mtid     -> Fill (mRequdens, mMtidal);
                h2d_Requdens_Rtid     -> Fill (mRequdens, mRtidal);
                h2d_Requdens_Mdelta   -> Fill (mRequdens, mMdelta);
                h2d_Requdens_Mequdens -> Fill (mRequdens, mMequdens);
                h2d_Requdens_J        -> Fill (mRequdens, mJfactor);
                h2d_Requdens_JJcont   -> Fill (mRequdens, mJJcont);


                mTree -> Fill();
        }

      counter++;

  }



    file_output->Write();
    file_output->Close();
    printf("Output file was closed \n");

    //return 0;
}

