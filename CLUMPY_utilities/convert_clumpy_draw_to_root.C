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

#define PI              3.14159265359


char *PLOT_PATH = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/plots";

char *histoMain_Xaxis_title = (char*)"redshift, z";
char *histoMain_Yaxis_title = (char*)"H(z) [km sec^{-1} Mpc^{-1}]";

Float_t Latex1ScaleFactor(0.032);
Float_t Latex2ScaleFactor(0.025);
const Int_t rebinFactor = 1;
const Int_t YesLogScale = 0;
const Int_t XAxis_SetNdivisions   = 10+100*12;
const Int_t YAxis_SetNdivisions   = 10+100*6;

typedef struct {float l,beta, dist;} GALAC_COORD;

void  draw_fits(const char *file, 
                bool Xaxis_LOGSCALE, 
                bool Yaxis_LOGSCALE, 
                Float_t xMin_l,Float_t xMax_l,
                Float_t yMin_l,Float_t yMax_l,
                char *fileNameSuff);


double rho_DM(float DMdensity, float rS, float alpha, float beta, float gamma, float radius, bool convert_flag){

    double GeVcm3 = 37.96*(DMdensity/1e9);// convert Msun/kpc3 to GeV/cm3
    double rho_DM = 0;
    if (convert_flag)
      //rho_DM = GeVcm3/( std::pow(radius/rS, gamma) * std::pow( 1 + std::pow(radius/rS, alpha),  (beta-gamma)/alpha)  );
      rho_DM = GeVcm3 * std::exp(-(2/alpha) * (std::pow(radius/rS, alpha) - 1) );
    else
      //rho_DM = DMdensity/( std::pow(radius/rS, gamma) * std::pow( 1 + std::pow(radius/rS, alpha),  (beta-gamma)/alpha)  );
      rho_DM = DMdensity * std::exp(-(2/alpha) * (std::pow(radius/rS, alpha) - 1) );


    return rho_DM;
}

double accreted_mass(float DMdensity, float sigmaDM_over_sigmaCrit, float time, bool convert_flag){

    //time - formation time from 13.7 Gyr to 0.49 Gyr
    //sigmaDM_over_sigmaCrit -
    //                        sigma_dm > 1e-45 cm-2
    //                        sigma_crit ~= 6x1e-46 cm-2
    //1 GeV = 1.79 x 10^-27 kg * c2
    //1 kg = 0.558659 x 10^27 GeV/c2
    //1 Msun = 1.989 x 10^30 kg
    //1 Msun = 1.111172751 x 10^57 GeV/c2

    double GeV2Msun = (1/1.111172751)*1e-57;
    double f = 0.45*sigmaDM_over_sigmaCrit;
    double GeVcm3 = 0;
    if (convert_flag)
       GeVcm3 = 37.96*(DMdensity/1e9);// convert Msun/kpc3 to GeV/cm3
    else
        GeVcm3 = DMdensity;

    double Macc = 1.3e43 * (GeVcm3/0.3) * time * f * GeV2Msun;
    return Macc;
}


double  accreted_rate(float DMdensity, float sigmaDM_over_sigmaCrit, float DM_particle_mass, float DMmass, float DMradius, bool convert_flag){
    /*
    DM_particle_mass in GeV
    sigmaDM_over_sigmaCrit -
                            sigma_dm > 1e-45 cm-2
                            sigma_crit ~= 6x1e-46 cm-2
    1 GeV = 1.111172751*10e-57 Msun
    */
    double GeV2Msun= (1/1.111172751)*1e-57;
    double velocity = 10;
    double f = 0.45*sigmaDM_over_sigmaCrit;
    double GeVcm3 = 0;
    if (convert_flag)
       GeVcm3 = 37.96*(DMdensity/1e9);// convert Msun/kpc3 to GeV/cm3
    else
        GeVcm3 = DMdensity;

    //double rate = 1.1e27 * (GeVcm3/0.3) * (220/velocity) * (1e3/(DM_particle_mass*GeV2Msun)) * DMmass * DMradius * f;
    double rate = 1.1e27 * (GeVcm3/0.3) * (220/velocity) * (1e3/(DM_particle_mass)) * DMmass * DMradius * f;

    return rate;
}

int FLAG=1;
//int FLAG=2;//for rs outputanalysis


void convert_clumpy_draw_to_root(){
  
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed666_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed1234_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed2241_maxNSubclumps50_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed2233_minSubClumpM_m07_maxSubClumpM_01_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed2131_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed2131_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed2131_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed3333_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed3333_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed1332_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed1332_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed1335_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed1335_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed5550_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed5550_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed5440_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed5440_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed1010_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed1010_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed1088_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed1088_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed5188_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed5188_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed7788_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed7788_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed9111_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed9111_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed7061_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed7061_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed9083_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_massive_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed9083_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_massive_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4444_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4444_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed3120_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed3120_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4446_minSubClumpM_m08_MmaxFrac_0.2_minNClumps10_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4446_minSubClumpM_m08_MmaxFrac_0.2_minNClumps10_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4448_minSubClumpM_m03_MmaxFrac_0.2_minNClumps10_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4448_minSubClumpM_m03_MmaxFrac_0.2_minNClumps10_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4450_minSubClumpM_m04_MmaxFrac_0.6_minNClumps10_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4450_minSubClumpM_m04_MmaxFrac_0.6_minNClumps10_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4452_minSubClumpM_m04_MmaxFrac_0.5_minNClumps10_mDM100GeV_nClumps400_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4452_minSubClumpM_m04_MmaxFrac_0.5_minNClumps10_mDM500GeV_nClumps400_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",

    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4456_minSubClumpM_m07_MmaxFrac_0.6_minNClumps11_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4457_minSubClumpM_m08_MmaxFrac_0.4_minNClumps11_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4457_minSubClumpM_m08_MmaxFrac_0.4_minNClumps11_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4459_minSubClumpM_m08_MmaxFrac_0.4_minNClumps11_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4459_minSubClumpM_m08_MmaxFrac_0.4_minNClumps11_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4462_minSubClumpM_m08_MmaxFrac_0.45_minNClumps11_mDM500GeV_nClumps400_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    //draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4465_minSubClumpM_m07_MmaxFrac_0.35_minNClumps11_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
    draw_fits("/mnt/dampe_data/CLUMPY/output/annihil_seed4465_minSubClumpM_m07_MmaxFrac_0.35_minNClumps11_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.drawn",
              false, false,
              0, 2.5,
              0, 300,
              //(char*)"clumpy_annihil_gal2D_LOS180_0_FOVdiameter360.0deg_nside1024");
              (char*)"clumpy_annihil_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048");
              //(char*)"clumpy_annihil_rs01_seed666_gamma052D_FOV8deg_nside512");
}

void draw_fits(const char *file, 
               bool Xaxis_LOGSCALE, 
               bool Yaxis_LOGSCALE, 
               Float_t xMin_l,Float_t xMax_l,
               Float_t yMin_l,Float_t yMax_l,
               char *fileNameSuff)
{
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed2245_maxNSubclumps150_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed1234_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed666_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seedBaobab_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed2241_maxNSubclumps50_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed2233_minSubClumpM_m07_maxSubClumpM01_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed2131_minSubClumpM_m07_maxSubClumpM01_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed2131_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed3333_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed3333_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed1332_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed1332_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed1335_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed1335_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed5550_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed5440_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed5440_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed1010_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed1010_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed5550_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048__.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_ZHAO_rs01_seed666_gamma052D_FOV8deg_nside512.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed1088_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed1088_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed5188_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed5188_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed7788_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed7788_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed9111_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed9111_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed7061_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed7061_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed9083_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_massive_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed9083_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM500GeV_massive_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4444_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4444_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed3120_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed3120_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4446_minSubClumpM_m08_MmaxFrac_02_minNClumps10_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4446_minSubClumpM_m08_MmaxFrac_02_minNClumps10_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4448_minSubClumpM_m03_MmaxFrac_0.2_minNClumps10_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4448_minSubClumpM_m03_MmaxFrac_0.2_minNClumps10_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4450_minSubClumpM_m04_MmaxFrac_0.6_minNClumps10_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4450_minSubClumpM_m04_MmaxFrac_0.6_minNClumps10_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4452_minSubClumpM_m04_MmaxFrac_0.5_minNClumps10_mDM100GeV_nClumps400_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4452_minSubClumpM_m04_MmaxFrac_0.5_minNClumps10_mDM500GeV_nClumps400_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";

    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4456_minSubClumpM_m07_MmaxFrac_0.6_minNClumps11_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4457_minSubClumpM_m08_MmaxFrac_0.4_minNClumps11_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4457_minSubClumpM_m08_MmaxFrac_0.4_minNClumps11_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4459_minSubClumpM_m08_MmaxFrac_0.4_minNClumps11_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4459_minSubClumpM_m08_MmaxFrac_0.4_minNClumps11_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";

    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4462_minSubClumpM_m08_MmaxFrac_0.45_minNClumps11_mDM500GeV_nClumps400_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    //const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4465_minSubClumpM_m07_MmaxFrac_0.35_minNClumps11_mDM500GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";
    const char *output_fname = (char*)"/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/ntuple_annihil_seed4465_minSubClumpM_m07_MmaxFrac_0.35_minNClumps11_mDM100GeV_nClumps200_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root";

    TFile *file_output = new TFile (output_fname,"RECREATE");

   float mGlong=0;
   float mGlat=0;
   float mDist=0;
   float mDgal=0;

   int mDMclumpID=0;

   float mZ=0;
   float mRdelta=0;
   float mRhoDM_v0=0;
   float mRhoDM_v1=0;
   float mRhoDM_v00=0;
   float mRhoDM_v11=0;
   float mRhoDM=0;
   float mRhos=0;
   float mRhos2=0;
   float mRs=0;
   float m_alphaS=0;
   float m_log10_alphaS=0;


   float mParam1=0;
   float mParam2=0;
   float mParam3=0;
   float mJfactor=0;
   float mJJcont=0;
   float m_log10Jfactor=0;
   float m_log10JJcont=0;

   float mMdelta=0;
   float mMtidal=0;
   float mRtidal=0;
   float mMequdens=0;
   float mRequdens=0;

   float mRho_dm_r00001=0;
   float mRho_dm_r0001=0;
   float mRho_dm_r001=0;
   float mRho_dm_r01=0;
   float mRho_dm_r0=0;
   float mRho_dm_r1=0;
   float mRho_dm_r10=0;
   float mRho_dm_r100=0;

   float mRho_dm__r00001=0;
   float mRho_dm__r0001=0;
   float mRho_dm__r001=0;
   float mRho_dm__r01=0;
   float mRho_dm__r0=0;
   float mRho_dm__r1=0;
   float mRho_dm__r10=0;
   float mRho_dm__r100=0;

   float mMacc_r00001=0;
   float mMacc_r0001=0;
   float mMacc_r001=0;
   float mMacc_r01=0;
   float mMacc_r0=0;
   float mMacc_r1=0;
   float mMacc_r10=0;
   float mMacc_r100=0;

   float mCc_r00001=0;
   float mCc_r0001=0;
   float mCc_r001=0;
   float mCc_r01=0;
   float mCc_r0=0;
   float mCc_r1=0;
   float mCc_r10=0;
   float mCc_r100=0;

   float mRtid_rS=0;
   float mReq_rS=0;
   float mRdelta_rS=0;
   float mDist_rS=0;
   float mDgal_rS=0;

   float mMdelta_tidal=0;
   float mMdelta_eqdens=0;

   float mMtidalMhalo=0;
   float mMeqdensMhalo=0;


   float mRvirial_rS=0;
   float mDELTA=0;
   float mCdelta=0;

  static GALAC_COORD clump_location;

  auto mTree = new TTree("CLUMPY_output", "CLUMPY simulation results");
  mTree->Branch("Clump",                      &clump_location,  "l/F:beta/F:dist/F");
  mTree->Branch("DMclumpID",                  &mDMclumpID,      "mDMclumpID/I");
  mTree->Branch("GalacLongitude_deg",         &mGlong,          "mGlong/F");
  mTree->Branch("GalacLatitude_deg",          &mGlat,           "mGlat/F");
  mTree->Branch("Distance_kpc",               &mDist,           "mDist/F");
  mTree->Branch("DM_Dgal_kpc",                &mDgal,           "mDgal/F");

  mTree->Branch("z",                          &mZ,              "mZ/F");
  mTree->Branch("alpha_S_deg",                &m_alphaS,        "m_alphaS/F");
  mTree->Branch("log10_alpha_S_deg",          &m_log10_alphaS,  "m_log10_alphaS/F");

  mTree->Branch("DMprof_rhoDM_coarse_Msolkpc3", &mRhoDM_v0,       "mRhoDM_v0/F");
  mTree->Branch("DMprof_rhoDM_coarse_GeVcm3",   &mRhoDM_v1,       "mRhoDM_v1/F");
  mTree->Branch("DMprof_rhoDM_smooth_Msolkpc3", &mRhoDM_v00,      "mRhoDM_v00/F");
  mTree->Branch("DMprof_rhoDM_smooth_GeVcm3",   &mRhoDM_v11,      "mRhoDM_v11/F");
  //mTree->Branch("DMprof_rhoDM_Msolkpc3",      &mRhoDM,          "mRhoDM/F");
  mTree->Branch("DMprof_rhos_Msolkpc3",       &mRhos,           "mRhos/F");
  mTree->Branch("DMprof_rhos2_Msolkpc3_2",    &mRhos2,          "mRhos2/F");
  mTree->Branch("DMprof_ScaleRadius_rs_kpc",  &mRs,             "mRs/F");
  mTree->Branch("DMprof_p1",                  &mParam1,         "mParam1/F");
  mTree->Branch("DMprof_p2",                  &mParam2,         "mParam2/F");
  mTree->Branch("DMprof_p3",                  &mParam3,         "mParam3/F");

  mTree->Branch("DM_OverdensityFactor",       &mDELTA,          "mDELTA/F");
  mTree->Branch("DM_C_Rvirial_over_rS",       &mRvirial_rS,     "mRvirial_rS/F");
  mTree->Branch("DM_Rtidal_over_rS",          &mRtid_rS,        "mRtid_rS/F");
  mTree->Branch("DM_Requdens_over_rS",        &mReq_rS,         "mReq_rS/F");
  mTree->Branch("DM_Rdelta_over_rS",          &mRdelta_rS,      "mRdelta_rS/F");
  mTree->Branch("DM_dist_over_rS",            &mDist_rS,        "mDist_rS/F");
  mTree->Branch("DM_Dgal_over_rS",            &mDgal_rS,        "mDgal_rS/F");

  mTree->Branch("DMprof_rhoDM_r00001_GeVcm3",   &mRho_dm_r00001,        "mRho_dm_r00001/F");
  mTree->Branch("DMprof_rhoDM_r0001_GeVcm3",    &mRho_dm_r0001,         "mRho_dm_r0001/F");
  mTree->Branch("DMprof_rhoDM_r001_GeVcm3",     &mRho_dm_r001,          "mRho_dm_r001/F");
  mTree->Branch("DMprof_rhoDM_r01_GeVcm3",      &mRho_dm_r01,           "mRho_dm_r01/F");
  mTree->Branch("DMprof_rhoDM_r0_GeVcm3",       &mRho_dm_r0,            "mRho_dm_r0/F");
  mTree->Branch("DMprof_rhoDM_r1_GeVcm3",       &mRho_dm_r1,            "mRho_dm_r1/F");
  mTree->Branch("DMprof_rhoDM_r10_GeVcm3",      &mRho_dm_r10,           "mRho_dm_r10/F");
  mTree->Branch("DMprof_rhoDM_r100_GeVcm3",     &mRho_dm_r100,          "mRho_dm_r100/F");

  mTree->Branch("DMprof_rhoDM_r00001_Msolkpc3",   &mRho_dm__r00001,        "mRho_dm__r00001/F");
  mTree->Branch("DMprof_rhoDM_r0001_Msolkpc3",    &mRho_dm__r0001,         "mRho_dm__r0001/F");
  mTree->Branch("DMprof_rhoDM_r001_Msolkpc3",     &mRho_dm__r001,          "mRho_dm__r001/F");
  mTree->Branch("DMprof_rhoDM_r01_Msolkpc3",      &mRho_dm__r01,           "mRho_dm__r01/F");
  mTree->Branch("DMprof_rhoDM_r0_Msolkpc3",       &mRho_dm__r0,            "mRho_dm__r0/F");
  mTree->Branch("DMprof_rhoDM_r1_Msolkpc3",       &mRho_dm__r1,            "mRho_dm__r1/F");
  mTree->Branch("DMprof_rhoDM_r10_Msolkpc3",      &mRho_dm__r10,           "mRho_dm__r10/F");
  mTree->Branch("DMprof_rhoDM_r100_Msolkpc3",     &mRho_dm__r100,          "mRho_dm__r100/F");

  mTree->Branch("DM_Maccreted_r00001_Msun",   &mMacc_r00001,        "mMacc_r00001/F");
  mTree->Branch("DM_Maccreted_r0001_Msun",    &mMacc_r0001,         "mMacc_r0001/F");
  mTree->Branch("DM_Maccreted_r001_Msun",     &mMacc_r001,          "mMacc_r001/F");
  mTree->Branch("DM_Maccreted_r01_Msun",      &mMacc_r01,           "mMacc_r01/F");
  mTree->Branch("DM_Maccreted_r0_Msun",       &mMacc_r0,            "mMacc_r0/F");
  mTree->Branch("DM_Maccreted_r1_Msun",       &mMacc_r1,            "mMacc_r1/F");
  mTree->Branch("DM_Maccreted_r10_Msun",      &mMacc_r10,           "mMacc_r10/F");
  mTree->Branch("DM_Maccreted_r100_Msun",     &mMacc_r100,          "mMacc_r100/F");

  mTree->Branch("DM_AccRate_r00001_Hz",   &mCc_r00001,        "mCc_r00001/F");
  mTree->Branch("DM_AccRate_r0001_Hz",    &mCc_r0001,         "mCc_r0001/F");
  mTree->Branch("DM_AccRate_r001_Hz",     &mCc_r001,          "mCc_r001/F");
  mTree->Branch("DM_AccRate_r01_Hz",      &mCc_r01,           "mCc_r01/F");
  mTree->Branch("DM_AccRate_r0_Hz",       &mCc_r0,            "mCc_r0/F");
  mTree->Branch("DM_AccRate_r1_Hz",       &mCc_r1,            "mCc_r1/F");
  mTree->Branch("DM_AccRate_r10_Hz",      &mCc_r10,           "mCc_r10/F");
  mTree->Branch("DM_AccRate_r100_Hz",     &mCc_r100,          "mCc_r100/F");

  mTree->Branch("DM_concentration",         &mCdelta,        "mCdelta/F");

  mTree->Branch("DM_Jfactor_GeV2cm5",         &mJfactor,        "mJfactor/F");
  mTree->Branch("DM_JJgalbkg_GeV2cm5",        &mJJcont,         "mJJcont/F");
  mTree->Branch("DM_log10_Jfactor_GeV2cm5",   &m_log10Jfactor,  "m_log10Jfactor/F");
  mTree->Branch("DM_log10_JJgalbkg_GeV2cm5",  &m_log10JJcont,   "m_log10JJcont/F");

  mTree->Branch("DM_Mdelta_Msol",             &mMdelta,         "mMdelta/F");
  mTree->Branch("Rdelta_kpc",                 &mRdelta,         "mRdelta/F");
  mTree->Branch("DM_Mtidal_Msol",             &mMtidal,         "mMtidal/F");
  mTree->Branch("DM_Rtidal_kpc",              &mRtidal,         "mRtidal/F");
  mTree->Branch("DM_Mequdens_Msol",           &mMequdens,       "mMequdens/F");
  mTree->Branch("DM_Requdens_kpc",            &mRequdens,       "mRequdens/F");
  mTree->Branch("DM_Mdelta_Mtidal",           &mMdelta_tidal,   "mMdelta_tidal/F");
  mTree->Branch("DM_Mdelta_Meqdens",          &mMdelta_eqdens,  "mMdelta_eqdens/F");
  mTree->Branch("DM_Mtidal_Mhalo",            &mMtidalMhalo,    "mMtidalMhalo/F");
  mTree->Branch("DM_Meqdens_Mhalo",           &mMeqdensMhalo,   "mMeqdensMhalo/F");


  TProfile *prf_r_Macc               = new TProfile ("prf_r_Macc",           ";radius [kpc];#LT M_{acc} #GT [M_{sun}]",             10000, 1e-4, 1e2, 1e-10, 1e10);
  TProfile *prf_r00001_rhoDM_Macc    = new TProfile ("prf_r00001_rhoDM_Macc", ";#rho_{DM} [GeV/cm3];#LT M_{acc} #GT [M_{sun}]",     10000, 1e-4, 1e9, 1e-10, 1e10);
  TProfile *prf_r0001_rhoDM_Macc     = new TProfile ("prf_r0001_rhoDM_Macc",  ";#rho_{DM} [GeV/cm3];#LT M_{acc} #GT [M_{sun}]",     10000, 1e-4, 1e9, 1e-10, 1e10);
  TProfile *prf_r001_rhoDM_Macc      = new TProfile ("prf_r001_rhoDM_Macc",   ";#rho_{DM} [GeV/cm3];#LT M_{acc} #GT [M_{sun}]",     10000, 1e-4, 1e9, 1e-10, 1e10);
  TProfile *prf_r01_rhoDM_Macc       = new TProfile ("prf_r01_rhoDM_Macc",    ";#rho_{DM} [GeV/cm3];#LT M_{acc} #GT [M_{sun}]",     10000, 1e-4, 1e9, 1e-10, 1e10);
  TProfile *prf_r0_rhoDM_Macc        = new TProfile ("prf_r0_rhoDM_Macc",     ";#rho_{DM} [GeV/cm3];#LT M_{acc} #GT [M_{sun}]",     10000, 1e-4, 1e9, 1e-10, 1e10);
  TProfile *prf_r1_rhoDM_Macc        = new TProfile ("prf_r1_rhoDM_Macc",     ";#rho_{DM} [GeV/cm3];#LT M_{acc} #GT [M_{sun}]",     10000, 1e-4, 1e9, 1e-10, 1e10);
  TProfile *prf_r10_rhoDM_Macc       = new TProfile ("prf_r10_rhoDM_Macc",    ";#rho_{DM} [GeV/cm3];#LT M_{acc} #GT [M_{sun}]",     10000, 1e-4, 1e9, 1e-10, 1e10);
  TProfile *prf_r100_rhoDM_Macc      = new TProfile ("prf_r100_rhoDM_Macc",   ";#rho_{DM} [GeV/cm3];#LT M_{acc} #GT [M_{sun}]",     10000, 1e-4, 1e9, 1e-10, 1e10);



  TProfile *prf_rhos_Mdelta     = new TProfile ("rhos_Mdelta",     ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT M_{#Delta} #GT [M_{sun}]",     10000, 0, 1.5e9, 0, 2000);
  TProfile *prf_rhos_Mtid       = new TProfile ("rhos_Mtid",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT M_{tidal} #GT [M_{sun}]",      10000, 0, 1.5e9, 0, 2000);
  TProfile *prf_rhos_Rtid       = new TProfile ("rhos_Rtid",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT R_{tidal} #GT [kpc]",          10000, 0, 1.5e9, 0, 1);
  TProfile *prf_rhos_Mequdens   = new TProfile ("rhos_Mequdens",   ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT M_{equdens} #GT [M_{sun}]",    10000, 0, 1.5e9, 0, 2000);
  TProfile *prf_rhos_Requdens   = new TProfile ("rhos_Requdens",   ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT R_{equdens} #GT [kpc]",        10000, 0, 1.5e9, 0, 1);
  TProfile *prf_rhos_J          = new TProfile ("rhos_J",          ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        10000, 0, 1.5e9, 0, 1e16);
  TProfile *prf_rhos_JJcont     = new TProfile ("rhos_JJcont",     ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 10000, 0, 1.5e9, 0, 100);
  TProfile *prf_rhos_dist       = new TProfile ("rhos_dist",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];#LT dist #GT [kpc]",                       10000, 0, 1.5e9, 0, 40);

  TProfile *prf_rs_Mdelta     = new TProfile ("rS_Mdelta",     ";scaled radius r_{s} [kpc];#LT M_{#Delta} #GT [M_{sun}]",         40000, 0, 4e-3, 0, 2000);
  TProfile *prf_rs_Mtid       = new TProfile ("rS_Mtid",       ";scaled radius r_{s} [kpc];#LT M_{tidal} #GT [M_{sun}]",          40000, 0, 4e-3, 0, 2000);
  TProfile *prf_rs_Rtid       = new TProfile ("rS_Rtid",       ";scaled radius r_{s} [kpc];#LT R_{tidal} #GT [kpc]",              40000, 0, 4e-3, 0, 1);
  TProfile *prf_rs_Mequdens   = new TProfile ("rS_Mequdens",   ";scaled radius r_{s} [kpc];#LT M_{equdens} #GT [M_{sun}]",        40000, 0, 4e-3, 0, 2000);
  TProfile *prf_rs_Requdens   = new TProfile ("rS_Requdens",   ";scaled radius r_{s} [kpc];#LT R_{equdens} #GT [kpc]",            40000, 0, 4e-3, 0, 1);
  TProfile *prf_rs_J          = new TProfile ("rS_J",          ";scaled radius r_{s} [kpc];#LT J-factor #GT [GeV^{2}cm^{-5}]",        40000, 0, 4e-3, 0, 1e16);
  TProfile *prf_rs_JJcont     = new TProfile ("rS_JJcont",     ";scaled radius r_{s} [kpc];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 40000, 0, 4e-3, 0, 100);
  TProfile *prf_rs_rhos       = new TProfile ("rS_rhos",       ";scaled radius r_{s} [kpc];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       40000, 0, 4e-3, 0, 2e9);
  TProfile *prf_rs_rhos2      = new TProfile ("rS_rhos2",      ";scaled radius r_{s} [kpc];#LT #rho_{s}^{2} #GT [(M_{sun}/kpc^{3})^{2}]",   40000, 0, 4e-3, 0, 1e18);
  TProfile *prf_rs_dist       = new TProfile ("rS_dist",       ";scaled radius r_{s} [kpc];#LT dist #GT [kpc]",                             40000, 0, 4e-3, 0, 40);

  TProfile *prf_rrs_dist      = new TProfile ("prf_rrs_dist",   ";r/r_{s};#LT distance #GT [kpc]",                             40000, 0, 100, 0, 400);
  TProfile *prf_rrs_dgal      = new TProfile ("prf_rrs_dgal",   ";r/r_{s};#LT distance.gal #GT [kpc]",                         40000, 0, 100, 0, 400);


  TProfile *prf_alphas_Mdelta     = new TProfile ("aS_Mdelta",     ";#alpha_{s}=asin(r_{s}/d) [deg];#LT M_{#Delta} #GT [M_{sun}]",                   10000, 0, 1e-1, 0, 2000);
  TProfile *prf_alphas_Mtid       = new TProfile ("aS_Mtid",       ";#alpha_{s}=asin(r_{s}/d) [deg];#LT M_{tidal} #GT [M_{sun}]",                    10000, 0, 1e-1, 0, 2000);
  TProfile *prf_alphas_Rtid       = new TProfile ("aS_Rtid",       ";#alpha_{s}=asin(r_{s}/d) [deg];#LT R_{tidal} #GT [kpc]",                        10000, 0, 1e-1, 0, 1);
  TProfile *prf_alphas_Mequdens   = new TProfile ("aS_Mequdens",   ";#alpha_{s}=asin(r_{s}/d) [deg];#LT M_{equdens} #GT [M_{sun}]",                  10000, 0, 1e-1, 0, 2000);
  TProfile *prf_alphas_Requdens   = new TProfile ("aS_Requdens",   ";#alpha_{s}=asin(r_{s}/d) [deg];#LT R_{equdens} #GT [kpc]",                      10000, 0, 1e-1, 0, 1);
  TProfile *prf_alphas_J          = new TProfile ("aS_J",          ";#alpha_{s}=asin(r_{s}/d) [deg];#LT J-factor #GT [GeV^{2}cm^{-5}]",              10000, 0, 1e-1, 0, 1e16);
  TProfile *prf_alphas_JJcont     = new TProfile ("aS_JJcont",     ";#alpha_{s}=asin(r_{s}/d) [deg];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]",       10000, 0, 1e-1, 0, 100);
  TProfile *prf_alphas_rhos       = new TProfile ("aS_rhos",       ";#alpha_{s}=asin(r_{s}/d) [deg];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",             10000, 0, 1e-1, 0, 2e9);
  TProfile *prf_alphas_rhos2      = new TProfile ("aS_rhos2",      ";#alpha_{s}=asin(r_{s}/d) [deg];#LT #rho_{s}^{2} #GT [(M_{sun}/kpc^{3})^{2}]",   10000, 0, 1e-1, 0, 1e18);
  TProfile *prf_alphas_dist       = new TProfile ("aS_dist",       ";#alpha_{s}=asin(r_{s}/d) [deg];#LT dist #GT [kpc]",                             10000, 0, 1e-1, 0, 40);



  TProfile *prf_Glong_Mdelta     = new TProfile ("Glong_Mdelta",    ";Galactic longitude [deg];#LT M_{#Delta} #GT [M_{sun}]",     10000, -200, 200, 0, 2000);
  TProfile *prf_Glong_Mtid       = new TProfile ("Glong_Mtid",      ";Galactic longitude [deg];#LT M_{tidal} #GT [M_{sun}]",      10000, -200, 200, 0, 2000);
  TProfile *prf_Glong_Rtid       = new TProfile ("Glong_Rtid",      ";Galactic longitude [deg];#LT R_{tidal} #GT [kpc]",          10000, -200, 200, 0, 1);
  TProfile *prf_Glong_Mequdens   = new TProfile ("Glong_Mequdens",  ";Galactic longitude [deg];#LT M_{equdens} #GT [M_{sun}]",    10000, -200, 200, 0, 2000);
  TProfile *prf_Glong_Requdens   = new TProfile ("Glong_Requdens",  ";Galactic longitude [deg];#LT R_{equdens} #GT [kpc]",        10000, -200, 200, 0, 1);
  TProfile *prf_Glong_J          = new TProfile ("Glong_J",         ";Galactic longitude [deg];#LT J-factor #GT [GeV^{2}cm^{-5}]",        10000, -200, 200, 0, 1e16);
  TProfile *prf_Glong_JJcont     = new TProfile ("Glong_JJcont",    ";Galactic longitude [deg];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 10000, -200, 200, 0, 100);

  TProfile *prf_Glat_Mdelta     = new TProfile ("Glat_Mdelta",    ";Galactic latitude [deg];#LT M_{#Delta} #GT [M_{sun}]",     10000, -100, 100, 0, 2000);
  TProfile *prf_Glat_Mtid       = new TProfile ("Glat_Mtid",      ";Galactic latitude [deg];#LT M_{tidal} #GT [M_{sun}]",      10000, -100, 100, 0, 2000);
  TProfile *prf_Glat_Rtid       = new TProfile ("Glat_Rtid",      ";Galactic latitude [deg];#LT R_{tidal} #GT [kpc]",          10000, -100, 100, 0, 1);
  TProfile *prf_Glat_Mequdens   = new TProfile ("Glat_Mequdens",  ";Galactic latitude [deg];#LT M_{equdens} #GT [M_{sun}]",    10000, -100, 100, 0, 2000);
  TProfile *prf_Glat_Requdens   = new TProfile ("Glat_Requdens",  ";Galactic latitude [deg];#LT R_{equdens} #GT [kpc]",        10000, -100, 100, 0, 1);
  TProfile *prf_Glat_J          = new TProfile ("Glat_J",         ";Galactic latitude [deg];#LT J-factor #GT [GeV^{2}cm^{-5}]",        10000, -100, 100, 0, 1e16);
  TProfile *prf_Glat_JJcont     = new TProfile ("Glat_JJcont",    ";Galactic latitude [deg];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 10000, -100, 100, 0, 100);

  TProfile *prf_dist_Mdelta     = new TProfile ("dist_Mdelta",    ";Distance [kpc];#LT M_{#Delta} #GT [M_{sun}]",             20000, 0, 250, 0, 10000);
  TProfile *prf_dist_Mtid       = new TProfile ("dist_Mtid",      ";Distance [kpc];#LT M_{tidal} #GT [M_{sun}]",              20000, 0, 250, 0, 10000);
  TProfile *prf_dist_Mequdens   = new TProfile ("dist_Mequdens",  ";Distance [kpc];#LT M_{equdens} #GT [M_{sun}]",            20000, 0, 250, 0, 10000);
  TProfile *prf_dist_Rdelta     = new TProfile ("dist_Rdelta",    ";Distance [kpc];outbound of DM halo, #LT R_{#Delta} #GT [kpc]", 20000, 0, 250, 0, 200);
  TProfile *prf_dist_Rtid       = new TProfile ("dist_Rtid",      ";Distance [kpc];#LT R_{tidal} #GT [kpc]",                  20000, 0, 250, 0, 1);
  TProfile *prf_dist_Requdens   = new TProfile ("dist_Requdens",  ";Distance [kpc];#LT R_{equdens} #GT [kpc]",                20000, 0, 250, 0, 1);
  TProfile *prf_dist_J          = new TProfile ("dist_J",         ";Distance [kpc];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 250, 0, 1e16);
  TProfile *prf_dist_JJcont     = new TProfile ("dist_JJcont",    ";Distance [kpc];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 250, 0, 100);
  TProfile *prf_dist_rrs        = new TProfile ("dist_rrs",       ";Distance [kpc];#LT r/r_{S} #GT ",                         20000, 0, 250, 0, 1e10);
  TProfile *prf_dist_rs         = new TProfile ("dist_rs",        ";Distance [kpc];#LT r_{S} #GT [kpc]",                      20000, 0, 250, 0, 10);
  TProfile *prf_dist_rhos       = new TProfile ("dist_rhos",      ";Distance [kpc];#LT #rho_{S} #GT [M_{sun}/kpc^{3}]",       20000, 0, 250, 0, 1e10);
  TProfile *prf_dist_rhoDM_v0   = new TProfile ("dist_rhoDM_v0",   ";coarsed Distance [kpc];#LT #rho_{DM} #GT [M_{sun}/kpc^{3}]",  20000, 0, 250, 0, 1e10);
  TProfile *prf_dist_rhoDM_v1   = new TProfile ("dist_rhoDM_v1",   ";coarsed Distance [kpc];#LT #rho_{DM} #GT [GeV/cm^{3}]",       20000, 0, 250, 0, 1e10);
  TProfile *prf_dist_rhoDM_v00  = new TProfile ("dist_rhoDM_v00",  ";smoothed Distance [kpc];#LT #rho_{DM} #GT [M_{sun}/kpc^{3}]", 20000, 0, 250, 0, 1e10);
  TProfile *prf_dist_rhoDM_v11  = new TProfile ("dist_rhoDM_v11",  ";smoothed Distance [kpc];#LT #rho_{DM} #GT [GeV/cm^{3}]",      20000, 0, 250, 0, 1e10);


  TProfile *prf_dgal_Mdelta     = new TProfile ("dgal_Mdelta",    ";DistanceGal [kpc];#LT M_{#Delta} #GT [M_{sun}]",             20000, 0, 250, 0, 10000);
  TProfile *prf_dgal_Mtid       = new TProfile ("dgal_Mtid",      ";DistanceGal [kpc];#LT M_{tidal} #GT [M_{sun}]",              20000, 0, 250, 0, 10000);
  TProfile *prf_dgal_Mequdens   = new TProfile ("dgal_Mequdens",  ";DistanceGal [kpc];#LT M_{equdens} #GT [M_{sun}]",            20000, 0, 250, 0, 10000);
  TProfile *prf_dgal_Rdelta     = new TProfile ("dgal_Rdelta",    ";DistanceGal [kpc];outbound of DM halo, #LT R_{#Delta} #GT [kpc]", 20000, 0, 250, 0, 200);
  TProfile *prf_dgal_Rtid       = new TProfile ("dgal_Rtid",      ";DistanceGal [kpc];#LT R_{tidal} #GT [kpc]",                  20000, 0, 250, 0, 1);
  TProfile *prf_dgal_Requdens   = new TProfile ("dgal_Requdens",  ";DistanceGal [kpc];#LT R_{equdens} #GT [kpc]",                20000, 0, 250, 0, 1);
  TProfile *prf_dgal_J          = new TProfile ("dgal_J",         ";DistanceGal [kpc];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 250, 0, 1e16);
  TProfile *prf_dgal_JJcont     = new TProfile ("dgal_JJcont",    ";DistanceGal [kpc];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 250, 0, 100);
  TProfile *prf_dgal_rrs        = new TProfile ("dgal_rrs",       ";DistanceGal [kpc];#LT r/r_{S} #GT ",                         20000, 0, 250, 0, 1e10);
  TProfile *prf_dgal_rs         = new TProfile ("dgal_rs",        ";DistanceGal [kpc];#LT r_{S} #GT [kpc]",                      20000, 0, 250, 0, 10);
  TProfile *prf_dgal_rhos       = new TProfile ("dgal_rhos",      ";DistanceGal [kpc];#LT #rho_{S} #GT [M_{sun}/kpc^{3}]",       20000, 0, 250, 0, 1e10);
  TProfile *prf_dgal_rhoDM_v0   = new TProfile ("dgal_rhoDM_v0",  ";DistanceGal [kpc];#LT #rho_{DM} #GT [M_{sun}/kpc^{3}]",      20000, 0, 250, 0, 1e10);
  TProfile *prf_dgal_rhoDM_v1   = new TProfile ("dgal_rhoDM_v1",  ";DistanceGal [kpc];#LT #rho_{DM} #GT [GeV/cm^{3}]",           20000, 0, 250, 0, 1e10);
  TProfile *prf_dgal_rhoDM_v00  = new TProfile ("dgal_rhoDM_v00", ";DistanceGal [kpc];#LT #rho_{DM} #GT [M_{sun}/kpc^{3}]",      20000, 0, 250, 0, 1e10);
  TProfile *prf_dgal_rhoDM_v11  = new TProfile ("dgal_rhoDM_v11", ";DistanceGal [kpc];#LT #rho_{DM} #GT [GeV/cm^{3}]",           20000, 0, 250, 0, 1e10);


  TProfile *prf_Mdelta_rhos       = new TProfile ("Mdelta_rhos",       ";M_{#Delta} [M_{sun}];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]", 20000, 0, 200, 0, 2e9);
  TProfile *prf_Mdelta_rs         = new TProfile ("Mdelta_rs",         ";M_{#Delta} [M_{sun}];#LT r_{s} #GT [kpc]",                20000, 0, 200, 0, 1e-3);
  TProfile *prf_Mdelta_Mtid       = new TProfile ("Mdelta_Mtid",       ";M_{#Delta} [M_{sun}];#LT M_{tidal} #GT [M_{sun}]",        20000, 0, 200, 0, 2000);
  TProfile *prf_Mdelta_Rtid       = new TProfile ("Mdelta_Rtid",       ";M_{#Delta} [M_{sun}];#LT R_{tidal} #GT [kpc]",            20000, 0, 200, 0, 1);
  TProfile *prf_Mdelta_Mequdens   = new TProfile ("Mdelta_Mequdens",   ";M_{#Delta} [M_{sun}];#LT M_{equdens} #GT [M_{sun}]",      20000, 0, 200, 0, 2000);
  TProfile *prf_Mdelta_Requdens   = new TProfile ("Mdelta_Requdens",   ";M_{#Delta} [M_{sun}];#LT R_{equdens} #GT [kpc]",          20000, 0, 200, 0, 1);
  TProfile *prf_Mdelta_J          = new TProfile ("Mdelta_J",          ";M_{#Delta} [M_{sun}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 200, 0, 1e16);
  TProfile *prf_Mdelta_JJcont     = new TProfile ("Mdelta_JJcont",     ";M_{#Delta} [M_{sun}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 200, 0, 100);

  TProfile *prf_Mtid_rhos       = new TProfile ("Mtid_rhos",      ";M_{tidal} [M_{sun}];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       20000, 0, 800, 0, 2e9);
  TProfile *prf_Mtid_rs         = new TProfile ("Mtid_rs",        ";M_{tidal} [M_{sun}];#LT r_{s} #GT [kpc]",                      20000, 0, 800, 0, 1e-3);
  TProfile *prf_Mtid_Mdelta     = new TProfile ("Mtid_Mdelta",    ";M_{tidal} [M_{sun}];#LT M_{#Delta} #GT [M_{sun}]",             20000, 0, 800, 0, 2000);
  TProfile *prf_Mtid_Rtid       = new TProfile ("Mtid_Rtid",      ";M_{tidal} [M_{sun}];#LT R_{tidal} #GT [kpc]",                  20000, 0, 800, 0, 1);
  TProfile *prf_Mtid_Mequdens   = new TProfile ("Mtid_Mequdens",  ";M_{tidal} [M_{sun}];#LT M_{equdens} #GT [M_{sun}]",            20000, 0, 800, 0, 2000);
  TProfile *prf_Mtid_Requdens   = new TProfile ("Mtid_Requdens",  ";M_{tidal} [M_{sun}];#LT R_{equdens} #GT [kpc]",                20000, 0, 800, 0, 1);
  TProfile *prf_Mtid_J          = new TProfile ("Mtid_J",         ";M_{tidal} [M_{sun}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 800, 0, 1e16);
  TProfile *prf_Mtid_JJcont     = new TProfile ("Mtid_JJcont",    ";M_{tidal} [M_{sun}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 800, 0, 100);

  TProfile *prf_Rtid_rhos       = new TProfile ("Rtid_rhos",      ";R_{tidal} [kpc];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",        40000, 0.002, 0.05, 0, 2e9);
  TProfile *prf_Rtid_rs         = new TProfile ("Rtid_rs",        ";R_{tidal} [kpc];#LT r_{s} #GT [kpc]",                       40000, 0.002, 0.05, 0, 1e-3);
  TProfile *prf_Rtid_Mdelta     = new TProfile ("Rtid_Mdelta",    ";R_{tidal} [kpc];#LT M_{#Delta} #GT [M_{sun}]",              40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Rtid_Mtid       = new TProfile ("Rtid_Mtid",      ";R_{tidal} [kpc];#LT M_{tidal} #GT [M_{sun}]",               40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Rtid_Mequdens   = new TProfile ("Rtid_Mequdens",  ";R_{tidal} [kpc];#LT M_{equdens} #GT [M_{sun}]",             40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Rtid_Requdens   = new TProfile ("Rtid_Requdens",  ";R_{tidal} [kpc];#LT R_{equdens} #GT [kpc]",                 40000, 0.002, 0.05, 0, 1);
  TProfile *prf_Rtid_J          = new TProfile ("Rtid_J",         ";R_{tidal} [kpc];#LT J-factor #GT [GeV^{2}cm^{-5}]",         40000, 0.002, 0.05, 0, 1e16);
  TProfile *prf_Rtid_JJcont     = new TProfile ("Rtid_JJcont",    ";R_{tidal} [kpc];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]",  40000, 0.002, 0.05, 0, 100);

  TProfile *prf_Mequdens_rhos       = new TProfile ("Mequdens_rhos",      ";M_{equdens} [M_{sun}];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       40000, 0, 500, 0, 2e9);
  TProfile *prf_Mequdens_rs         = new TProfile ("Mequdens_rs",        ";M_{equdens} [M_{sun}];#LT r_{s} #GT [kpc]",                      40000, 0, 500, 0, 1e-3);
  TProfile *prf_Mequdens_Mdelta     = new TProfile ("Mequdens_Mdelta",    ";M_{equdens} [M_{sun}];#LT M_{#Delta} #GT [M_{sun}]",             40000, 0, 500, 0, 2000);
  TProfile *prf_Mequdens_Requdens   = new TProfile ("Mequdens_Requdens",  ";M_{equdens} [M_{sun}];#LT R_{equdens} #GT [kpc]",                40000, 0, 500, 0, 1);
  TProfile *prf_Mequdens_Mtid       = new TProfile ("Mequdens_Mtid",      ";M_{equdens} [M_{sun}];#LT M_{tidal} #GT [M_{sun}]",              40000, 0, 500, 0, 2000);
  TProfile *prf_Mequdens_Rtid       = new TProfile ("Mequdens_Rtid",      ";M_{equdens} [M_{sun}];#LT R_{tidal} #GT [kpc]",                  40000, 0, 500, 0, 2000);
  TProfile *prf_Mequdens_J          = new TProfile ("Mequdens_J",         ";M_{equdens} [M_{sun}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        40000, 0, 500, 0, 1e16);
  TProfile *prf_Mequdens_JJcont     = new TProfile ("Mequdens_JJcont",    ";M_{equdens} [M_{sun}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 40000, 0, 500, 0, 100);

  TProfile *prf_Requdens_rhos       = new TProfile ("Requdens_rhos",      ";R_{equdens} [kpc];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       40000, 0.002, 0.05, 0, 2e9);
  TProfile *prf_Requdens_rs         = new TProfile ("Requdens_rs",        ";R_{equdens} [kpc];#LT r_{s} #GT [kpc]",                      40000, 0.002, 0.05, 0, 1e-3);
  TProfile *prf_Requdens_Mdelta     = new TProfile ("Requdens_Mdelta",    ";R_{equdens} [kpc];#LT M_{#Delta} #GT [M_{sun}]",             40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Requdens_Mequdens   = new TProfile ("Requdens_Mequdens",  ";R_{equdens} [kpc];#LT M_{equdens} #GT [M_{sun}]",            40000, 0.002, 0.05, 0, 2000);
  TProfile *prf_Requdens_Rtid       = new TProfile ("Requdens_Rtid",      ";R_{equdens} [kpc];#LT R_{tidal} #GT [kpc]",                  40000, 0.002, 0.05, 0, 1);
  TProfile *prf_Requdens_Mtid       = new TProfile ("Requdens_Mtid",      ";R_{equdens} [kpc];#LT M_{tidal} #GT [M_{sun}]",              40000, 0.002, 0.05, 0, 1000);
  TProfile *prf_Requdens_J          = new TProfile ("Requdens_J",         ";R_{equdens} [kpc];#LT J-factor #GT [GeV^{2}cm^{-5}]",        40000, 0.002, 0.05, 0, 1e16);
  TProfile *prf_Requdens_JJcont     = new TProfile ("Requdens_JJcont",    ";R_{equdens} [kpc];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 40000, 0.002, 0.05, 0, 100);


  TProfile *prf_MequdensDist_rhos       = new TProfile ("MequdensDist_rhos",      ";M_{equdens}/d^{2} [M_{sun} kpc^{-2}];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       20000, 0, 1000, 0, 2e9);
  TProfile *prf_MequdensDist_rs         = new TProfile ("MequdensDist_rs",        ";M_{equdens}/d^{2} [M_{sun} kpc^{-2}];#LT r_{s} #GT [kpc]",                      20000, 0, 1000, 0, 1e-3);
  TProfile *prf_MequdensDist_Mdelta     = new TProfile ("MequdensDist_Mdelta",    ";M_{equdens}/d^{2} [M_{sun} kpc^{-2}];#LT M_{#Delta} #GT [M_{sun}]",             20000, 0, 1000, 0, 2000);
  TProfile *prf_MequdensDist_Requdens   = new TProfile ("MequdensDist_Requdens",  ";M_{equdens}/d^{2} [M_{sun} kpc^{-2}];#LT R_{equdens} #GT [kpc]",                20000, 0, 1000, 0, 1);
  TProfile *prf_MequdensDist_Mtid       = new TProfile ("MequdensDist_Mtid",      ";M_{equdens}/d^{2} [M_{sun} kpc^{-2}];#LT M_{tidal} #GT [M_{sun}]",              20000, 0, 1000, 0, 2000);
  TProfile *prf_MequdensDist_Rtid       = new TProfile ("MequdensDist_Rtid",      ";M_{equdens}/d^{2} [M_{sun} kpc^{-2}];#LT R_{tidal} #GT [kpc]",                  20000, 0, 1000, 0, 2000);
  TProfile *prf_MequdensDist_J          = new TProfile ("MequdensDist_J",         ";M_{equdens}/d^{2} [M_{sun} kpc^{-2}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 1000, 0, 1e16);
  TProfile *prf_MequdensDist_JJcont     = new TProfile ("MequdensDist_JJcont",    ";M_{equdens}/d^{2} [M_{sun} kpc^{-2}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 1000, 0, 100);

  TProfile *prf_MtidDist_rhos       = new TProfile ("MtidDist_rhos",      ";M_{tidal}/d^{2} [M_{sun} kpc^{-2}];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]",       20000, 0, 1000, 0, 2e9);
  TProfile *prf_MtidDist_rs         = new TProfile ("MtidDist_rs",        ";M_{tidal}/d^{2} [M_{sun} kpc^{-2}];#LT r_{s} #GT [kpc]",                      20000, 0, 1000, 0, 1e-3);
  TProfile *prf_MtidDist_Mdelta     = new TProfile ("MtidDist_Mdelta",    ";M_{tidal}/d^{2} [M_{sun} kpc^{-2}];#LT M_{#Delta} #GT [M_{sun}]",             20000, 0, 1000, 0, 2000);
  TProfile *prf_MtidDist_Rtid       = new TProfile ("MtidDist_Rtid",      ";M_{tidal}/d^{2} [M_{sun} kpc^{-2}];#LT R_{tidal} #GT [kpc]",                  20000, 0, 1000, 0, 1);
  TProfile *prf_MtidDist_Mequdens   = new TProfile ("MtidDist_Mequdens",  ";M_{tidal}/d^{2} [M_{sun} kpc^{-2}];#LT M_{equdens} #GT [M_{sun}]",            20000, 0, 1000, 0, 2000);
  TProfile *prf_MtidDist_Requdens   = new TProfile ("MtidDist_Requdens",  ";M_{tidal}/d^{2} [M_{sun} kpc^{-2}];#LT R_{equdens} #GT [kpc]",                20000, 0, 1000, 0, 1);
  TProfile *prf_MtidDist_J          = new TProfile ("MtidDist_J",         ";M_{tidal}/d^{2} [M_{sun} kpc^{-2}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 1000, 0, 1e16);
  TProfile *prf_MtidDist_JJcont     = new TProfile ("MtidDist_JJcont",    ";M_{tidal}/d^{2} [M_{sun} kpc^{-2}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 1000, 0, 100);


  TProfile *prf_MdeltaDist_rhos       = new TProfile ("MdeltaDist_rhos",       ";M_{#Delta}/d^{2} [M_{sun} kpc^{-2}];#LT #rho_{s} #GT [M_{sun}/kpc^{3}]", 20000, 0, 1000, 0, 2e9);
  TProfile *prf_MdeltaDist_rs         = new TProfile ("MdeltaDist_rs",         ";M_{#Delta}/d^{2} [M_{sun} kpc^{-2}];#LT r_{s} #GT [kpc]",                20000, 0, 1000, 0, 1e-3);
  TProfile *prf_MdeltaDist_Mtid       = new TProfile ("MdeltaDist_Mtid",       ";M_{#Delta}/d^{2} [M_{sun} kpc^{-2}];#LT M_{tidal} #GT [M_{sun}]",        20000, 0, 1000, 0, 2000);
  TProfile *prf_MdeltaDist_Rtid       = new TProfile ("MdeltaDist_Rtid",       ";M_{#Delta}/d^{2} [M_{sun} kpc^{-2}];#LT R_{tidal} #GT [kpc]",            20000, 0, 1000, 0, 1);
  TProfile *prf_MdeltaDist_Mequdens   = new TProfile ("MdeltaDist_Mequdens",   ";M_{#Delta}/d^{2} [M_{sun} kpc^{-2}];#LT M_{equdens} #GT [M_{sun}]",      20000, 0, 1000, 0, 2000);
  TProfile *prf_MdeltaDist_Requdens   = new TProfile ("MdeltaDist_Requdens",   ";M_{#Delta}/d^{2} [M_{sun} kpc^{-2}];#LT R_{equdens} #GT [kpc]",          20000, 0, 1000, 0, 1);
  TProfile *prf_MdeltaDist_J          = new TProfile ("MdeltaDist_J",          ";M_{#Delta}/d^{2} [M_{sun} kpc^{-2}];#LT J-factor #GT [GeV^{2}cm^{-5}]",        20000, 0, 1000, 0, 1e16);
  TProfile *prf_MdeltaDist_JJcont     = new TProfile ("MdeltaDist_JJcont",     ";M_{#Delta}/d^{2} [M_{sun} kpc^{-2}];#LT J/J_{continuum} #GT [GeV^{2}cm^{-5}]", 20000, 0, 1000, 0, 100);





  TH2F *h2d_rhos_Mdelta     = new TH2F ("h2d_rhos_Mdelta",     ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];M_{#Delta} [M_{sun}]",             4000, 0, 1.0e9, 1000, 0, 1000);
  TH2F *h2d_rhos_Mtid       = new TH2F ("h2d_rhos_Mtid",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];M_{tidal} [M_{sun}]",              4000, 0, 1.0e9, 1000, 0, 800);
  TH2F *h2d_rhos_Rtid       = new TH2F ("h2d_rhos_Rtid",       ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];R_{tidal} [kpc]",                  4000, 0, 1.0e9, 1000, 1e-8, 5e-2);
  TH2F *h2d_rhos_Mequdens   = new TH2F ("h2d_rhos_Mequdens",   ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];M_{equdens} [M_{sun}]",            4000, 0, 1.0e9, 1000, 0, 1000);
  TH2F *h2d_rhos_Requdens   = new TH2F ("h2d_rhos_Requdens",   ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];R_{equdens} [kpc]",                4000, 0, 1.0e9, 1000, 1e-8, 5e-2);
  TH2F *h2d_rhos_J          = new TH2F ("h2d_rhos_J",          ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];J-factor [GeV^{2}cm^{-5}]",        4000, 0, 1.0e9, 1000, 1e14, 1e16);
  TH2F *h2d_rhos_JJcont     = new TH2F ("h2d_rhos_JJcont",     ";DM saturation dens. #rho_{s} [M_{sun}/kpc^{3}];J/J_{continuum} [GeV^{2}cm^{-5}]", 4000, 0, 1.0e9, 1000, 0, 20);

  TH2F *h2d_rs_Mdelta     = new TH2F ("h2d_rs_Mdelta",     ";scaled radius r_{s} [kpc];M_{#Delta} [M_{sun}]",             4000, 0, 4e-3, 1000, 0, 1000);
  TH2F *h2d_rs_Mtid       = new TH2F ("h2d_rs_Mtid",       ";scaled radius r_{s} [kpc];M_{tidal} [M_{sun}]",              4000, 0, 4e-3, 1000, 0, 800);
  TH2F *h2d_rs_Rtid       = new TH2F ("h2d_rs_Rtid",       ";scaled radius r_{s} [kpc];R_{tidal} [kpc]",                  4000, 0, 4e-3, 1000, 1e-8, 5e-2);
  TH2F *h2d_rs_Mequdens   = new TH2F ("h2d_rs_Mequdens",   ";scaled radius r_{s} [kpc];M_{equdens} [M_{sun}]",            4000, 0, 4e-3, 1000, 0, 1000);
  TH2F *h2d_rs_Requdens   = new TH2F ("h2d_rs_Requdens",   ";scaled radius r_{s} [kpc];R_{equdens} [kpc]",                4000, 0, 4e-3, 1000, 1e-8, 5e-2);
  TH2F *h2d_rs_J          = new TH2F ("h2d_rs_J",          ";scaled radius r_{s} [kpc];J-factor [GeV^{2}cm^{-5}]",        4000, 0, 4e-3, 1000, 1e14, 1e16);
  TH2F *h2d_rs_JJcont     = new TH2F ("h2d_rs_JJcont",     ";scaled radius r_{s} [kpc];J/J_{continuum} [GeV^{2}cm^{-5}]", 4000, 0, 4e-3, 1000, 0, 20);


  TH2F *h2d_Glong_Mdelta     = new TH2F ("h2d_Glong_Mdelta",     ";Galactic longitude [deg];M_{#Delta} [M_{sun}]",              5000, -240, 240, 1000, 0, 1000);
  TH2F *h2d_Glong_Mtid       = new TH2F ("h2d_Glong_Mtid",       ";Galactic longitude [deg];M_{tidal} [M_{sun}]",               5000, -240, 240, 1000, 0, 800);
  TH2F *h2d_Glong_Rtid       = new TH2F ("h2d_Glong_Rtid",       ";Galactic longitude [deg];R_{tidal} [kpc]",                   5000, -240, 240, 1000, 1e-8, 5e-2);
  TH2F *h2d_Glong_Mequdens   = new TH2F ("h2d_Glong_Mequdens",   ";Galactic longitude [deg];M_{equdens} [M_{sun}]",             5000, -240, 240, 1000, 0, 1000);
  TH2F *h2d_Glong_Requdens   = new TH2F ("h2d_Glong_Requdens",   ";Galactic longitude [deg];R_{equdens} [kpc]",                 5000, -240, 240, 1000, 1e-8, 5e-2);
  TH2F *h2d_Glong_J          = new TH2F ("h2d_Glong_J",          ";Galactic longitude [deg];J-factor [GeV^{2}cm^{-5}]",         5000, -240, 240, 1000, 1e14, 1e16);
  TH2F *h2d_Glong_JJcont     = new TH2F ("h2d_Glong_JJcont",     ";Galactic longitude [deg];J/J_{continuum} [GeV^{2}cm^{-5}]",  5000, -240, 240, 1000, 0, 20);

  TH2F *h2d_Glat_Mdelta     = new TH2F ("h2d_Glat_Mdelta",     ";Galactic latitude [deg];M_{#Delta} [M_{sun}]",              5000, -120, 120, 1000, 0, 1000);
  TH2F *h2d_Glat_Mtid       = new TH2F ("h2d_Glat_Mtid",       ";Galactic latitude [deg];M_{tidal} [M_{sun}]",               5000, -120, 120, 1000, 0, 800);
  TH2F *h2d_Glat_Rtid       = new TH2F ("h2d_Glat_Rtid",       ";Galactic latitude [deg];R_{tidal} [kpc]",                   5000, -120, 120, 1000, 1e-8, 5e-2);
  TH2F *h2d_Glat_Mequdens   = new TH2F ("h2d_Glat_Mequdens",   ";Galactic latitude [deg];M_{equdens} [M_{sun}]",             5000, -120, 120, 1000, 0, 1000);
  TH2F *h2d_Glat_Requdens   = new TH2F ("h2d_Glat_Requdens",   ";Galactic latitude [deg];R_{equdens} [kpc]",                 5000, -120, 120, 1000, 1e-8, 5e-2);
  TH2F *h2d_Glat_J          = new TH2F ("h2d_Glat_J",          ";Galactic latitude [deg];J-factor [GeV^{2}cm^{-5}]",         5000, -120, 120, 1000, 1e14, 1e16);
  TH2F *h2d_Glat_JJcont     = new TH2F ("h2d_Glat_JJcont",     ";Galactic latitude [deg];J/J_{continuum} [GeV^{2}cm^{-5}]",  5000, -120, 120, 1000, 0, 20);

  TH2F *h2d_dist_Mdelta     = new TH2F ("h2d_dist_Mdelta",     ";Distance [kpc];M_{#Delta} [M_{sun}]",               10000, 0, 250, 1000, 0, 1000);
  TH2F *h2d_dist_Mtid       = new TH2F ("h2d_dist_Mtid",       ";Distance [kpc];M_{tidal} [M_{sun}]",                10000, 0, 250, 1000, 0, 800);
  TH2F *h2d_dist_Rtid       = new TH2F ("h2d_dist_Rtid",       ";Distance [kpc];R_{tidal} [kpc]",                    10000, 0, 250, 1000, 1e-8, 5e-2);
  TH2F *h2d_dist_Mequdens   = new TH2F ("h2d_dist_Mequdens",   ";Distance [kpc];M_{equdens} [M_{sun}]",              10000, 0, 250, 1000, 0, 1000);
  TH2F *h2d_dist_Requdens   = new TH2F ("h2d_dist_Requdens",   ";Distance [kpc];R_{equdens} [kpc]",                  10000, 0, 250, 1000, 1e-8, 5e-2);
  TH2F *h2d_dist_J          = new TH2F ("h2d_dist_J",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}]",          10000, 0, 250, 1000, 1e14, 1e16);
  TH2F *h2d_dist_JJcont     = new TH2F ("h2d_dist_JJcont",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}]",   10000, 0, 250, 1000, 0, 20);

  TH2F *h2d_dgal_Mdelta     = new TH2F ("h2d_dgal_Mdelta",     ";dgalance [kpc];M_{#Delta} [M_{sun}]",               10000, 0, 250, 1000, 0, 1000);
  TH2F *h2d_dgal_Mtid       = new TH2F ("h2d_dgal_Mtid",       ";dgalance [kpc];M_{tidal} [M_{sun}]",                10000, 0, 250, 1000, 0, 800);
  TH2F *h2d_dgal_Rtid       = new TH2F ("h2d_dgal_Rtid",       ";dgalance [kpc];R_{tidal} [kpc]",                    10000, 0, 250, 1000, 1e-8, 5e-2);
  TH2F *h2d_dgal_Mequdens   = new TH2F ("h2d_dgal_Mequdens",   ";dgalance [kpc];M_{equdens} [M_{sun}]",              10000, 0, 250, 1000, 0, 1000);
  TH2F *h2d_dgal_Requdens   = new TH2F ("h2d_dgal_Requdens",   ";dgalance [kpc];R_{equdens} [kpc]",                  10000, 0, 250, 1000, 1e-8, 5e-2);
  TH2F *h2d_dgal_J          = new TH2F ("h2d_dgal_J",          ";dgalance [kpc];J-factor [GeV^{2}cm^{-5}]",          10000, 0, 250, 1000, 1e14, 1e16);
  TH2F *h2d_dgal_JJcont     = new TH2F ("h2d_dgal_JJcont",     ";dgalance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}]",   10000, 0, 250, 1000, 0, 20);

  TH2F *h2d_dgal_rs          = new TH2F ("h2d_dgal_rs",          ";dgalance [kpc];r_{S} [kpc]",        10000, 0, 250, 1000, 0, 1e1);
  TH2F *h2d_dgal_dist_rS     = new TH2F ("h2d_dgal_dist_rS",     ";dgalance [kpc];dist/r_{s}",         10000, 0, 250, 5000, 0, 1e5);
  TH2F *h2d_dgal_dgal_rS     = new TH2F ("h2d_dgal_dgal_rS",     ";dgalance [kpc];distgal/r_{s}",      10000, 0, 250, 5000, 0, 1e6);
  TH2F *h2d_dgal_Rvir_rS     = new TH2F ("h2d_dgal_Rvir_rS",     ";dgalance [kpc];R_{vir}/r_{S}",      10000, 0, 250, 5000, 0, 1e7);
  TH2F *h2d_dgal_Rtid_rS     = new TH2F ("h2d_dgal_Rtid_rS",     ";dgalance [kpc];R_{tidal}/r_{S}",    10000, 0, 250, 1000, 0, 100);
  TH2F *h2d_dgal_Req_rS      = new TH2F ("h2d_dgal_Req_rS",      ";dgalance [kpc];R_{eqdens}/r_{S}",   10000, 0, 250, 1000, 0, 500);
  TH2F *h2d_dgal_Rdelta_rS   = new TH2F ("h2d_dgal_Rdelta_rS",   ";dgalance [kpc];R_{#Delta}/r_{S}",   10000, 0, 250, 1000, 0, 200);




  TH2F *h2d_dist_J_cut0          = new TH2F ("h2d_dist_J_cut0",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}:  0-50]",            10000, 0, 250, 1000, 1e14, 1e16);
  TH2F *h2d_dist_J_cut1          = new TH2F ("h2d_dist_J_cut1",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}:  50-100]",          10000, 0, 250, 1000, 1e14, 1e16);
  TH2F *h2d_dist_J_cut2          = new TH2F ("h2d_dist_J_cut2",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 100-200]",          10000, 0, 250, 1000, 1e14, 1e16);
  TH2F *h2d_dist_J_cut3          = new TH2F ("h2d_dist_J_cut3",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 200-300]",          10000, 0, 250, 1000, 1e14, 1e16);
  TH2F *h2d_dist_J_cut4          = new TH2F ("h2d_dist_J_cut4",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 300-400]",          10000, 0, 250, 1000, 1e14, 1e16);
  TH2F *h2d_dist_JJcont_cut0     = new TH2F ("h2d_dist_JJcont_cut0",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 0-100]",     10000, 0, 250, 1000, 0, 20);
  TH2F *h2d_dist_JJcont_cut1     = new TH2F ("h2d_dist_JJcont_cut1",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 0-100]",     10000, 0, 250, 1000, 0, 20);
  TH2F *h2d_dist_JJcont_cut2     = new TH2F ("h2d_dist_JJcont_cut2",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 100-200]",   10000, 0, 250, 1000, 0, 20);
  TH2F *h2d_dist_JJcont_cut3     = new TH2F ("h2d_dist_JJcont_cut3",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 200-300]",   10000, 0, 250, 1000, 0, 20);
  TH2F *h2d_dist_JJcont_cut4     = new TH2F ("h2d_dist_JJcont_cut4",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}: 300-400]",   10000, 0, 250, 1000, 0, 20);




  TH2F *h2d_Mdelta_rhos       = new TH2F ("h2d_Mdelta_rhos",       ";M_{#Delta} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",        20000, 0, 200, 1000, 1e8, 2e9);
  TH2F *h2d_Mdelta_rs         = new TH2F ("h2d_Mdelta_rs",         ";M_{#Delta} [M_{sun}];r_{s} [kpc]",                       20000, 0, 200, 1000, 0, 1e-3);
  TH2F *h2d_Mdelta_Mtid       = new TH2F ("h2d_Mdelta_Mtid",       ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}]",               20000, 0, 200, 1000, 0, 800);
  TH2F *h2d_Mdelta_Rtid       = new TH2F ("h2d_Mdelta_Rtid",       ";M_{#Delta} [M_{sun}];R_{tidal} [kpc]",                   20000, 0, 200, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mdelta_Mequdens   = new TH2F ("h2d_Mdelta_Mequdens",   ";M_{#Delta} [M_{sun}];M_{equdens} [M_{sun}]",             20000, 0, 200, 1000, 0, 1000);
  TH2F *h2d_Mdelta_Requdens   = new TH2F ("h2d_Mdelta_Requdens",   ";M_{#Delta} [M_{sun}];R_{equdens} [kpc]",                 20000, 0, 200, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mdelta_J          = new TH2F ("h2d_Mdelta_J",          ";M_{#Delta} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",         20000, 0, 200, 1000, 1e14, 1e16);
  TH2F *h2d_Mdelta_JJcont     = new TH2F ("h2d_Mdelta_JJcont",     ";M_{#Delta} [M_{sun}];J/J_{continuum} [GeV^{2}cm^{-5}]",  20000, 0, 200, 1000, 0, 20);

  TH2F *h2d_Mtid_rhos       = new TH2F ("h2d_Mtid_rhos",       ";M_{tidal} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",       20000, 0, 1000, 1000, 1e8, 2e9);
  TH2F *h2d_Mtid_rs         = new TH2F ("h2d_Mtid_rs",         ";M_{tidal} [M_{sun}];r_{s} [kpc]",                      20000, 0, 1000, 1000, 0, 1e-3);
  TH2F *h2d_Mtid_Mdelta     = new TH2F ("h2d_Mtid_Mdelta",     ";M_{tidal} [M_{sun}];M_{#Delta} [M_{sun}]",             20000, 0, 1000, 1000, 0, 800);
  TH2F *h2d_Mtid_Rtid       = new TH2F ("h2d_Mtid_Rtid",       ";M_{tidal} [M_{sun}];R_{tidal} [kpc]",                  20000, 0, 1000, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mtid_Mequdens   = new TH2F ("h2d_Mtid_Mequdens",   ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}]",            20000, 0, 1000, 1000, 0, 1000);
  TH2F *h2d_Mtid_Requdens   = new TH2F ("h2d_Mtid_Requdens",   ";M_{tidal} [M_{sun}];R_{equdens} [kpc]",                20000, 0, 1000, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mtid_J          = new TH2F ("h2d_Mtid_J",          ";M_{tidal} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",        20000, 0, 1000, 1000, 1e14, 1e16);
  TH2F *h2d_Mtid_JJcont     = new TH2F ("h2d_Mtid_JJcont",     ";M_{tidal} [M_{sun}];J/J_{continuum} [GeV^{2}cm^{-5}]", 20000, 0, 1000, 1000, 0, 20);

  TH2F *h2d_Rtid_rhos       = new TH2F ("h2d_Rtid_rhos",       ";R_{tidal} [kpc];#rho_{s} [M_{sun}/kpc^{3}]",        40000, 0.002, 0.05, 1000, 1e8, 2e9);
  TH2F *h2d_Rtid_rs         = new TH2F ("h2d_Rtid_rs",         ";R_{tidal} [kpc];r_{s} [kpc]",                       40000, 0.002, 0.05, 1000, 0, 1e-3);
  TH2F *h2d_Rtid_Mdelta     = new TH2F ("h2d_Rtid_Mdelta",     ";R_{tidal} [kpc];M_{#Delta} [M_{sun}]",              40000, 0.002, 0.05, 1000, 0, 800);
  TH2F *h2d_Rtid_Mtid       = new TH2F ("h2d_Rtid_Mtid",       ";R_{tidal} [kpc];M_{tidal} [M_{sun}]",               40000, 0.002, 0.05, 1000, 0, 1000);
  TH2F *h2d_Rtid_Mequdens   = new TH2F ("h2d_Rtid_Mequdens",   ";R_{tidal} [kpc];M_{equdens} [M_{sun}]",             40000, 0.002, 0.05, 1000, 0, 1000);
  TH2F *h2d_Rtid_Requdens   = new TH2F ("h2d_Rtid_Requdens",   ";R_{tidal} [kpc];R_{equdens} [kpc]",                 40000, 0.002, 0.05, 1000, 1e-8, 5e-2);
  TH2F *h2d_Rtid_J          = new TH2F ("h2d_Rtid_J",          ";R_{tidal} [kpc];J-factor [GeV^{2}cm^{-5}]",         40000, 0.002, 0.05, 1000, 1e14, 1e16);
  TH2F *h2d_Rtid_JJcont     = new TH2F ("h2d_Rtid_JJcont",     ";R_{tidal} [kpc];J/J_{continuum} [GeV^{2}cm^{-5}]",  40000, 0.002, 0.05, 1000, 0, 20);

  TH2F *h2d_Requdens_rhos       = new TH2F ("h2d_Requdens_rhos",       ";R_{equdensal} [kpc];#rho_{s} [M_{sun}/kpc^{3}]",       40000, 0.002, 0.05, 1000, 1e8, 2e9);
  TH2F *h2d_Requdens_rs         = new TH2F ("h2d_Requdens_rs",         ";R_{equdensal} [kpc];r_{s} [kpc]",                      40000, 0.002, 0.05, 1000, 0, 1e-3);
  TH2F *h2d_Requdens_Mdelta     = new TH2F ("h2d_Requdens_Mdelta",     ";R_{equdensal} [kpc];M_{#Delta} [M_{sun}]",             40000, 0.002, 0.05, 1000, 0, 1000);
  TH2F *h2d_Requdens_Mequdens   = new TH2F ("h2d_Requdens_Mequdens",   ";R_{equdensal} [kpc];M_{equdens} [M_{sun}]",            40000, 0.002, 0.05, 1000, 0, 1000);
  TH2F *h2d_Requdens_Rtid       = new TH2F ("h2d_Requdens_Rtid",       ";R_{equdensal} [kpc];R_{tidal} [kpc]",                  40000, 0.002, 0.05, 1000, 1e-8, 5e-2);
  TH2F *h2d_Requdens_Mtid       = new TH2F ("h2d_Requdens_Mtid",       ";R_{equdensal} [kpc];M_{tidal} [M_{sun}]",              40000, 0.002, 0.05, 1000, 0, 800);
  TH2F *h2d_Requdens_J          = new TH2F ("h2d_Requdens_J",          ";R_{equdensal} [kpc];J-factor [GeV^{2}cm^{-5}]",        40000, 0.002, 0.05, 1000, 1e14, 1e16);
  TH2F *h2d_Requdens_JJcont     = new TH2F ("h2d_Requdens_JJcont",     ";R_{equdensal} [kpc];J/J_{continuum} [GeV^{2}cm^{-5}]", 40000, 0.002, 0.05, 1000, 0, 20);

  TH2F *h2d_Mequdens_rhos       = new TH2F ("h2d_Mequdens_rhos",       ";M_{equdensal} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",       40000, 0, 500, 1000, 1e8, 2e9);
  TH2F *h2d_Mequdens_rs         = new TH2F ("h2d_Mequdens_rs",         ";M_{equdensal} [M_{sun}];r_{s} [kpc]",                      40000, 0, 500, 1000, 0, 1e-3);
  TH2F *h2d_Mequdens_Mdelta     = new TH2F ("h2d_Mequdens_Mdelta",     ";M_{equdensal} [M_{sun}];M_{#Delta} [M_{sun}]",             40000, 0, 500, 1000, 0, 1000);
  TH2F *h2d_Mequdens_Requdens   = new TH2F ("h2d_Mequdens_Requdens",   ";M_{equdensal} [M_{sun}];R_{equdens} [kpc]",                40000, 0, 500, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mequdens_Mtid       = new TH2F ("h2d_Mequdens_Mtid",       ";M_{equdensal} [M_{sun}];M_{tidal} [M_{sun}]",              40000, 0, 500, 1000, 0, 800);
  TH2F *h2d_Mequdens_Rtid       = new TH2F ("h2d_Mequdens_Rtid",       ";M_{equdensal} [M_{sun}];R_{tidal} [kpc]",                  40000, 0, 500, 1000, 1e-8, 5e-2);
  TH2F *h2d_Mequdens_J          = new TH2F ("h2d_Mequdens_J",          ";M_{equdensal} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",        40000, 0, 500, 1000, 1e14, 1e16);
  TH2F *h2d_Mequdens_JJcont     = new TH2F ("h2d_Mequdens_JJcont",     ";M_{equdensal} [M_{sun}];J/J_{continuum} [GeV^{2}cm^{-5}]", 40000, 0, 500, 1000, 0, 20);


  TProfile2D *prf2d_dist_J_Mdelta        = new TProfile2D ("prf2d_dist_J_Mdelta",         ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{#Delta} [M_{sun}]",               5000, 0, 40, 1000, 1e14, 1e16, 0, 1000);
  TProfile2D *prf2d_dist_J_Mtidal        = new TProfile2D ("prf2d_dist_J_Mtidal",         ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}]",                 5000, 0, 40, 1000, 1e14, 1e16, 0, 800);
  TProfile2D *prf2d_dist_J_Mequdens      = new TProfile2D ("prf2d_dist_J_Mequdens",       ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];M_{equdens}[M_{sun}]",               5000, 0, 40, 1000, 1e14, 1e16, 0, 800);
  TProfile2D *prf2d_dist_J_Rtidal        = new TProfile2D ("prf2d_dist_J_Rtidal",         ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];R_{tidal}[kpc]",                     5000, 0, 40, 1000, 1e14, 1e16, 1e-8, 5e-2);
  TProfile2D *prf2d_dist_J_Requdens      = new TProfile2D ("prf2d_dist_J_Requdens",       ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];R_{equdens}[kpc]",                   5000, 0, 40, 1000, 1e14, 1e16, 1e-8, 5e-2);
  TProfile2D *prf2d_dist_J_rhos          = new TProfile2D ("prf2d_dist_J_rhos",           ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];#rho_{s} [M_{sun}/kpc^{3}]",         5000, 0, 40, 1000, 1e14, 1e16, 0, 1.0e9);
  TProfile2D *prf2d_dist_J_Glong         = new TProfile2D ("prf2d_dist_J_Glong",          ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];Galactic longitude [deg]",           5000, 0, 40, 1000, 1e14, 1e16, -240, 240);
  TProfile2D *prf2d_dist_J_Glat          = new TProfile2D ("prf2d_dist_J_Glat",           ";Distance [kpc];J-factor [GeV^{2}cm^{-5}];Galactic latitude [deg]",            5000, 0, 40, 1000, 1e14, 1e16, -120, 120);

  TProfile2D *prf2d_dist_JJcont_Mtidal   = new TProfile2D ("prf2d_dist_JJcont_Mtidal",    ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{tidal}[M_{sun}]",          5000, 0, 40, 1000, 0, 20, 0, 800);
  TProfile2D *prf2d_dist_JJcont_Mequdens = new TProfile2D ("prf2d_dist_JJcont_Mequdens",  ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{equdens}[M_{sun}]",        5000, 0, 40, 1000, 0, 20, 0, 800);
  TProfile2D *prf2d_dist_JJcont_Mdelta   = new TProfile2D ("prf2d_dist_JJcont_Mdelta",    ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];M_{#Delta} [M_{sun}]",        5000, 0, 40, 1000, 0, 20, 0, 1000);
  TProfile2D *prf2d_dist_JJcont_Rtidal   = new TProfile2D ("prf2d_dist_JJcont_Rtidal",    ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];R_{tidal}[kpc]",              5000, 0, 40, 1000, 0, 20, 1e-8, 5e-2);
  TProfile2D *prf2d_dist_JJcont_Requdens = new TProfile2D ("prf2d_dist_JJcont_Requdens",  ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];R_{equdens}[kpc]",            5000, 0, 40, 1000, 0, 20, 1e-8, 5e-2);
  TProfile2D *prf2d_dist_JJcont_rhos     = new TProfile2D ("prf2d_dist_JJcont_rhos",      ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];#rho_{s} [M_{sun}/kpc^{3}]",  5000, 0, 40, 1000, 0, 20, 0, 1.0e9);
  TProfile2D *prf2d_dist_JJcont_Glong    = new TProfile2D ("prf2d_dist_JJcont_Glong",     ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];Galactic longitude [deg]",    5000, 0, 40, 1000, 0, 20, -240, 240);
  TProfile2D *prf2d_dist_JJcont_Glat     = new TProfile2D ("prf2d_dist_JJcont_Glat",      ";Distance [kpc];J/J_{continuum} [GeV^{2}cm^{-5}];Galactic latitude [deg]",     5000, 0, 40, 1000, 0, 20, -120, 120);

  TProfile2D *prf2d_Glon_Glat_Mdelta     = new TProfile2D ("prf2d_Glon_Glat_Mdelta",      ";Galactic longitude [deg];Galactic latitude [deg];M_{#Delta} [M_{sun}]",       1000, -200, 200, 1000, -100, 100, 0, 1000);
  TProfile2D *prf2d_Glon_Glat_Mtidal     = new TProfile2D ("prf2d_Glon_Glat_Mtidal",      ";Galactic longitude [deg];Galactic latitude [deg];M_{tidal} [M_{sun}]",        1000, -200, 200, 1000, -100, 100, 0, 1000);
  TProfile2D *prf2d_Glon_Glat_Mequdens   = new TProfile2D ("prf2d_Glon_Glat_Mequdens",    ";Galactic longitude [deg];Galactic latitude [deg];M_{equdens} [M_{sun}]",      1000, -200, 200, 1000, -100, 100, 0, 1000);
  TProfile2D *prf2d_Glon_Glat_Rtidal     = new TProfile2D ("prf2d_Glon_Glat_Rtidal",      ";Galactic longitude [deg];Galactic latitude [deg];R_{tidal} [kpc]",            1000, -200, 200, 1000, -100, 100, 1e-8, 5e-2);
  TProfile2D *prf2d_Glon_Glat_Requdens   = new TProfile2D ("prf2d_Glon_Glat_Requdens",    ";Galactic longitude [deg];Galactic latitude [deg];R_{equdens} [kpc]",          1000, -200, 200, 1000, -100, 100, 1e-8, 5e-2);
  TProfile2D *prf2d_Glon_Glat_J          = new TProfile2D ("prf2d_Glon_Glat_J",           ";Galactic longitude [deg];Galactic latitude [deg];J-factor [GeV^{2}cm^{-5}]",  1000, -200, 200, 1000, -100, 100, 1e14, 1e16);
  TProfile2D *prf2d_Glon_Glat_JJcont     = new TProfile2D ("prf2d_Glon_Glat_JJcont",      ";Galactic longitude [deg];Galactic latitude [deg];J/J_{cont} [GeV^{2}cm^{-5}]",1000, -200, 200, 1000, -100, 100, 0, 20);
  TProfile2D *prf2d_Glon_Glat_rhos       = new TProfile2D ("prf2d_Glon_Glat_rhos",        ";Galactic longitude [deg];Galactic latitude [deg];#rho_{s} [M_{sun}/kpc^{3}]", 1000, -200, 200, 1000, -100, 100, 0, 1.5e9);
  TProfile2D *prf2d_Glon_Glat_rs         = new TProfile2D ("prf2d_Glon_Glat_rs",          ";Galactic longitude [deg];Galactic latitude [deg];r_{s} [kpc]",                1000, -200, 200, 1000, -100, 100, 0, 4e-3);

  TProfile2D *prf2d_Mdelta_Mtidal_rs         = new TProfile2D ("prf2d_Mdelta_Mtidal_rs",          ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 800, 0, 4e-3);
  TProfile2D *prf2d_Mdelta_Mtidal_rhos       = new TProfile2D ("prf2d_Mdelta_Mtidal_rhos",        ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 800, 0, 1.5e9);
  TProfile2D *prf2d_Mdelta_Mtidal_Rtidal     = new TProfile2D ("prf2d_Mdelta_Mtidal_Rtidal",      ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 800, 1e-8, 5e-2);
  TProfile2D *prf2d_Mdelta_Mtidal_Mequdens   = new TProfile2D ("prf2d_Mdelta_Mtidal_Mequdens",    ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];M_{equdens} [M_{sun}]",         1000,  0, 800, 1000, 0, 800, 0, 800);
  TProfile2D *prf2d_Mdelta_Mtidal_Requdens   = new TProfile2D ("prf2d_Mdelta_Mtidal_Requdens",    ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 800, 1e-8, 5e-2);
  TProfile2D *prf2d_Mdelta_Mtidal_Jfact      = new TProfile2D ("prf2d_Mdelta_Mtidal_Jfact",       ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 800, 1e14, 1e16);
  TProfile2D *prf2d_Mdelta_Mtidal_JJ         = new TProfile2D ("prf2d_Mdelta_Mtidal_JJ",          ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 800, 0, 20);
  TProfile2D *prf2d_Mdelta_Mtidal_dist       = new TProfile2D ("prf2d_Mdelta_Mtidal_dist",        ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];distance [kpc]",                1000,  0, 800, 1000, 0, 800, 0, 40);
  TProfile2D *prf2d_Mdelta_Mtidal_Glong      = new TProfile2D ("prf2d_Mdelta_Mtidal_Glong",       ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 800, -200, 200);
  TProfile2D *prf2d_Mdelta_Mtidal_Glat       = new TProfile2D ("prf2d_Mdelta_Mtidal_Glat",        ";M_{#Delta} [M_{sun}];M_{tidal} [M_{sun}];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 800, -100, 100);

  TProfile2D *prf2d_Mtidal_Mequdens_rs         = new TProfile2D ("prf2d_Mtidal_Mequdens_rs",          ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 800, 0, 4e-3);
  TProfile2D *prf2d_Mtidal_Mequdens_rhos       = new TProfile2D ("prf2d_Mtidal_Mequdens_rhos",        ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 800, 0, 1.5e9);
  TProfile2D *prf2d_Mtidal_Mequdens_Mdelta     = new TProfile2D ("prf2d_Mtidal_Mequdens_Mdelta",      ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];M_{#Delta} [M_{sun}]",          1000,  0, 800, 1000, 0, 800, 0, 800);
  TProfile2D *prf2d_Mtidal_Mequdens_Rtidal     = new TProfile2D ("prf2d_Mtidal_Mequdens_Rtidal",      ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 800, 1e-8, 5e-2);
  TProfile2D *prf2d_Mtidal_Mequdens_Requdens   = new TProfile2D ("prf2d_Mtidal_Mequdens_Requdens",    ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 800, 1e-8, 5e-2);
  TProfile2D *prf2d_Mtidal_Mequdens_Jfact      = new TProfile2D ("prf2d_Mtidal_Mequdens_Jfact",       ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 800, 1e14, 1e16);
  TProfile2D *prf2d_Mtidal_Mequdens_JJ         = new TProfile2D ("prf2d_Mtidal_Mequdens_JJ",          ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 800, 0, 20);
  TProfile2D *prf2d_Mtidal_Mequdens_dist       = new TProfile2D ("prf2d_Mtidal_Mequdens_dist",        ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];distance [kpc]",                1000,  0, 800, 1000, 0, 800, 0, 40);
  TProfile2D *prf2d_Mtidal_Mequdens_Glong      = new TProfile2D ("prf2d_Mtidal_Mequdens_Glong",       ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 800, -200, 200);
  TProfile2D *prf2d_Mtidal_Mequdens_Glat       = new TProfile2D ("prf2d_Mtidal_Mequdens_Glat",        ";M_{tidal} [M_{sun}];M_{equdens} [M_{sun}];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 800, -100, 100);



  TProfile2D *prf2d_Mdelta_dist_rs         = new TProfile2D ("prf2d_Mdelta_dist_rs",          ";M_{#Delta} [M_{sun}];Distance [kpc];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 250, 0, 4e-3);
  TProfile2D *prf2d_Mdelta_dist_rhos       = new TProfile2D ("prf2d_Mdelta_dist_rhos",        ";M_{#Delta} [M_{sun}];Distance [kpc];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 250, 0, 1.5e9);
  TProfile2D *prf2d_Mdelta_dist_Rtidal     = new TProfile2D ("prf2d_Mdelta_dist_Rtidal",      ";M_{#Delta} [M_{sun}];Distance [kpc];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mdelta_dist_Mequdens   = new TProfile2D ("prf2d_Mdelta_dist_Mequdens",    ";M_{#Delta} [M_{sun}];Distance [kpc];M_{equdens} [M_{sun}]",         1000,  0, 800, 1000, 0, 250, 0, 800);
  TProfile2D *prf2d_Mdelta_dist_Requdens   = new TProfile2D ("prf2d_Mdelta_dist_Requdens",    ";M_{#Delta} [M_{sun}];Distance [kpc];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mdelta_dist_Jfact      = new TProfile2D ("prf2d_Mdelta_dist_Jfact",       ";M_{#Delta} [M_{sun}];Distance [kpc];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 250, 1e14, 1e16);
  TProfile2D *prf2d_Mdelta_dist_JJ         = new TProfile2D ("prf2d_Mdelta_dist_JJ",          ";M_{#Delta} [M_{sun}];Distance [kpc];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 250, 0, 20);
  TProfile2D *prf2d_Mdelta_dist_Glong      = new TProfile2D ("prf2d_Mdelta_dist_Glong",       ";M_{#Delta} [M_{sun}];Distance [kpc];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 250, -200, 200);
  TProfile2D *prf2d_Mdelta_dist_Glat       = new TProfile2D ("prf2d_Mdelta_dist_Glat",        ";M_{#Delta} [M_{sun}];Distance [kpc];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 250, -100, 100);

  TProfile2D *prf2d_Mtidal_dist_rs         = new TProfile2D ("prf2d_Mtidal_dist_rs",          ";M_{tidal} [M_{sun}];Distance [kpc];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 250, 0, 4e-3);
  TProfile2D *prf2d_Mtidal_dist_rhos       = new TProfile2D ("prf2d_Mtidal_dist_rhos",        ";M_{tidal} [M_{sun}];Distance [kpc];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 250, 0, 1.5e9);
  TProfile2D *prf2d_Mtidal_dist_Mdelta     = new TProfile2D ("prf2d_Mtidal_dist_Mdelta",      ";M_{tidal} [M_{sun}];Distance [kpc];M_{#Delta} [M_{sun}]",          1000,  0, 800, 1000, 0, 250, 0, 800);
  TProfile2D *prf2d_Mtidal_dist_Rtidal     = new TProfile2D ("prf2d_Mtidal_dist_Rtidal",      ";M_{tidal} [M_{sun}];Distance [kpc];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mtidal_dist_Requdens   = new TProfile2D ("prf2d_Mtidal_dist_Requdens",    ";M_{tidal} [M_{sun}];Distance [kpc];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mtidal_dist_Jfact      = new TProfile2D ("prf2d_Mtidal_dist_Jfact",       ";M_{tidal} [M_{sun}];Distance [kpc];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 250, 1e14, 1e16);
  TProfile2D *prf2d_Mtidal_dist_JJ         = new TProfile2D ("prf2d_Mtidal_dist_JJ",          ";M_{tidal} [M_{sun}];Distance [kpc];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 250, 0, 20);
  TProfile2D *prf2d_Mtidal_dist_Glong      = new TProfile2D ("prf2d_Mtidal_dist_Glong",       ";M_{tidal} [M_{sun}];Distance [kpc];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 250, -200, 200);
  TProfile2D *prf2d_Mtidal_dist_Glat       = new TProfile2D ("prf2d_Mtidal_dist_Glat",        ";M_{tidal} [M_{sun}];Distance [kpc];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 250, -100, 100);

  TProfile2D *prf2d_Mequdens_dist_rs         = new TProfile2D ("prf2d_Mequdens_dist_rs",          ";M_{tidal} [M_{sun}];Distance [kpc];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 250, 0, 4e-3);
  TProfile2D *prf2d_Mequdens_dist_rhos       = new TProfile2D ("prf2d_Mequdens_dist_rhos",        ";M_{tidal} [M_{sun}];Distance [kpc];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 250, 0, 1.5e9);
  TProfile2D *prf2d_Mequdens_dist_Mdelta     = new TProfile2D ("prf2d_Mequdens_dist_Mdelta",      ";M_{tidal} [M_{sun}];Distance [kpc];M_{#Delta} [M_{sun}]",          1000,  0, 800, 1000, 0, 250, 0, 800);
  TProfile2D *prf2d_Mequdens_dist_Rtidal     = new TProfile2D ("prf2d_Mequdens_dist_Rtidal",      ";M_{tidal} [M_{sun}];Distance [kpc];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mequdens_dist_Requdens   = new TProfile2D ("prf2d_Mequdens_dist_Requdens",    ";M_{tidal} [M_{sun}];Distance [kpc];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mequdens_dist_Jfact      = new TProfile2D ("prf2d_Mequdens_dist_Jfact",       ";M_{tidal} [M_{sun}];Distance [kpc];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 250, 1e14, 1e16);
  TProfile2D *prf2d_Mequdens_dist_JJ         = new TProfile2D ("prf2d_Mequdens_dist_JJ",          ";M_{tidal} [M_{sun}];Distance [kpc];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 250, 0, 20);
  TProfile2D *prf2d_Mequdens_dist_Glong      = new TProfile2D ("prf2d_Mequdens_dist_Glong",       ";M_{tidal} [M_{sun}];Distance [kpc];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 250, -200, 200);
  TProfile2D *prf2d_Mequdens_dist_Glat       = new TProfile2D ("prf2d_Mequdens_dist_Glat",        ";M_{tidal} [M_{sun}];Distance [kpc];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 250, -100, 100);


  TProfile2D *prf2d_Mdelta_dgal_rs         = new TProfile2D ("prf2d_Mdelta_dgal_rs",          ";M_{#Delta} [M_{sun}];DistGal [kpc];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 250, 0, 4e-3);
  TProfile2D *prf2d_Mdelta_dgal_rhos       = new TProfile2D ("prf2d_Mdelta_dgal_rhos",        ";M_{#Delta} [M_{sun}];DistGal [kpc];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 250, 0, 1.5e9);
  TProfile2D *prf2d_Mdelta_dgal_Rtidal     = new TProfile2D ("prf2d_Mdelta_dgal_Rtidal",      ";M_{#Delta} [M_{sun}];DistGal [kpc];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mdelta_dgal_Mequdens   = new TProfile2D ("prf2d_Mdelta_dgal_Mequdens",    ";M_{#Delta} [M_{sun}];DistGal [kpc];M_{equdens} [M_{sun}]",         1000,  0, 800, 1000, 0, 250, 0, 800);
  TProfile2D *prf2d_Mdelta_dgal_Requdens   = new TProfile2D ("prf2d_Mdelta_dgal_Requdens",    ";M_{#Delta} [M_{sun}];DistGal [kpc];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mdelta_dgal_Jfact      = new TProfile2D ("prf2d_Mdelta_dgal_Jfact",       ";M_{#Delta} [M_{sun}];DistGal [kpc];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 250, 1e14, 1e16);
  TProfile2D *prf2d_Mdelta_dgal_JJ         = new TProfile2D ("prf2d_Mdelta_dgal_JJ",          ";M_{#Delta} [M_{sun}];DistGal [kpc];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 250, 0, 20);
  TProfile2D *prf2d_Mdelta_dgal_Glong      = new TProfile2D ("prf2d_Mdelta_dgal_Glong",       ";M_{#Delta} [M_{sun}];DistGal [kpc];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 250, -200, 200);
  TProfile2D *prf2d_Mdelta_dgal_Glat       = new TProfile2D ("prf2d_Mdelta_dgal_Glat",        ";M_{#Delta} [M_{sun}];DistGal [kpc];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 250, -100, 100);

  TProfile2D *prf2d_Mtidal_dgal_rs         = new TProfile2D ("prf2d_Mtidal_dgal_rs",          ";M_{tidal} [M_{sun}];DistGal [kpc];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 250, 0, 4e-3);
  TProfile2D *prf2d_Mtidal_dgal_rhos       = new TProfile2D ("prf2d_Mtidal_dgal_rhos",        ";M_{tidal} [M_{sun}];DistGal [kpc];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 250, 0, 1.5e9);
  TProfile2D *prf2d_Mtidal_dgal_Mdelta     = new TProfile2D ("prf2d_Mtidal_dgal_Mdelta",      ";M_{tidal} [M_{sun}];DistGal [kpc];M_{#Delta} [M_{sun}]",          1000,  0, 800, 1000, 0, 250, 0, 800);
  TProfile2D *prf2d_Mtidal_dgal_Rtidal     = new TProfile2D ("prf2d_Mtidal_dgal_Rtidal",      ";M_{tidal} [M_{sun}];DistGal [kpc];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mtidal_dgal_Requdens   = new TProfile2D ("prf2d_Mtidal_dgal_Requdens",    ";M_{tidal} [M_{sun}];DistGal [kpc];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mtidal_dgal_Jfact      = new TProfile2D ("prf2d_Mtidal_dgal_Jfact",       ";M_{tidal} [M_{sun}];DistGal [kpc];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 250, 1e14, 1e16);
  TProfile2D *prf2d_Mtidal_dgal_JJ         = new TProfile2D ("prf2d_Mtidal_dgal_JJ",          ";M_{tidal} [M_{sun}];DistGal [kpc];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 250, 0, 20);
  TProfile2D *prf2d_Mtidal_dgal_Glong      = new TProfile2D ("prf2d_Mtidal_dgal_Glong",       ";M_{tidal} [M_{sun}];DistGal [kpc];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 250, -200, 200);
  TProfile2D *prf2d_Mtidal_dgal_Glat       = new TProfile2D ("prf2d_Mtidal_dgal_Glat",        ";M_{tidal} [M_{sun}];DistGal [kpc];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 250, -100, 100);

  TProfile2D *prf2d_Mequdens_dgal_rs         = new TProfile2D ("prf2d_Mequdens_dgal_rs",          ";M_{tidal} [M_{sun}];DistGal [kpc];r_{s} [kpc]",                   1000,  0, 800, 1000, 0, 250, 0, 4e-3);
  TProfile2D *prf2d_Mequdens_dgal_rhos       = new TProfile2D ("prf2d_Mequdens_dgal_rhos",        ";M_{tidal} [M_{sun}];DistGal [kpc];#rho_{s} [M_{sun}/kpc^{3}]",    1000,  0, 800, 1000, 0, 250, 0, 1.5e9);
  TProfile2D *prf2d_Mequdens_dgal_Mdelta     = new TProfile2D ("prf2d_Mequdens_dgal_Mdelta",      ";M_{tidal} [M_{sun}];DistGal [kpc];M_{#Delta} [M_{sun}]",          1000,  0, 800, 1000, 0, 250, 0, 800);
  TProfile2D *prf2d_Mequdens_dgal_Rtidal     = new TProfile2D ("prf2d_Mequdens_dgal_Rtidal",      ";M_{tidal} [M_{sun}];DistGal [kpc];R_{tidal} [kpc]",               1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mequdens_dgal_Requdens   = new TProfile2D ("prf2d_Mequdens_dgal_Requdens",    ";M_{tidal} [M_{sun}];DistGal [kpc];R_{equdens} [kpc]",             1000,  0, 800, 1000, 0, 250, 1e-8, 5e-2);
  TProfile2D *prf2d_Mequdens_dgal_Jfact      = new TProfile2D ("prf2d_Mequdens_dgal_Jfact",       ";M_{tidal} [M_{sun}];DistGal [kpc];J-factor [GeV^{2}cm^{-5}]",     1000,  0, 800, 1000, 0, 250, 1e14, 1e16);
  TProfile2D *prf2d_Mequdens_dgal_JJ         = new TProfile2D ("prf2d_Mequdens_dgal_JJ",          ";M_{tidal} [M_{sun}];DistGal [kpc];J/J_{cont} [GeV^{2}cm^{-5}]",   1000,  0, 800, 1000, 0, 250, 0, 20);
  TProfile2D *prf2d_Mequdens_dgal_Glong      = new TProfile2D ("prf2d_Mequdens_dgal_Glong",       ";M_{tidal} [M_{sun}];DistGal [kpc];Gal.longitude [deg]",           1000,  0, 800, 1000, 0, 250, -200, 200);
  TProfile2D *prf2d_Mequdens_dgal_Glat       = new TProfile2D ("prf2d_Mequdens_dgal_Glat",        ";M_{tidal} [M_{sun}];DistGal [kpc];Gal.latitude [deg]",            1000,  0, 800, 1000, 0, 250, -100, 100);



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
  std::ifstream ifstr(file, std::ifstream::in);
  if (ifstr.fail()) {
    std::cerr << "\n Could not open " << file << "\n\n";
    return;
  }
  else {
     std::cout << "File is successfully opened.\n" << std::endl;
     //cout << "ifstream.rdbuf  = " << ifstr.rdbuf()   << endl;
     std::cout << "ifstream.get    = " << ifstr.get()   << std::endl;
  }

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
        if (counter>1){
          //break;

          if (line.length() > 184 || line.length() < 170)
            continue;

          if (counter>50000 && counter<500001){
          //if (counter>50000){
                //std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
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
          else if (counter>500000){
                //std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
                //std::getline(ss,dummy,' ');      //std::cout << "dummy: " << dummy << std::endl;
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

                        if (DMclumID.length() == 0){
                           std::getline(ss,DMclumID,' ');      //std::cout << "Glong: " << Glong << " --> length:" << Glong.length()  << std::endl;

                        }
                    }
                  }
                 }
                }

            }
                if (DMclumID.length() == 0)
                  continue;

                std::cout << "DMclumID: " << DMclumID << "  length:" << DMclumID.length() << std::endl;

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
                //float dm_glong = std::stof(Glong);
                //float dm_glat  = std::stof(Glat);
                //float dm_dist  = std::stof(DM_clump_dist);

                mDMclumpID          = dm_ID;
                clump_location.l    = std::stof(Glat);
                clump_location.beta = std::stof(Glong);
                clump_location.dist = std::stof(DM_clump_dist);

                //if (std::stof(Glong) < 200 && std::stof(Glong) < -200)
                //if (Glong.length() != 0)

                if (std::abs(std::stof(Glong))<1000)
                  mGlong = (double)std::stof(Glong);
                //mGlong = clump_location.beta;

                if (std::abs(std::stof(Glat))<1000)
                  mGlat = (double)std::stof(Glat);
                 //mGlat = clump_location.l;

                if (std::stof(DM_clump_dist)<1000)
                // && std::stof(DM_clump_dist)>0.0)
                //if (DM_clump_dist.length() != 0 && std::stof(DM_clump_dist) > 0.0)
                if (std::stof(DM_clump_dist) > 0.0)
                mDist  = std::abs((double)std::stof(DM_clump_dist));
                //mDist  = clump_location.dist;

                mZ      = std::stof(z);
                mRdelta = std::stof(Rdelta);

                if (std::abs(std::stof(rhos))<1e10  && std::stof(rhos)>0 ){
                  mRhos   = std::stof(rhos);
                  mRhos2  = std::stof(rhos)*std::stof(rhos);
                }


                if (rs.length() != 0)
                  mRs     = (double)std::stof(rs);

                m_alphaS = std::asin(mRs/mDist);
                if (std::asin(mRs/mDist) > 0.)
                  m_log10_alphaS = std::log10(std::asin(mRs/mDist));


                //float dm_z      = std::stof(z);
                //float dm_Rdelta = std::stof(Rdelta);
                //float dm_rhos   = std::stof(rhos);
                //float dm_rs     = std::stof(rs);

                //float dm_p1     = std::stof(DMprof_p1);
                //float dm_p2     = std::stof(DMprof_p2);
                //float dm_p3     = std::stof(DMprof_p3);

                //float dm_Jfactor = std::stof(DM_J);
                //float dm_JJcont  = std::stof(DM_JJcont);

//                 float dm_Mdelta = 0.;//std::stof(Mdelta);
//                 float dm_Mtid   = 0.;//std::stof(Mtid);
//                 float dm_Rtid   = 0.;//std::stof(Rtid);
//
//                 float dm_Mequdens = 0.;//std::stof(Mequdens);
//                 float dm_Requdens = 0.;//std::stof(Requdens);
//                 float dm_Dgal     = 0.;//std::stof(Dgal);
                float dm_Rtidal;
                if (FLAG==1){
//                   dm_Mdelta = std::stof(Mdelta);
//                   dm_Mtid   = std::stof(Mtid);
//                   dm_Rtid   = std::stof(Rtid);
//
//                   dm_Mequdens = std::stof(Mequdens);
//                   dm_Requdens = std::stof(Requdens);
//                   dm_Dgal     = std::stof(Dgal);
                 dm_Rtidal  = std::stof(Rtid);

                 mMdelta    = std::stof(Mdelta);
                 mMtidal    = std::stof(Mtid);
                 //mRtidal    = std::stof(Rtid);
                 mRtidal    = dm_Rtidal;
                 mMequdens  = std::stof(Mequdens);
                 mRequdens  = std::stof(Requdens);
                 if (std::stof(Dgal)<1000 && std::stof(Dgal)>0.0)
                 //if (std::stof(Dgal) > 0.0)
                 mDgal      = std::abs((double) std::stof(Dgal) );
                 //mDgal      = std::stof(Dgal);

                 if (std::stof(Mequdens) > 0 && std::stof(Mdelta) > 0 ) mMeqdensMhalo  = std::stof(Mequdens) / std::stof(Mdelta);
                 if (std::stof(Mtid) > 0 && std::stof(Mdelta) > 0 )     mMtidalMhalo   = std::stof(Mtid) / std::stof(Mdelta);
                 if (std::stof(Mtid) > 0 && std::stof(Mdelta) > 0 )     mMdelta_tidal  = std::stof(Mdelta) / std::stof(Mtid);
                 if (std::stof(Mequdens) > 0 && std::stof(Mdelta) > 0 ) mMdelta_eqdens = std::stof(Mdelta) / std::stof(Mequdens);

                }

                float dm_Glong_wrtHaloCenter   = 0.0;
                float dm_Glat_wrtHaloCenter    = 0.0;
                if (FLAG==2){
                   dm_Glong_wrtHaloCenter   = std::stof(Glong_wrtHaloCenter);
                   dm_Glat_wrtHaloCenter    = std::stof(Glat_wrtHaloCenter);
                }



                std::cout << "dm_glong:  " << mGlong << std::endl;
                std::cout << "dm_glat:   " << mGlat << std::endl;
                std::cout << "dm_dist:   " << mDist << std::endl;
                std::cout << "dm_Rdelta: " << mRdelta << std::endl;
                std::cout << "dm_Rtidal: " << dm_Rtidal << std::endl;
                //std::cout << "dm_rhos:   " << dm_rhos << std::endl;
                //std::cout << "dm_rs:     " << dm_rs << std::endl;
                //mDist  = 0.001;


                if (mRs!=0){
                  mRvirial_rS = 260 / mRs;
                  mRtid_rS   = std::stof(Rtid)/mRs;
                  mReq_rS    = std::stof(Requdens)/mRs;
                  mRdelta_rS = std::stof(Rdelta)/mRs;
                  mDist_rS   = mDist/mRs;
                  mDgal_rS   = std::stof(Dgal)/mRs;
                }

                mDELTA  = (3*std::stof(Mdelta)) / (4* PI * mRhos * std::pow(mRdelta, 3));



                mParam1    = std::stof(DMprof_p1);
                mParam2    = std::stof(DMprof_p2);
                mParam3    = std::stof(DMprof_p3);
                mJfactor   = std::stof(DM_J);
                mJJcont    = std::stof(DM_JJcont);
                if (std::stof(DM_J)>0.)
                  m_log10Jfactor = std::log10(std::stof(DM_J));

                if (std::stof(DM_JJcont)>0.)
                  m_log10JJcont = std::log10(std::stof(DM_JJcont));
                //m_log10JJcont = std::log10(std::stof(DM_JJcont));
               // The concentration , at a given characteristic overdensity , is defined to be
               mCdelta = mRdelta/mRs;
               //where is the Rdelta - radius of the DM halo for which the density equals this overdensity,
               //and rS - is the position where the slope of the DM halo density profile reaches criticality

                if ( rho_DM(std::stof(rhos), std::stof(rs), std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(DM_clump_dist), false) > 0
                  && rho_DM(std::stof(rhos), std::stof(rs), std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(DM_clump_dist), false) < 1e30 )
                  mRhoDM_v0 = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(DM_clump_dist), false);

                if ( rho_DM(std::stof(rhos), std::stof(rs), std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(DM_clump_dist), true) > 0
                  && rho_DM(std::stof(rhos), std::stof(rs), std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(DM_clump_dist), true) < 1e30 )
                  mRhoDM_v1 = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(DM_clump_dist), true);

                if ( rho_DM(std::stof(rhos), std::stof(rs), std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(Dgal), false) > 0
                  && rho_DM(std::stof(rhos), std::stof(rs), std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(Dgal), false) < 1e30 )
                  mRhoDM_v00 = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(Dgal), false);

                if ( rho_DM(std::stof(rhos), std::stof(rs), std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(Dgal), true) > 0
                  && rho_DM(std::stof(rhos), std::stof(rs), std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(Dgal), true) < 1e30 )
                  mRhoDM_v11 = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), (double)std::stof(Dgal), true);


                mRho_dm_r00001 = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-4, true);
                mRho_dm_r0001  = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-3, true);
                mRho_dm_r001   = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-2, true);
                mRho_dm_r01    = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-1, true);
                mRho_dm_r0     = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 0, true);
                mRho_dm_r1     = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1, true);
                mRho_dm_r10    = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 10, true);
                mRho_dm_r100   = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 100, true);

                prf_r00001_rhoDM_Macc -> Fill (mRho_dm_r00001, mMdelta);
                prf_r0001_rhoDM_Macc  -> Fill (mRho_dm_r0001, mMdelta);
                prf_r001_rhoDM_Macc   -> Fill (mRho_dm_r001, mMdelta);
                prf_r01_rhoDM_Macc    -> Fill (mRho_dm_r01, mMdelta);
                prf_r0_rhoDM_Macc     -> Fill (mRho_dm_r0, mMdelta);
                prf_r1_rhoDM_Macc     -> Fill (mRho_dm_r1, mMdelta);
                prf_r10_rhoDM_Macc    -> Fill (mRho_dm_r10, mMdelta);
                prf_r100_rhoDM_Macc   -> Fill (mRho_dm_r100, mMdelta);

                mRho_dm__r00001 = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-4, false);
                mRho_dm__r0001  = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-3, false);
                mRho_dm__r001   = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-2, false);
                mRho_dm__r01    = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-1, false);
                mRho_dm__r0     = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 0, false);
                mRho_dm__r1     = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1, false);
                mRho_dm__r10    = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 10, false);
                mRho_dm__r100   = rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 100, false);

                for ( double sigmaDMOM=1e-2; sigmaDMOM<=1 ; sigmaDMOM+=0.05){// sigmaDMOM = landa with an increment of 0.05
                  mMacc_r00001 = accreted_mass(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-4, true), sigmaDMOM, 13.7, false);
                  mMacc_r0001  = accreted_mass(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-3, true), sigmaDMOM, 13.7, false);
                  mMacc_r001   = accreted_mass(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-2, true), sigmaDMOM, 13.7, false);
                  mMacc_r01    = accreted_mass(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-1, true), sigmaDMOM, 13.7, false);
                  mMacc_r0     = accreted_mass(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 0, true), sigmaDMOM, 13.7, false);
                  mMacc_r1     = accreted_mass(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1, true), sigmaDMOM, 13.7, false);
                  mMacc_r10    = accreted_mass(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 10, true), sigmaDMOM, 13.7, false);
                  mMacc_r100   = accreted_mass(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 100, true), sigmaDMOM, 13.7, false);

                  prf_r_Macc   -> Fill (1e-4, mMacc_r00001);
                  prf_r_Macc   -> Fill (1e-3, mMacc_r0001);
                  prf_r_Macc   -> Fill (1e-2, mMacc_r001);
                  prf_r_Macc   -> Fill (1e-1, mMacc_r01);
                  //prf_r_Macc   -> Fill (0, mMacc_r0);
                  prf_r_Macc   -> Fill (1, mMacc_r1);
                  prf_r_Macc   -> Fill (10, mMacc_r10);
                  prf_r_Macc   -> Fill (100, mMacc_r100);

                  mCc_r00001 = std::log10( accreted_rate(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-4, false), sigmaDMOM, 100, mMdelta, mRdelta,  false) );
                  mCc_r0001  = std::log10( accreted_rate(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-3, false), sigmaDMOM, 100, mMdelta, mRdelta,  false) );
                  mCc_r001   = std::log10( accreted_rate(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-2, false), sigmaDMOM, 100, mMdelta, mRdelta,  false) );
                  mCc_r01    = std::log10( accreted_rate(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1e-1, false), sigmaDMOM, 100, mMdelta, mRdelta, false) );
                  mCc_r0     = std::log10( accreted_rate(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 0, false),   sigmaDMOM, 100, mMdelta, mRdelta,  false) );
                  mCc_r1     = std::log10( accreted_rate(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 1, false),   sigmaDMOM, 100, mMdelta, mRdelta,  false) );
                  mCc_r10    = std::log10( accreted_rate(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 10, false),  sigmaDMOM, 100, mMdelta, mRdelta,  false) );
                  mCc_r100   = std::log10( accreted_rate(rho_DM(std::stof(rhos),   std::stof(rs),  std::stof(DMprof_p1), std::stof(DMprof_p2), std::stof(DMprof_p3), 100, false), sigmaDMOM, 100, mMdelta, mRdelta,  false) );
                }


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

                prf_rrs_dist    -> Fill (mRvirial_rS, mDist);
                prf_rrs_dgal    -> Fill (mRvirial_rS, mDgal);



                prf_alphas_Mdelta   -> Fill (std::asin(mRs/mDist), mMdelta);
                prf_alphas_Mtid     -> Fill (std::asin(mRs/mDist), mMtidal);
                prf_alphas_Rtid     -> Fill (std::asin(mRs/mDist), mRtidal);
                prf_alphas_Mequdens -> Fill (std::asin(mRs/mDist), mMequdens);
                prf_alphas_Requdens -> Fill (std::asin(mRs/mDist), mRequdens);
                prf_alphas_J        -> Fill (std::asin(mRs/mDist), mJfactor);
                prf_alphas_JJcont   -> Fill (std::asin(mRs/mDist), mJJcont);
                prf_alphas_rhos     -> Fill (std::asin(mRs/mDist), mRhos);
                prf_alphas_rhos2    -> Fill (std::asin(mRs/mDist), mRhos2);
                prf_alphas_dist     -> Fill (std::asin(mRs/mDist), mDist);

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

                prf_dist_Rdelta   -> Fill (mDist, mRdelta);
                prf_dist_Mdelta   -> Fill (mDist, mMdelta);
                prf_dist_Mtid     -> Fill (mDist, mMtidal);
                prf_dist_Rtid     -> Fill (mDist, mRtidal);
                prf_dist_Mequdens -> Fill (mDist, mMequdens);
                prf_dist_Requdens -> Fill (mDist, mRequdens);
                prf_dist_J        -> Fill (mDist, mJfactor);
                prf_dist_JJcont   -> Fill (mDist, mJJcont);
                prf_dist_rrs      -> Fill (mDist, mRvirial_rS);
                prf_dist_rs       -> Fill (mDist, mRs);
                prf_dist_rhos     -> Fill (mDist, mRhos);
                prf_dist_rhoDM_v0  -> Fill (mDist, mRhoDM_v0);
                prf_dist_rhoDM_v1  -> Fill (mDist, mRhoDM_v1);
                prf_dist_rhoDM_v00 -> Fill (mDist, mRhoDM_v00);
                prf_dist_rhoDM_v11 -> Fill (mDist, mRhoDM_v11);


                prf_dgal_Rdelta   -> Fill (mDgal, mRdelta);
                prf_dgal_Mdelta   -> Fill (mDgal, mMdelta);
                prf_dgal_Mtid     -> Fill (mDgal, mMtidal);
                prf_dgal_Rtid     -> Fill (mDgal, mRtidal);
                prf_dgal_Mequdens -> Fill (mDgal, mMequdens);
                prf_dgal_Requdens -> Fill (mDgal, mRequdens);
                prf_dgal_J        -> Fill (mDgal, mJfactor);
                prf_dgal_JJcont   -> Fill (mDgal, mJJcont);
                prf_dgal_rrs      -> Fill (mDgal, mRvirial_rS);
                prf_dgal_rs       -> Fill (mDgal, mRs);
                prf_dgal_rhos     -> Fill (mDgal, mRhos);
                prf_dgal_rhoDM_v0  -> Fill (mDgal, mRhoDM_v0);
                prf_dgal_rhoDM_v1  -> Fill (mDgal, mRhoDM_v1);
                prf_dgal_rhoDM_v00 -> Fill (mDgal, mRhoDM_v00);
                prf_dgal_rhoDM_v11 -> Fill (mDgal, mRhoDM_v11);




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


                prf_MequdensDist_rhos     -> Fill (mMequdens/(mDist*mDist), mRhos);
                prf_MequdensDist_rs       -> Fill (mMequdens/(mDist*mDist), mRs);
                prf_MequdensDist_Mdelta   -> Fill (mMequdens/(mDist*mDist), mMdelta);
                prf_MequdensDist_Requdens -> Fill (mMequdens/(mDist*mDist), mRequdens);
                prf_MequdensDist_Rtid     -> Fill (mMequdens/(mDist*mDist), mRtidal);
                prf_MequdensDist_Mtid     -> Fill (mMequdens/(mDist*mDist), mMtidal);
                prf_MequdensDist_J        -> Fill (mMequdens/(mDist*mDist), mJfactor);
                prf_MequdensDist_JJcont   -> Fill (mMequdens/(mDist*mDist), mJJcont);

                prf_MtidDist_rhos     -> Fill (mMtidal/(mDist*mDist), mRhos);
                prf_MtidDist_rs       -> Fill (mMtidal/(mDist*mDist), mRs);
                prf_MtidDist_Mequdens -> Fill (mMtidal/(mDist*mDist), mMequdens);
                prf_MtidDist_Requdens -> Fill (mMtidal/(mDist*mDist), mRequdens);
                prf_MtidDist_Rtid     -> Fill (mMtidal/(mDist*mDist), mRtidal);
                prf_MtidDist_J        -> Fill (mMtidal/(mDist*mDist), mJfactor);
                prf_MtidDist_JJcont   -> Fill (mMtidal/(mDist*mDist), mJJcont);

                prf_MdeltaDist_rhos     -> Fill (mMdelta/(mDist*mDist), mRhos);
                prf_MdeltaDist_rs       -> Fill (mMdelta/(mDist*mDist), mRs);
                prf_MdeltaDist_Mtid     -> Fill (mMdelta/(mDist*mDist), mMtidal);
                prf_MdeltaDist_Rtid     -> Fill (mMdelta/(mDist*mDist), mRtidal);
                prf_MdeltaDist_Mequdens -> Fill (mMdelta/(mDist*mDist), mMequdens);
                prf_MdeltaDist_Requdens -> Fill (mMdelta/(mDist*mDist), mRequdens);
                prf_MdeltaDist_J        -> Fill (mMdelta/(mDist*mDist), mJfactor);
                prf_MdeltaDist_JJcont   -> Fill (mMdelta/(mDist*mDist), mJJcont);



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

                h2d_dgal_Mdelta   -> Fill (mDgal, mMdelta);
                h2d_dgal_Mtid     -> Fill (mDgal, mMtidal);
                h2d_dgal_Rtid     -> Fill (mDgal, mRtidal);
                h2d_dgal_Mequdens -> Fill (mDgal, mMequdens);
                h2d_dgal_Requdens -> Fill (mDgal, mRequdens);
                h2d_dgal_J        -> Fill (mDgal, mJfactor);
                h2d_dgal_JJcont   -> Fill (mDgal, mJJcont);

                h2d_dgal_rs         -> Fill (mDgal, mRs);
                h2d_dgal_dist_rS    -> Fill (mDgal, mDist_rS);
                h2d_dgal_dgal_rS    -> Fill (mDgal, mDgal_rS);
                h2d_dgal_Rvir_rS    -> Fill (mDgal, mRvirial_rS);
                h2d_dgal_Rtid_rS    -> Fill (mDgal, mRtid_rS);
                h2d_dgal_Req_rS     -> Fill (mDgal, mReq_rS);
                h2d_dgal_Rdelta_rS  -> Fill (mDgal, mRdelta_rS);





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
                prf2d_dist_JJcont_Mtidal   -> Fill (mDist, mJJcont,  mMtidal);
                prf2d_dist_JJcont_Mequdens -> Fill (mDist, mJJcont,  mMequdens);
                prf2d_dist_JJcont_Mdelta   -> Fill (mDist, mJJcont,  mMdelta);
                prf2d_dist_JJcont_Rtidal   -> Fill (mDist, mJJcont,  mRtidal);
                prf2d_dist_JJcont_Requdens -> Fill (mDist, mJJcont,  mRequdens);
                prf2d_dist_JJcont_rhos     -> Fill (mDist, mJJcont,  mRhos);
                prf2d_dist_JJcont_Glong    -> Fill (mDist, mJJcont,  mGlong);
                prf2d_dist_JJcont_Glat     -> Fill (mDist, mJJcont,  mGlat);

                prf2d_Glon_Glat_Mdelta    -> Fill (mGlong, mGlat, mMdelta);
                prf2d_Glon_Glat_Mtidal    -> Fill (mGlong, mGlat, mMtidal);
                prf2d_Glon_Glat_Rtidal    -> Fill (mGlong, mGlat, mRtidal);
                prf2d_Glon_Glat_Mequdens  -> Fill (mGlong, mGlat, mMequdens);
                prf2d_Glon_Glat_Requdens  -> Fill (mGlong, mGlat, mRequdens);
                prf2d_Glon_Glat_J         -> Fill (mGlong, mGlat, mJfactor);
                prf2d_Glon_Glat_JJcont    -> Fill (mGlong, mGlat, mJJcont);
                prf2d_Glon_Glat_rhos      -> Fill (mGlong, mGlat, mRhos);
                prf2d_Glon_Glat_rs        -> Fill (mGlong, mGlat, mRs);



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





                prf2d_Mdelta_dist_rs       -> Fill (mMdelta, mDist, mRs);
                prf2d_Mdelta_dist_rhos     -> Fill (mMdelta, mDist, mRhos);
                prf2d_Mdelta_dist_Rtidal   -> Fill (mMdelta, mDist, mRtidal);
                prf2d_Mdelta_dist_Mequdens -> Fill (mMdelta, mDist, mMequdens);
                prf2d_Mdelta_dist_Requdens -> Fill (mMdelta, mDist, mRequdens);
                prf2d_Mdelta_dist_Jfact    -> Fill (mMdelta, mDist, mJfactor);
                prf2d_Mdelta_dist_JJ       -> Fill (mMdelta, mDist, mJJcont);
                prf2d_Mdelta_dist_Glong    -> Fill (mMdelta, mDist, mGlong);
                prf2d_Mdelta_dist_Glat     -> Fill (mMdelta, mDist, mGlat);

                prf2d_Mtidal_dist_rs       -> Fill (mMtidal, mDist, mRs);
                prf2d_Mtidal_dist_rhos     -> Fill (mMtidal, mDist, mRhos);
                prf2d_Mtidal_dist_Mdelta   -> Fill (mMtidal, mDist, mMdelta);
                prf2d_Mtidal_dist_Rtidal   -> Fill (mMtidal, mDist, mRtidal);
                prf2d_Mtidal_dist_Requdens -> Fill (mMtidal, mDist, mRequdens);
                prf2d_Mtidal_dist_Jfact    -> Fill (mMtidal, mDist, mJfactor);
                prf2d_Mtidal_dist_JJ       -> Fill (mMtidal, mDist, mJJcont);
                prf2d_Mtidal_dist_Glong    -> Fill (mMtidal, mDist, mGlong);
                prf2d_Mtidal_dist_Glat     -> Fill (mMtidal, mDist, mGlat);

                prf2d_Mequdens_dist_rs       -> Fill (mMequdens, mDist, mRs);
                prf2d_Mequdens_dist_rhos     -> Fill (mMequdens, mDist, mRhos);
                prf2d_Mequdens_dist_Mdelta   -> Fill (mMequdens, mDist, mMdelta);
                prf2d_Mequdens_dist_Rtidal   -> Fill (mMequdens, mDist, mRtidal);
                prf2d_Mequdens_dist_Requdens -> Fill (mMequdens, mDist, mRequdens);
                prf2d_Mequdens_dist_Jfact    -> Fill (mMequdens, mDist, mJfactor);
                prf2d_Mequdens_dist_JJ       -> Fill (mMequdens, mDist, mJJcont);
                prf2d_Mequdens_dist_Glong    -> Fill (mMequdens, mDist, mGlong);
                prf2d_Mequdens_dist_Glat     -> Fill (mMequdens, mDist, mGlat);


                prf2d_Mdelta_dgal_rs       -> Fill (mMdelta, mDgal, mRs);
                prf2d_Mdelta_dgal_rhos     -> Fill (mMdelta, mDgal, mRhos);
                prf2d_Mdelta_dgal_Rtidal   -> Fill (mMdelta, mDgal, mRtidal);
                prf2d_Mdelta_dgal_Mequdens -> Fill (mMdelta, mDgal, mMequdens);
                prf2d_Mdelta_dgal_Requdens -> Fill (mMdelta, mDgal, mRequdens);
                prf2d_Mdelta_dgal_Jfact    -> Fill (mMdelta, mDgal, mJfactor);
                prf2d_Mdelta_dgal_JJ       -> Fill (mMdelta, mDgal, mJJcont);
                prf2d_Mdelta_dgal_Glong    -> Fill (mMdelta, mDgal, mGlong);
                prf2d_Mdelta_dgal_Glat     -> Fill (mMdelta, mDgal, mGlat);

                prf2d_Mtidal_dgal_rs       -> Fill (mMtidal, mDgal, mRs);
                prf2d_Mtidal_dgal_rhos     -> Fill (mMtidal, mDgal, mRhos);
                prf2d_Mtidal_dgal_Mdelta   -> Fill (mMtidal, mDgal, mMdelta);
                prf2d_Mtidal_dgal_Rtidal   -> Fill (mMtidal, mDgal, mRtidal);
                prf2d_Mtidal_dgal_Requdens -> Fill (mMtidal, mDgal, mRequdens);
                prf2d_Mtidal_dgal_Jfact    -> Fill (mMtidal, mDgal, mJfactor);
                prf2d_Mtidal_dgal_JJ       -> Fill (mMtidal, mDgal, mJJcont);
                prf2d_Mtidal_dgal_Glong    -> Fill (mMtidal, mDgal, mGlong);
                prf2d_Mtidal_dgal_Glat     -> Fill (mMtidal, mDgal, mGlat);

                prf2d_Mequdens_dgal_rs       -> Fill (mMequdens, mDgal, mRs);
                prf2d_Mequdens_dgal_rhos     -> Fill (mMequdens, mDgal, mRhos);
                prf2d_Mequdens_dgal_Mdelta   -> Fill (mMequdens, mDgal, mMdelta);
                prf2d_Mequdens_dgal_Rtidal   -> Fill (mMequdens, mDgal, mRtidal);
                prf2d_Mequdens_dgal_Requdens -> Fill (mMequdens, mDgal, mRequdens);
                prf2d_Mequdens_dgal_Jfact    -> Fill (mMequdens, mDgal, mJfactor);
                prf2d_Mequdens_dgal_JJ       -> Fill (mMequdens, mDgal, mJJcont);
                prf2d_Mequdens_dgal_Glong    -> Fill (mMequdens, mDgal, mGlong);
                prf2d_Mequdens_dgal_Glat     -> Fill (mMequdens, mDgal, mGlat);
        }


      counter++;
      mTree -> Fill(); //Fill the tree. For each event==line
  }



    file_output->Write(); // Save all objects in this file
    file_output->Close(); // Close the file. Note that this is automatically done when you leave
    printf("Output file was closed \n");

    //return 0;
}

