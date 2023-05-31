#import matplotlib.pyplot as plt
#matplotlib.use("Qt5Agg")
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#from matplotlib import cmap

from matplotlib import ticker

import uproot
#import ROOT
#from ROOT import TDirectory
from root_numpy import fill_hist, fill_profile
from root_numpy import root2array, tree2array
from root_numpy import array2tree, array2root



import os
from os import path
import math
import healpy as hp
import numpy as np
import pandas as pd

import randoms as rg

from healpy.newvisufunc import projview, newprojplot
# classic healpy mollweide projections plot with graticule and axis labels
from NPTFit import create_mask as cm # Module for creating masks

import astropy.units as u
import astropy.coordinates as apycoords
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
#from astropy.utils.compat import argparse

import argparse

from astrotools import auger, coord, skymap
from astrotools.skymap import PlotSkyPatch

from gammapy.data import EventList
from gammapy.datasets import MapDataset
from gammapy.irf import PSFMap, EDispKernelMap
from gammapy.maps import Map, MapAxis, WcsGeom
from gammapy.modeling.models import (
    PowerLawSpectralModel,
    PointSpatialModel,
    SkyModel,
    TemplateSpatialModel,
    PowerLawNormSpectralModel,
    Models,
    create_fermi_isotropic_diffuse_model,
)
from gammapy.modeling import Fit



def create(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-Nf', default=1, metavar='Nf', help='Nrandom/Nobject', type=int)
    parser.add_argument('-i', '--inputname', default=None, type=str, metavar='data', help='Input data file')
    parser.add_argument('-o', '--outputname', default=None, type=str, metavar='outfile', help='Output random file')
    parser.add_argument('-rsd', default=True, type=bool, metavar='rsd', help='Include RSD')
    parser.add_argument('-view', action='store_true', help='Visualize the random map')
    parser.add_argument('-nside', default=32, type=int, help='Mask resolution')
    args = parser.parse_args(args)
    name, ext = path.splitext(args.outputname)
    print (args.nside)

    if(ext!='fits'):
        fmt='ascii'
    else:
        fmt='fits'
    rg.make_random(args.nside,args.inputname,args.Nf,outfile=args.outputname,viewmap=args.view,fmt=fmt,rsd=args.rsd)



def cat2hpx(lon, lat, nside, radec=True):
    """
    Convert a catalogue to a HEALPix map of number counts per resolution
    element.

    Parameters
    ----------
    lon, lat : (ndarray, ndarray)
        Coordinates of the sources in degree. If radec=True, assume input is in the icrs
        coordinate system. Otherwise assume input is glon, glat

    nside : int
        HEALPix nside of the target map

    radec : bool
        Switch between R.A./Dec and glon/glat as input coordinate system.

    Return
    ------
    hpx_map : ndarray
        HEALPix map of the catalogue number counts in Galactic coordinates

    """

    npix = hp.nside2npix(nside)

    if radec:
        eq = SkyCoord(lon, lat, 'icrs', unit='deg')
        l, b = eq.galactic.l.value, eq.galactic.b.value
    else:
        l, b = lon, lat

    # conver to theta, phi
    theta = np.radians(90. - b)
    phi = np.radians(l)

    # convert to HEALPix indices
    indices = hp.ang2pix(nside, theta, phi)

    idx, counts = np.unique(indices, return_counts=True)

    # fill the fullsky map
    hpx_map = np.zeros(npix, dtype=int)
    hpx_map[idx] = counts

    return hpx_map


def mollweid_vec(lat_start,long_start,lat_end,long_end,color,width,headwidth,headlength):
    proj = hp.projector.MollweideProj()
    plt.annotate('', xy=(proj.ang2xy(lat_end, long_end,lonlat=True)), xytext=(proj.ang2xy(lat_start, long_start,lonlat=True)),   arrowprops=dict(color=color,width=width,headwidth=headwidth,headlength=headlength))


def view_observed_gsm(self, logged=False, show=False, **kwargs):
    """ View the GSM (Mollweide), with below-horizon area masked. """
    sky = self.observed_sky
    if logged:
        sky = np.log2(sky)

        # Get RA and DEC of zenith
        ra_rad, dec_rad = self.radec_of(0, np.pi / 2)
        ra_deg  = ra_rad / np.pi * 180
        dec_deg = dec_rad / np.pi * 180

        # Apply rotation
        derotate = hp.Rotator(rot=[ra_deg, dec_deg])
        g0, g1 = derotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]

        coordrotate = hp.Rotator(coord=['C', 'G'], inv=True)
        g0, g1 = coordrotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]

        hp.mollview(sky, coord='G', **kwargs)

        if show:
            plt.show()

        return sky

def view(self, idx=0, logged=False):
    """ View generated map using healpy's mollweide projection.
    Parameters
    ----------
    idx: int
        index of map to view. Only required if you generated maps at
        multiple frequencies.
    logged: bool
        Take the log of the data before plotting. Defaults to False.

    """
    if self.generated_map_data is None:
        raise RuntimeError("No GSM map has been generated yet. Run generate() first.")

    if self.generated_map_data.ndim == 2:
        gmap = self.generated_map_data[idx]
        freq = self.generated_map_freqs[idx]
    else:
        gmap = self.generated_map_data
        freq = self.generated_map_freqs

    if logged:
        gmap = np.log2(gmap)

    hp.mollview(gmap, coord='G',
                title='Global Sky Model %s %s' % (str(freq), self.unit))
    plt.show()


def _interp_beam(self,beam_file,pol,chan):
    line0=open(beam_file).readlines()[0]
    if 'dBi' in line0:
        db=True
    else:
        db=False

    print('db='+str(db))
    data=n.loadtxt(beam_file,skiprows=2);
    theta,phi=hp.pix2ang(self.nSide,range(self.nPix))
    theta=n.round(n.degrees(theta)).astype(int)
    phi=n.round(n.degrees(phi)).astype(int)
    if db:
        self.data[pol,chan,:]=10**((data[:,2].squeeze().reshape(360,181))[phi,theta]/10.)
    else:
        self.data[pol,chan,:]=(data[:,2].squeeze().reshape(360,181))[phi,theta]

    self.data[pol,chan,:]/=self.data[pol,chan,:].flatten().max();
    #print('rotatexy='+str(self.rotatexy))
    if self.rotatexz:
        self.data[pol,chan,:]=rotateBeam(self.data[pol,chan,:].flatten(),rot=[0,-90,0])
    if self.rotatexy:
        self.data[pol,chan,:]=rotateBeam(self.data[pol,chan,:].flatten(),rot=[-90,0,90])
    if(self.invert):
        self.data[pol,chan,:]=rotateBeam(self.data[pol,chan,:].flatten(),rot=[0,180,0])
    if DEBUG:
        hp.mollview(self.data[pol,chan,:])
        plt.show()
    self.data[pol,chan,theta>90.]=0.


def plot_sky_projection_healpy_count(Sliced_Halo_data,nside):
    HEALPix_mode = read_data_bool(tag_name = 'HEALPix_Cartesian',file_name = 'parameters/Output_Parameters.xml')
    HEALPix_grat = read_data_bool(tag_name = 'HEALPix_Graticule',file_name = 'parameters/Output_Parameters.xml')
    fdir  = './Output/plots/HEALPix/'
    Sl_n  = len(Sliced_Halo_data)
    Z_max = max(Sliced_Halo_data[Sl_n-1].Z_red[:])
    rc('text',usetex=True)
    for k in range(Sl_n):
        pix     = zeros(12*nside**2)
        n       = len(Sliced_Halo_data[k].RA[:])
        for i in range(n):
           j    = hp.ang2pix(nside,DtoR*(90.0-Sliced_Halo_data[k].DEC[i]),DtoR*Sliced_Halo_data[k].RA[i])
           pix[j] += 1
        clf()
        if (HEALPix_mode):
           hp.cartview(pix)
        else:
           hp.mollview(pix)
        if (HEALPix_grat):
           hp.graticule()
        fname = 'sky_projection_HEALPix_Count_%i_%i.pdf'%(nside,(k+1))
        title(r'sky projection for redshift between %0.2f and %0.2f'%(Sliced_Halo_data[k].z_min,Sliced_Halo_data[k].z_max),fontsize = 20)
        print ('Saving plot', fname)
        savefig(fdir+fname,bbox_inches='tight')
    rc('text',usetex=False)
    close()
    return 0

def plot_map(n, limit, radius='small'):
    """
    Plot the sky area there there are fewer than N stars brighter than limit
    in a radius of either 'small' (0.7 deg) or 'large' (0.9 deg).  For speed, limit
    should be one of the mag limits already evaluated
    limits = [10.13, 10.02, 9.92, 10.38, 10.29, 10.2, 10.0, 10.1, 10.3]
    Makes a mollweide visualization using healpy.
    """
    if limit not in limits:
        raise ValueError("not calculated for limit {}".format(limit))
    if radius == 'small':
        if not os.path.exists('small_{}.npy'.format(limit)):
            make_maps()
        if radius == 'small':
            field = np.load('small_{}.npy'.format(limit))
        else:
            field = np.load('big_{}.npy'.format(limit))
        if len(field) != healpy.nside2npix(NSIDE):
            make_maps()
            if radius == 'small':
                field = np.load('small_{}.npy'.format(limit))
            else:
                field = np.load('big_{}.npy'.format(limit))
        map = (field > n).astype('int')
        hp.mollview(map, xsize=2000)

def column(array, i):
    return [row[i] for row in array]

def read_hpx_maps(fns):
    '''Read in one or more healpix maps and add them together. Must input
    an array of strings even if only inputting a single map.

    Parameters
    ----------
    fns : list of strings
        The filenames for the healpix maps to read in.

    Returns
    -------
    hpx_map: array-like
        A healpix map that is the sum of the Healpix maps in the input files.

    Notes
    -----
    The nside of the output map will be the nside of the file map in the list.
    Every other map will be upgraded or downgraded that that nside value.
    '''

    hpx_map = hp.read_map(fns[0], verbose=False)
    nside   = hp.npix2nside(len(hpx_map))
    for fn_tmp in fns[1:]:
        tmp_map = hp.read_map(fn_tmp, verbose=False)
        hpx_map += hp.ud_grade(tmp_map, nside)

    return hpx_map

l_ene=[]
#Neutron stars coordinates
#Galactic latitude - GB, Galactic latitude or Declination, in degrees, scalar or vector
#Galactic longitude - GL, Galactic longitude or Right Ascension, in degrees, scalar or vector
l_GB_GL_dist = [
[-7.9138,  2.7882,  8.400],
[48.3400,  19.8500,  1.051],
[2.1220,  49.9680,  4.167],
[-27.3200,  65.0300,  10.400],
[-12.0210,  104.9310,  3.159],
[-17.1369,  184.1245,  0.522],
[-1.1870,  168.2747,  1.900],
[-35.0400,  244.5100,  12.100],
[-4.5000,  245.2400,  1.100],
[-4.5000,  245.2400,  1.100],
[-3.8600,  295.7900,  10.000],
[54.2800,  80.8100,  0.964],
[0.9500,  6.5000,  0.730],
[2.877,  9.966,  7.400],
[-2.20605091,  5.83588493,  3.000],
[0.4400,  12.8200,  5.700],
[24.7350,  72.8300,  4.356],
[15.6120,  53.3430,  0.915],
[-1.0100,  37.3400,  7.000],
[0.1500,  41.6000,  7.400],
[0.1900,  45.2500,  5.600],
[-16.8700,  19.9800,  2.004],
[-46.0753,  62.0185,  0.268],
[-32.52909252,  276.3349498,  7.900],
[-43.55935025,  300.4149551,  9.000],
[-43.55935025,  300.4149551,  9.000],
[-11.29080163,  356.8502712,  2.500],
[2.163693898,  327.4196998,  2.000],
[2.163693898,  327.4196998,  6.400],
[37.52303328,  58.14903389,  6.600],
[2.173484363,  347.754429,  1.800],
[2.173484363,  347.754429,  1.800],
[0.335517538,  292.0903507,  5.700],
[-0.3539,  351.4972,  7.500],
[-43.804,  303.514,  59.700],
[-32.52909252,  276.3349498,  7.900],
[3.929852619,  263.0582983,  1.900],
[3.929852619,  263.0582983,  1.900],
[-0.8505,  330.9263,  5.800],
[37.52303328,  58.14903389,  6.600],
[75.4140,  311.3100,  0.709],
[15.9600,  350.9760,  1.800],
[-0.5970,  21.5870,  4.685],
[3.5010,  55.7770,  0.300],
[2.1554,  78.4883,  6.100],
[-11.31636463,  87.32819642,  8.000],
[-19.8109,  279.9777,  7.100],
[-0.0470,  359.9440,  8.300],
[-0.1750,  359.7880,  8.109],
[-0.0200,  0.1500,  8.217],
[-0.2330,  0.1260,  8.149],
[-0.0700,  6.8400,  4.000],
[16.8060,  44.6420,  3.030],
[0.3870,  12.9040,  4.492],
[-1.0010,  16.8050,  5.000],
[3.65256,  1.07299,  8.000],
[40.91351866,  58.99968688,  7.100],
[40.91268749,  58.99528387,  7.100],
[40.9127738,  58.97150791,  7.100],
[40.90884273,  59.00526155,  7.100],
[40.91029222,  59.00755123,  7.100],
[40.90293961,  59.02380103,  7.100],
[0.319238,  344.369163,  6.400],
[-17.21600527,  358.5988662,  0.167],
[2.830882823,  295.7666832,  8.178],
[-2.092058679,  31.07638514,  10.000],
[35.38,  243.652,  6.260],
[-9.33975534,  281.8354681,  10.000],
[6.7700,  20.7900,  7.800],
[3.0600,  42.2900,  1.200],
[-4.697,  59.197,  1.400],
[-30.0400,  169.9900,  1.300],
[-36.7736,  183.3368,  2.100],
[-41.9600,  253.3900,  0.157],
[-2.0100,  200.5700,  0.420],
[29.5990,  149.7300,  1.136],
[21.0900,  202.7300,  1.110],
[-5.7440,  283.6680,  4.04	],
[50.8600,  160.3500,  0.700],
[51.1000,  231.7900,  0.645],
[45.7800,  243.4900,  1.370],
[12.2540,  280.8510,  0.340],
[13.8000,  298.9700,  1.613],
[16.451,  344.09,  1.887],
[20.1900,  352.6400,  0.700],
[25.2200,  28.7500,  1.311],
[17.7400,  27.7200,  1.471],
[-11.9669,  338.1647,  2.400],
[21.6410,  37.8850,  1.667],
[1.6900,  3.8400,  6.900],
[1.663,  6.299,  3.232],
[-2.2060,  5.8360,  3.000],
[5.367,  44.875,  2.083],
[-19.6000,  359.7300,  1.140],
[1.795,  46.564,  1.496],
[-25.7300,  336.5250,  4.000],
[-9.1200,  30.0300,  1.111],
[4.7100,  69.2900,  6.940],
[2.5536,  66.8583,  6.500],
[-1.1687,  61.0975,  7.269],
[-8.675,  60.522,  2.162],
[-6.6200,  64.7500,  1.163],
[38.2710,  41.0510,  1.515],
[-16.0260,  48.6210,  1.399],
[-15.3100,  61.9200,  1.600],
[-14.0050,  103.3950,  0.863],
[-42.3600,  91.3610,  1.667],
[-3.926,  77.832,  5.100],
[1.302,  86.861,  4.120],
[-42.0800,  47.7800,  0.500],
[-43.006,  72.991,  0.971],
[-36.19939,  46.48299,  8.500],
[46.8060,  3.8590,  7.500],
[-44.9020,  305.8960,  4.690],
[-44.9030,  305.9000,  4.900],
[3.80168,  7.72912,  8.200],
[-5.00873502,  353.5321094,  8.200],
[0.6110,  8.3820,  0.760],
[-5.58068,  7.79821,  5.600],
[-47.44257,  143.961909,  3.900],
[1.6700,  3.8100,  6.900],
[1.6700,  3.8100,  6.900],
[1.6700,  3.8100,  6.900],
[1.6700,  3.8100,  6.900],
[1.6700,  3.8100,  6.900],
[1.6700,  3.8100,  6.900],
[1.6700,  3.8100,  6.900],

[1.94,  326.35,  2.9],#J1537-5312 \citep{Cameron:2020pin}
[-2.05, 325.08,  1.9],#J1547-5709 \citep{Cameron:2020pin}
[2.79,  335.82,  2.4],#J1618-4624 \citep{Cameron:2020pin}
[2.92,  357.04,  3.8],#J1727-2951 \citep{Cameron:2020pin}

[-0.34, 3.76,  10.3],#J1755-25 \citep{Ng:2015zza}
[-1.13, 302.20, 5.6],#J1244-6359 \citep{Ng:2015zza}
[-4.02, 291.42, 4.5],#J1101-6424 \citep{Ng:2015zza}
[20.19, 352.64, 1.2],#J1614-2230 \citep{Ng:2015zza}

[37.895, 128.289, 0.903], #J1125+7819 128.289   37.895    0.903 \citep{Lynch:2018zxo}
[31.763, 113.840, 3.035], #J1641+8049 113.840   31.763    3.035 \citep{Lynch:2018zxo}
[20.030, 97.946, 3.380]  #J1938+6604 97.946    20.030    3.380 \citep{Lynch:2018zxo}

]

RA_DEC_deg_dist_kpc=[
[275.9186959, -30.361146, 8.400],
[234.2915072, 11.93206496, 1.051],
[288.8666643, 16.10760744, 4.167],
[322.5050175, 12.1772803, 10.400],
[346.48268, 47.12926, 3.159],
[73.4392237, 15.9892517, 0.522],
[77.3824504, 38.02169, 1.900],
[78.5278863, -40.0469147, 12.100],
[114.4635351, -30.66130953, 1.100],
[114.4635351, -30.66130953, 1.100],
[175.279225, -65.7553092, 10.000],
[229.5699962, 49.07618089, 0.964],
[269.1943076, -22.866486, 0.730],
[269.2657683, -18.9009378, 7.400],
[271.8369634, -25.00053194, 3.000],
[272.979308, -17.61047, 5.700],
[274.1497265, 45.1760727, 4.356],
[277.394445, 24.9383869, 0.915],
[285.7741384, 3.455335864, 7.000],
[286.70358, 7.77386, 7.400],
[288.3710592, 11.034928, 5.600],
[292.623815, -18.862853, 2.004],
[335.524871, -1.62103475, 0.268],
[83.20648168, -66.37033409, 7.900],
[19.27144054, -73.44333745, 9.000],
[19.27144054, -73.44333745, 9.000],
[276.4450771, -37.10514704, 2.500],
[235.5973471, -52.38599372, 2.000],
[235.5973471, -52.38599372, 6.400],
[254.4575459, 35.34235738, 6.600],
[255.9865524, -37.84414259, 1.800],
[255.9865524, -37.84414259, 1.800],
[170.3128836, -60.62378618, 5.700],
[261.297496, -36.282665, 7.500],
[11.3965, -73.3175, 59.700],
[83.20648168, -66.37033409, 7.900],
[135.528587, -40.55469421, 1.900],
[135.528587, -40.55469421, 1.900],
[243.1793333, -52.4232222, 5.800],
[254.4575459, 35.34235738, 6.600],
[195.0149029, 12.68235336, 0.709],
[245.9092575, -26.5316025, 1.800],
[278.170278, -10.3591, 4.685],
[290.436729, 21.883958, 0.300],
[305.5, 40.7, 6.100],
[326.1714768, 38.32140738, 8.000],
[117.14084, -67.75153, 7.100],
[266.417359, -29.008304, 8.300],
[266.4492935, -29.20855, 8.109],
[266.513989, 28.8370514, 8.217],
[266.7077283, 28.9497194, 8.149],
[270.332562, -23.079066, 4.000],
[272.65533, 17.743717, 3.030],
[273.0660802, -17.5605197, 4.492],
[276.2623, -14.7816, 5.000],
[263.5560833, -26.0885, 8.000],
[250.4202925, 36.45416078, 7.100],
[250.4208645, 36.45076217, 7.100],
[250.4182977, 36.43291333, 7.100],
[250.4266468, 36.45783392, 7.100],
[250.4250921, 36.45971322, 7.100],
[250.4358575, 36.47111206, 7.100],
[255.2036833, -41.65596111, 6.400],
[284.149375, -37.91019444, 0.167],
[178.755525, -59.25375978, 8.178],
[283.8766912, -2.60464935, 10.000],
[148.0346725, -6.1231917, 6.260],
[140.6444842, -63.29482105, 10.000],
[271.207898, -7.59019, 7.800],
[284.4016262, 9.721442158, 1.200],
[299.9032078, 20.80420061, 1.400],
[54.4326079, 17.2541189, 1.300],
[57.18182917, 4.53651611, 2.100],
[69.31623406, -47.25253075, 0.157],
[95.34214317, 10.0440931, 0.420],
[115.1908032, 66.34264435, 1.136],
[117.7881472, 18.1273573, 1.110],
[148.8, -61.83333, 4.04	],
[153.1393269, 53.11729541, 0.700],
[155.741643, 10.031299, 0.645],
[155.9486967, 0.64467931, 1.370],
[161.4591034, -45.16502896, 0.340],
[186.994683, -48.8952058, 1.613],
[240.2162629, -30.8970535, 1.887],
[243.6521148, -22.5086934, 0.700],
[258.4563908, 7.793744379, 1.311],
[264.72486, 3.553018519, 1.471],
[265.1858792, -53.67821389, 2.400],
[265.3797697, 13.86225441, 1.667],
[267.02, -24.77917, 6.900],
[268.41603, -22.6783, 3.232],
[271.8369634, -25.000532, 3.000],
[283.4888259, 13.06223574, 2.083],
[287.4476345, -37.73738437, 1.140],
[287.5404229, 12.94040385, 1.496],
[287.9281484, -59.97413969, 4.000],
[289.7001363, -6.70969457, 1.111],
[296.6047079, 34.28741392, 6.940],
[297.3734894, 31.10105723, 6.500],
[297.6877653, 24.24915664, 7.269],
[304.2393479, 19.7976634, 2.162],
[304.8831208, 24.42091772, 1.163],
[250.0697711, 22.40245113, 1.515],
[304.3446051, 6.05154554, 1.399],
[310.8370056, 17.19136111, 1.600],
[345.6957433, 44.70613348, 0.863],
[349.2884849, 14.65868362, 1.667],
[311.2562713, 36.55039017, 5.100],
[313.4692835, 46.84769947, 4.120],
[326.460248, -7.83847452, 0.500],
[338.5961379, 6.191301758, 0.971],
[320.8105833, -5.79802778, 8.500],
[229.631075, 2.087631, 7.500],
[6.02793, -72.06855742, 4.690],
[14.4333, -72.0219, 4.900],
[267.2194583, -20.35958333, 8.200],
[267.5575067, -37.05304167, 8.200],
[270.5222316, -21.4010136, 0.760],
[276.1370417, -24.86983333, 5.600],
[27.28660833, 13.05901389, 3.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],

[234.40706108,  -53.20696028,  2.9],#J1537-5312 \citep{Cameron:2020pin}
[236.85052,     -57.15488053,  1.9],#J1547-5709 \citep{Cameron:2020pin}
[244.71989913,  -46.40970833,  2.4],#J1618-4624 \citep{Cameron:2020pin}
[261.751675,    -29.86133333,  3.8],#J1727-2951 \citep{Cameron:2020pin}

[268.775,      -25.88333333,  10.3],#J1755-25 \citep{Ng:2015zza}
[191.19872083, -63.9965,      5.6],#J1244-6359 \citep{Ng:2015zza}
[165.40496792, -64.41092556,  4.5],#J1101-6424 \citep{Ng:2015zza}
[243.65210458, -22.50863361,  1.2],#J1614-2230 \citep{Ng:2015zza}

[171.4994027, 78.33019874, 0.903], #J1125+7819 128.289   37.895    0.903 \citep{Lynch:2018zxo}
[250.3368257, 80.8313651, 3.035], #J1641+8049 113.840   31.763    3.035 \citep{Lynch:2018zxo}
[294.7371730, 80.8313651, 3.380]  #J1938+6604 97.946    20.030    3.380 \citep{Lynch:2018zxo}
]

NS_and_compStar_mass_and_unc=[
[1.58, 0.06, 0.06, 0, 0, 0],
[1.3332, 0.001, 0.001, 1.3452, 0.001, 0.001],
[1.4398, 0.002, 0.002, 1.3886, 0.002, 0.002],
[1.358, 0.01, 0.01, 1.354, 0.01, 0.01],
[1.38, 0.06, 0.1, 1.4, 0.1, 0.1],
[1.559, 0.005, 0.005, 1.174, 0.004, 0.004],
[1.34, 0.08, 0.08, 1.46, 0.08, 0.08],
[1.25, 0.05, 0.06, 1.22, 0.06, 0.05],
[1.3381, 0.0007, 0.0007, 1.2489, 0.0007, 0.0007],
[1.2489, 0.0007, 0.0007, 1.337, 0.005, 0.005],
[1.27, 0.01, 0.01, 1.02, 0.01, 0.01],
[1.56, 0.13, 0.44, 1.05, 0.45, 0.11],
[1.341, 0.007, 0.007, 1.23, 0.007, 0.007],
[1.3384, 0.09, 0.09, 1.3946, 0.09, 0.09],
[1.3655, 0.0021, 0.0021, 1.2068, 0.0022, 0.0016],
[1.5, 0.12, 0.4, 1.06, 0.45, 0.1],
[1.84, 0.11, 0.11, 0.193, 0.012, 0.012],
[1.306, 0.007, 0.007, 1.299, 0.007, 0.007],
[1.666, 0.01, 0.012, 1.033, 0.011, 0.008],
[1.291, 0.011, 0.011, 1.322, 0.011, 0.011],
[1.65, 0.05, 0.05, 1.24, 0.05, 0.05],
[1.29, 0.1, 0.1, 1.3, 0.1, 0.1],
[1.831, 0.01, 0.01, 1.3194, 0.04, 0.04],
[1.57, 0.11, 0.11, 1.29, 0.05, 0.05],
[1.21, 0.12, 0.12, 1.04, 0.09, 0.09],
[1.05, 0.09, 0.09, 15.5, 1.5, 1.5],
[0.97, 0.24, 0.24, 0.33, 0.05, 0.05],
[1.96, 0.36, 0.36, 0, 0, 0],
[1.02, 0.17, 0.17, 16.4, 5.2, 4],
[1.5, 0.3, 0.3, 2.3, 0.3, 0.3],
[2.44, 0.27, 0.27, 58, 11, 11],
[1.96, 0.19, 0.19, 0, 0, 0],
[1.57, 0.16, 0.16, 19.7, 4.3, 4.3],
[1.46, 0.38, 0.38, 13.6, 1.6, 1.6],
[1.58, 0.34, 0.34, 8.8, 1.8, 1.8],
[1.31, 0.14, 0.14, 15.6, 1.8, 1.8],
[1.88, 0.13, 0.13, 23.1, 0.2, 0.2],
[1.77, 0.08, 0.08, 0.35, 0.1, 0.1],
[1.74, 0.14, 0.14, 0, 0, 0],
[1.073, 0.36, 0.36, 0, 0, 0],
[1.4, 0.1, 0.1, 0, 0, 0],
[1.35, 0.1, 0.1, 0, 0, 0],
[1.4, 0.1, 0.1, 0, 0, 0],
[1.4, 0.1, 0.1, 0, 0, 0],
[0.832, 1.19, 0.0551, 0, 0, 0],
[1.78, 0.23, 0.23, 0.6, 0.13, 0.13],
[1.77, 0.4, 0.7, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[2.13, 0.04, 1.3, 0.65, 0.001, 0.001],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.3, 0.5, 0.5, 0, 0, 0],
[1.48, 0.21, 0.64, 0, 0, 0],
[1.48, 0.21, 0.64, 0, 0, 0],
[1.48, 0.21, 0.64, 0.1876, 0.1, 0.1],
[1.48, 0.21, 0.64, 0.2085, 0.1, 0.1],
[1.48, 0.21, 0.64, 0.0225, 0.1, 0.1],
[1.48, 0.21, 0.64, 0.1517, 0.1, 0.1],
[1.74, 0.3, 0.3, 16, 2, 2],
[1.7, 0.3, 0.3, 0, 0, 0],
[1.43, 0.26, 0.61, 0, 0, 0],
[1.41, 0.24, 0.24, 0, 0, 0],
[2.35, 0.17, 0.17, 0.032, 0.002, 0.002],
[1.44, 0.1, 0.1, 0.35, 0.03, 0.03],
[1.26, 0.08, 0.17, 0.36, 0.67, 0.15],
[1.37, 0.13, 0.1, 0.244, 0.014, 0.012],
[2.4, 0.12, 0.12, 0.035, 0.002, 0.002],
[1.4378, 0.0013, 0.0013, 0.19751, 0.00015, 0.00015],
[2.01, 0.04, 0.04, 0.172, 0.003, 0.003],
[1.76, 0.15, 0.15, 0.254, 0.014, 0.014],
[1.7, 0.1, 0.17, 0.76, 0.28, 0.07],
[2.08, 0.1, 0.1, 0.2, 0.1, 0.1],
[1.26, 0.14, 0.14, 0.12, 0.02, 0.02],
[1.71, 0.02, 0.02, 0.254, 0.002, 0.002],
[1.83, 0.11, 0.11, 0.16, 0.02, 0.02],
[1.7, 0.3, 0.3, 0.92, 0.05, 0.05],
[1.761, 0.16, 0.16, 0.241, 0.02, 0.02],
[1.35, 0.1, 0.1, 0.19, 0.1, 0.1],
[0.86, 1.52, 1.52, 0.167, 0.3, 0.3],
[2.5, 0.9, 0.7, 0.34, 0.09, 0.07],
[1.908, 0.016, 0.016, 0.493, 0.003, 0.003],
[1.35, 0.07, 0.07, 0.292, 0.011, 0.011],
[1.47, 0.07, 0.06, 0.181, 0.007, 0.005],
[1.53, 0.19, 0.19, 0.181, 0.1, 0.1],
[1.14, 0.43, 0.25, 0.22, 0.05, 0.04],
[1.91, 0.02, 0.1, 0.26, 0.1, 0.1],
[1.25, 0.2, 0.2, 0.49, 0.1, 0.1],
[1.34, 0.1, 0.1, 0.0092, 0.1, 0.1],
[1.4, 0.7, 0.7, 0.33, 0.37, 0.37],
[1.48, 0.03, 0.03, 0.208, 0.002, 0.002],
[1.6, 0.6, 0.6, 0.3, 0.34, 0.34],
[1.33, 0.11, 0.11, 0.18, 0.018, 0.018],
[1.18, 0.01, 0.09, 0.219, 0.01, 0.01],
[1.832, 0.0029, 0.0029, 0.2659, 0.003, 0.003],
[1.34, 0.17, 0.15, 0.81, 0.06, 0.05],
[1.496, 0.023, 0.023, 0.28, 0.005, 0.005],
[1, 0.5, 0.5, 0.43, 0.5, 0.5],
[1.33, 0.1, 0.1, 0.33, 0.1, 0.1],
[4.4, 2.9, 2, 0.6, 0.4, 0.2],
[2.4, 3.4, 1.4, 0.32, 0.44, 0.16],
[1.41, 0.2, 0.18, 0.175, 0.016, 0.015],
[5.3, 3.2, 3.6, 2.3, 1.7, 1.3],
[4.7, 3.4, 2.8, 0.7, 0.5, 0.4],
[1.33, 0.3, 0.28, 0.94, 0.14, 0.13],
[1.4, 0.21, 0.18, 0.86, 0.07, 0.06],
[1.3, 0.4, 0.5, 0.83, 0.06, 0.06],
[1.353, 0.014, 0.017, 0.298, 0.15, 0.12],
[1.46, 0.3, 0.39, 0.53, 0.28, 0.39],
[2.08, 0.19, 0.19, 0.21, 0.1, 0.1],
[1.61, 0.04, 0.04, 0.18, 0.086, 0.016],
[1.41, 0.1, 0.1, 0.15, 0.1, 0.1],
[2.92, 0.2, 0.2, 0.12, 0.1, 0.1],
[1.97, 0.015, 0.015, 0.54, 0.1, 0.1],
[1.24, 0.11, 0.11, 0.78, 0.04, 0.04],
[1.616, 0.007, 0.007, 0.27, 0.1, 0.1],
[1.34, 0.08, 0.08, 0.175, 0.01, 0.01],
[1.32, 0.1, 0.1, 0.49, 0.1, 0.1],
[1.874, 0.32, 0.068, 0.25, 0.1, 0.1],
[1.728, 0.066, 0.136, 0.35, 0.1, 0.1],
[1.73, 0.1, 0.1, 0.39, 0.1, 0.1],
[1.69, 0.1, 0.1, 0.25, 0.1, 0.1],
[1.6, 0.1, 0.1, 0.25, 0.1, 0.1],
[1.48, 0.1, 0.1, 0.22, 0.1, 0.1],

[1.4, 0.1, 0.1, 0.1339705,  0.0000003, 0.0000003],#J1537-5312 \citep{Cameron:2020pin}
[1.4, 0.1, 0.1, 0.20435198, 0.00000017, 0.00000017],#J1547-5709 \citep{Cameron:2020pin}
[1.4, 0.1, 0.1, 0.7044464,  0.0000006, 0.0000006],#J1618-4624 \citep{Cameron:2020pin}
[1.4, 0.1, 0.1, 0.01601,    0.00003, 0.00003],#J1727-2951 \citep{Cameron:2020pin}

[1.35, 0.1, 0.1, 0.48,  0.1, 0.1],#J1755-25 \citep{Ng:2015zza}
[1.35, 0.1, 0.1, 0.68,  0.1, 0.1],#J1244-6359 \citep{Ng:2015zza}
[1.35, 0.1, 0.1, 0.57,  0.1, 0.1],#J1101-6424 \citep{Ng:2015zza}
[1.35, 0.1, 0.1, 0.5,   0.006, .006],#J1614-2230 \citep{Ng:2015zza}

[1.35, 0.1, 0.1, 0.29, 0.1, 0.1], #J1125+7819 128.289   37.895    0.903 \citep{Lynch:2018zxo}
[1.35, 0.1, 0.1, 0.04, 0.1, 0.1], #J1641+8049 113.840   31.763    3.035 \citep{Lynch:2018zxo}
[1.35, 0.1, 0.1, 0.87, 0.1, 0.1]  #J1938+6604 97.946    20.030    3.380 \citep{Lynch:2018zxo}

]

NS_age_Gyr = [
[0.0,   '4U 1820-30'],
[0.0,   'B1534+12'],
[0.0,   'B1913+16'],
[0.0,   'B2127+11C'],
[0.0,   'B2303+46'],
[0.0,   'J0453+1559'],
[0.15291,   'J0509+3801'],
[0.0,   'J0514-4002A / NGC 1851'],
[0.0,   'J0737-3039A'],
[0.0,   'J0737-3039B'],
[0.0,   'J1141-6545'],
[0.0,   'J1518+4904'],
[0.4435,   'J1756-2251'],
[0.130,   'J1757-1854'],
[0.0,   'J1807-2500B'],
[2.341E-3,   'J1811-1736'],
[0.0,   'J1816+4510'],
[13.0,   'J1829+2456'],
[0.0,   'J1903+0327'],
[0.0,   'J1906+0746'],
[0.0,   'J1913+1102'],
[0.163,   'J1930-1852'],
[30.6,   'J2222-0137'],
[0.0,   'LMC X-4'],
[0.0,   'Sk 160/SMC X-1'],
[0.0,   'Sk 160/SMC X-1'],
[0.0,   '2A 1822-371/V* V691 CrA'],
[0.0,   '4U 1538-522/V* QV Nor'],
[0.0,   '4U 1538-522/V* QV Nor'],
[0.0,   '4U 1656+35/Her X-1'],
[0.0,   '4U 1700-377/HD 153919'],
[0.0,   '4U 1700-377/HD 153919'],
[0.0,   'Cen X-3/V* V779 Cen'],
[0.0,   'EXO 1722-363'],
[0.0,   'J0045-7319'],
[0.0,   'LMC X-4'],
[0.0,   'Vela X-1/V* GP Vel'],
[0.0,   'Vela X-1/V* GP Vel'],
[0.0,   '4U 1608-522'],
[0.0,   '4U 1656+35/Her X-1'],
[0.0,   'B1257+12'],
[0.0,   'B1620-26'],
[0.0,   'B1829-10'],
[0.0,   'B1919+21'],
[0.0,   'Cyg X-7'],
[0.0,   'Cyg X-2/V* V1341 Cyg'],
[0.0,   'EXO 0748-676'],
[0.0,   'J1745-2900'],
[0.0,   'J1745-2912'],
[0.0,   'J1746-2849'],
[0.0,   'J1746-2856'],
[0.0,   'J1801-2304'],
[0.0,   'J1810+1744'],
[0.0,   'J1812-1733'],
[0.0,   'J1825-1446'],
[0.0,   'KS 1731-260'],
[0.0,   'M13 A'],
[0.0,   'M13 C'],
[0.0,   'M13 B'],
[0.0,   'M13 D'],
[0.0,   'M13 E'],
[0.0,   'M13 F'],
[0.0,   'OAO 1657-415'],
[0.0,   'RX J1856–3754'],
[0.0,   'w Cen'],
[0.0,   'XTE J1855-026'],
[0.0,   'J0952-0607'],
[0.0,   '2S 0921-630/V* V395 Car'],
[0.0,   'B1802-07'],
[0.0,   'B1855+09'],
[0.0,   'B1957+20'],
[0.0,   'J0337+1715'],
[0.0,   'J0348+0432'],
[6.6,   'J0437-4715'],
[0.0,   'J0621+1002'],
[3.746,   'J0740+6620'],
[8.0,   'J0751+1807'],
[0.0,   'J0955-6150'],
[0.0,   'J1012+5307'],
[0.0,   'J1022+1001'],
[0.0,   'J1023+0038'],
[0.0,   'J1045-4509'],
[2.4,   'J1227-4853'],
[0.0,   'J1600-3053'],
[1.2,   'J1614-2230'],
[0.0,   'J1713+0747'],
[2.25,   'J1738+0333'],
[0.0,   'J1740-5340 (NGC 6397 (20))'],
[0.0,   'J1741+1351'],
[0.0,   'J1748-2446I'],
[0.0,   'J1753-2240'],
[0.0,   'J1807-2459B'],
[0.0,   'J1853+1303'],
[0.0,   'J1909-3744'],
[0.0,   'J1910+1256'],
[0.0,   'J1910-5959A'],
[0.0,   'J1918-0642'],
[0.0,   'J1946+3417'],
[2.2,   'J1949+3106'],
[3.4,   'J1950+2414'],
[0.0,   'J2016+1948'],
[0.0,   'J2019+2425'],
[0.0,   'J1640-2224'],
[0.0,   'J2017-0603'],
[0.0,   'J2043+1711'],
[0.0,   'J2302+4442'],
[0.0,   'J2317+1439'],
[0.85,   'J2045+3633'],
[1.15,   'J2053+4650'],
[10.4,   'J2145-0750'],
[8.8,   'J2234+0611'],
[0.0,   'XTE J2123-058/V* LZ Aqr '],
[0.0,   'B1516+02B/J1518+0204B (M5)'],
[0.0,   'J0024-7204H/B0021-72H (47 Tucanae)'],
[0.0,   'J0024-7204I/B0021-72I (47 Tucanae)'],
[0.0,   'J1748-2021B / NGC 6440B'],
[0.0,   'J1750-37A / NGC 6441A'],
[0.0,   'J1802-2124'],
[0.0,   'J1824-2452C / M28'],
[0.0,   'J1911-5958A / NGC 6752'],
[0.0,   'Ter 5ai / J1748-2446ai'],
[0.0,   'Ter 5I / J1748-2446I '],
[0.0,   'Ter 5J / J1748-2446J '],
[0.0,   'Ter 5U / J1748-2446U'],
[0.0,   'Ter 5W / J1748-2446W'],
[0.0,   'Ter 5X / J1748-2446X'],
[0.0,   'Ter 5Z / J1748-2446Z'],

[6.9,   'J1537-5312'],
[9.1,   'J1547-5709'],
[1.4,   'J1618-4624'],
[30.2,   'J1727-2951'],

[0.0,   'J1755-25'],
[0.520,   'J1244-6359'],
[44.0,   'J1101-6424'],
[5.2,   'J1614-2230'],

[0.0,   'J1125+7819'],
[0.0,   'J1641+8049'],
[0.0,   'J1938+6604']
]


NS_J2000_name = [
['NS-NS binaries',  '4U 1820-30'],
['NS-NS binaries',  'B1534+12'],
['NS-NS binaries',  'B1913+16'],
['NS-NS binaries',  'B2127+11C'],
['NS-NS binaries',  'B2303+46'],
['NS-NS binaries',  'J0453+1559'],
['NS-NS binaries',  'J0509+3801'],
['NS-NS binaries',  'J0514-4002A / NGC 1851'],
['NS-NS binaries',  'J0737-3039A'],
['NS-NS binaries',  'J0737-3039B'],
['NS-NS binaries',  'J1141-6545'],
['NS-NS binaries',  'J1518+4904'],
['NS-NS binaries',  'J1756-2251'],
['NS-NS binaries',  'J1757-1854'],
['NS-NS binaries',  'J1807-2500B'],
['NS-NS binaries',  'J1811-1736'],
['NS-NS binaries',  'J1816+4510'],
['NS-NS binaries',  'J1829+2456'],
['NS-NS binaries',  'J1903+0327'],
['NS-NS binaries',  'J1906+0746'],
['NS-NS binaries',  'J1913+1102'],
['NS-NS binaries',  'J1930-1852'],
['NS-NS binaries',  'J2222-0137'],
['NS-NS binaries',  'LMC X-4'],
['NS-NS binaries',  'Sk 160/SMC X-1'],
['NS in X-ray binaries',  'Sk 160/SMC X-1'],
['NS in X-ray binaries',  '2A 1822-371/V* V691 CrA'],
['NS in X-ray binaries',  '4U 1538-522/V* QV Nor'],
['NS in X-ray binaries',  '4U 1538-522/V* QV Nor'],
['NS in X-ray binaries',  '4U 1656+35/Her X-1'],
['NS in X-ray binaries',  '4U 1700-377/HD 153919'],
['NS in X-ray binaries',  '4U 1700-377/HD 153919'],
['NS in X-ray binaries',  'Cen X-3/V* V779 Cen'],
['NS in X-ray binaries',  'EXO 1722-363'],
['NS in X-ray binaries',  'J0045-7319'],
['NS in X-ray binaries',  'LMC X-4'],
['NS in X-ray binaries',  'Vela X-1/V* GP Vel'],
['WD-NS binaries',  'Vela X-1/V* GP Vel'],
['radio millisecond pulsars',  '4U 1608-522'],
['radio millisecond pulsars',  '4U 1656+35/Her X-1'],
['radio millisecond pulsars',  'B1257+12'],
['radio millisecond pulsars',  'B1620-26'],
['radio millisecond pulsars',  'B1829-10'],
['radio millisecond pulsars',  'B1919+21'],
['radio millisecond pulsars',  'Cyg X-7'],
['WD-NS binaries',  'Cyg X-2/V* V1341 Cyg'],
['radio millisecond pulsars',  'EXO 0748-676'],
['radio millisecond pulsars',  'J1745-2900'],
['radio millisecond pulsars',  'J1745-2912'],
['radio millisecond pulsars',  'J1746-2849'],
['radio millisecond pulsars',  'J1746-2856'],
['radio millisecond pulsars',  'J1801-2304'],
['radio millisecond pulsars',  'J1810+1744'],
['radio millisecond pulsars',  'J1812-1733'],
['radio millisecond pulsars',  'J1825-1446'],
['radio millisecond pulsars',  'KS 1731-260'],
['radio millisecond pulsars',  'M13 A'],
['radio millisecond pulsars',  'M13 C'],
['NS in X-ray binaries',  'M13 B'],
['NS in X-ray binaries',  'M13 D'],
['NS in X-ray binaries',  'M13 E'],
['NS in X-ray binaries',  'M13 F'],
['radio millisecond pulsars',  'OAO 1657-415'],
['radio millisecond pulsars',  'RX J1856–3754'],
['radio millisecond pulsars',  'w Cen'],
['radio millisecond pulsars',  'XTE J1855-026'],
['radio millisecond pulsars',  'J0952-0607'],
['WD-NS binaries',  '2S 0921-630/V* V395 Car'],
['WD-NS binaries',  'B1802-07'],
['WD-NS binaries',  'B1855+09'],
['WD-NS binaries',  'B1957+20'],
['WD-NS binaries',  'J0337+1715'],
['WD-NS binaries',  'J0348+0432'],
['WD-NS binaries',  'J0437-4715'],
['WD-NS binaries',  'J0621+1002'],
['WD-NS binaries',  'J0740+6620'],
['WD-NS binaries',  'J0751+1807'],
['WD-NS binaries',  'J0955-6150'],
['WD-NS binaries',  'J1012+5307'],
['WD-NS binaries',  'J1022+1001'],
['WD-NS binaries',  'J1023+0038'],
['WD-NS binaries',  'J1045-4509'],
['WD-NS binaries',  'J1227-4853'],
['WD-NS binaries',  'J1600-3053'],
['WD-NS binaries',  'J1614-2230'],
['WD-NS binaries',  'J1713+0747'],
['WD-NS binaries',  'J1738+0333'],
['WD-NS binaries',  'J1740-5340 (NGC 6397 (20))'],
['WD-NS binaries',  'J1741+1351'],
['WD-NS binaries',  'J1748-2446I'],
['WD-NS binaries',  'J1753-2240'],
['WD-NS binaries',  'J1807-2459B'],
['WD-NS binaries',  'J1853+1303'],
['WD-NS binaries',  'J1909-3744'],
['WD-NS binaries',  'J1910+1256'],
['WD-NS binaries',  'J1910-5959A'],
['WD-NS binaries',  'J1918-0642'],
['WD-NS binaries',  'J1946+3417'],
['WD-NS binaries',  'J1949+3106'],
['WD-NS binaries',  'J1950+2414'],
['WD-NS binaries',  'J2016+1948'],
['WD-NS binaries',  'J2019+2425'],
['WD-NS binaries',  'J1640-2224'],
['WD-NS binaries',  'J2017-0603'],
['WD-NS binaries',  'J2043+1711'],
['WD-NS binaries',  'J2302+4442'],
['WD-NS binaries',  'J2317+1439'],
['WD-NS binaries',  'J2045+3633'],
['WD-NS binaries',  'J2053+4650'],
['WD-NS binaries',  'J2145-0750'],
['WD-NS binaries',  'J2234+0611'],
['WD-NS binaries',  'XTE J2123-058/V* LZ Aqr '],
['WD-NS binaries GalCluster pulsar',  'B1516+02B/J1518+0204B (M5)'],
['WD-NS binaries GalCluster pulsar',  'J0024-7204H/B0021-72H (47 Tucanae)'],
['WD-NS binaries GalCluster pulsar',  'J0024-7204I/B0021-72I (47 Tucanae)'],
['WD-NS binaries GalCluster pulsar',  'J1748-2021B / NGC 6440B'],
['WD-NS binaries GalCluster pulsar',  'J1750-37A / NGC 6441A'],
['WD-NS binaries GalCluster pulsar',  'J1802-2124'],
['WD-NS binaries GalCluster pulsar',  'J1824-2452C / M28'],
['WD-NS binaries GalCluster pulsar',  'J1911-5958A / NGC 6752'],
['WD-NS binaries GalCluster pulsar',  'Ter 5ai / J1748-2446ai'],
['WD-NS binaries GalCluster pulsar',  'Ter 5I / J1748-2446I '],
['WD-NS binaries GalCluster pulsar',  'Ter 5J / J1748-2446J '],
['WD-NS binaries GalCluster pulsar',  'Ter 5U / J1748-2446U'],
['WD-NS binaries GalCluster pulsar',  'Ter 5W / J1748-2446W'],
['WD-NS binaries GalCluster pulsar',  'Ter 5X / J1748-2446X'],
['WD-NS binaries GalCluster pulsar',  'Ter 5Z / J1748-2446Z'],

['WD-NS binaries',  'J1537-5312'],
['WD-NS binaries',  'J1547-5709'],
['WD-NS binaries',  'J1618-4624'],
['WD-NS binaries',  'J1727-2951'],
['WD-NS binaries',  'J1755-25'],
['WD-NS binaries',  'J1244-6359'],
['WD-NS binaries',  'J1101-6424'],
['WD-NS binaries',  'J1614-2230'],

['WD-NS binaries',  'J1125+7819'],
['WD-NS binaries',  'J1641+8049'],
['WD-NS binaries',  'J1938+6604']
]

NS_radio_mag_param = [
[0, 0, 0, 0],
[26.3821, 39.77, 0.66, 9.70E+09],
[16.9405, 69.44, 15.62, 2.28E+10],
[32.7554, 69.22, 0, 1.25E+10],
[0.9378, 18.96, 0, 7.88E+11],
[21.8427, 0, 0, 2.95E+09],
[13.0648, 0, 0, 2.49E+10],
[200.3777, 40.99, 0, 5.97E+07],
[44.0541, 0, 1.94, 6.40E+09],
[0.3606, 0, 1.57, 1.59E+12],
[2.5387, 0, 240, 1.32E+12],
[24.4290, 7.43, 3.72, 1.07E+09],
[35.1351, 0, 0.32, 5.45E+09],
[46.5176, 0, 95.64, 7.61E+09],
[238.8810, 0, 0, 9.80E+08],
[9.5986, 0, 25.39, 9.80E+09],
[313.1749, 28.46, 0, 3.75E+08],
[24.3844, 0.25, 0, 1.48E+09],
[465.1352, 0, 63.7, 2.04E+08],
[6.9409, 49.28, 30.12, 1.73E+12],
[36.6502, 0, 1.02, 2.12E+09],
[5.3902, 0, 0, 5.85E+10],
[30.4712, 0, 0, 7.60E+08],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[1.079591939, 3564.09, 1069.23, 2.06E+12],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[160.8097, 10.06, 1.01, 8.53E+08],
[90.2873, 48.6, 5.18, 2.76E+09],
[3.0271, 74.63, 41.7, 1.19E+12],
[0.7478, 5.13, 1.69, 1.36E+12],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0.2657, 0, 0, 2.60E+14],
[5.3368, 0, 0, 0],
[0.6764, 0, 27.01, 4.38E+12],
[1.0579, 0, 0, 3.48E+12],
[2.4048, 0, 112, 6.93E+12],
[602.4096, 111.49, 0, 8.84E+07],
[1.8576, 0, 96.85, 7.36E+11],
[3.5817, 0, 57.17, 2.55E+12],
[0, 0, 0, 0],
[96.36223457, 0, 0, 0],
[268.6668662, 0, 0, 0],
[283.4409432, 0, 0, 0],
[320.6886957, 0, 0, 0],
[402.0937023, 0, 0, 0],
[332.9448049, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[0, 0, 0, 0],
[707.3144458, 0, 0, 8.31E+07],
[0, 0, 0, 0],
[43.2884, 188.6, 60.84, 3.32E+09],
[186.4941, 28.8, 5.76, 3.13E+08],
[622.1220305, 39.2, 0.78, 1.67E+08],
[365.9534, 0, 0, 2.22E+08],
[25.5606, 0, 0, 3.11E+09],
[173.6879, 13.52, 3.69, 5.81E+08],
[34.6574, 1.68, 0.34, 1.18E+09],
[346.5320, 42.98, 1.45, 1.90E+08],
[287.4579, 12.32, 3.94, 1.67E+08],
[500.2501, 0, 0, 1.70E+08],
[190.2678, 14.7, 1.86, 3.04E+08],
[60.7794, 31.22, 2.04, 8.54E+08],
[592.4215, 0, 0, 1.09E+08],
[133.7931, 1.73, 0.32, 3.68E+08],
[592.9878, 0, 0, 1.38E+08],
[277.9377, 0, 8.69, 1.87E+08],
[317.3789, 0, 0.59, 1.76E+08],
[218.8118, 11.68, 15.63, 2.00E+08],
[170.9374, 0, 1.45, 3.80E+08],
[0, 0, 0, 0],
[266.8692, 8.89, 2.58, 3.41E+08],
[104.4911, 0, 0, 0],
[10.51106825, 0, 1.57, 9.72E+09],
[238.8814, 0, 0, 5.94E+08],
[244.3913777, 22.22, 1.8, 1.91E+08],
[339.3157, 0, 3.25, 2.06E+08],
[200.6588053, 0, 1.16, 2.22E+08],
[306.1674, 0, 3.36, 9.93E+07],
[130.7895, 7.28, 0.72, 4.48E+08],
[315.4436, 0, 43.35, 1.01E+08],
[76.1140, 0, 12.83, 1.12E+09],
[232.3002, 0, 6.29, 2.88E+08],
[15.39873763, 15.43, 0, 5.16E+09],
[254.1603, 0, 0, 1.68E+08],
[316.1239842, 18.37, 1.06, 9.55E+07],
[345.2781365, 0, 0.35, 1.54E+08],
[420.1894, 4.44, 0.54, 1.13E+08],
[192.5919636, 8.42, 1.04, 2.72E+08],
[290.2546082, 20.28, 1.67, 9.26E+07],
[31.5638, 0, 0, 4.37E+09],
[79.4516, 0, 0, 1.49E+09],
[62.2959, 23.47, 5.26, 7.00E+08],
[279.5966, 1.23, 0.57, 2.10E+08],
[0, 0, 0, 0],
[125.8346, 28.12, 1.41, 0],
[311.4934, 0, 0, 0],
[0, 0, 0, 0],
[59.66541822, 0, 0, 0],
[8.960506248, 0, 0, 0],
[79.0664, 0, 0.44, 9.69E+08],
[240.442414, 0, 0, 0],
[0, 0, 0, 0],
[47.1068, 0, 0, 0],
[104.4911, 0, 0, 0],
[12.4474, 0, 0, 0],
[304.0308, 0, 0, 0],
[237.8019, 0, 0, 0],
[333.4156, 0, 0, 0],
[406.0765, 0, 0, 0],

[0.0331,  0, 0, 0],#J1537-5312 \citep{Cameron:2020pin}
[0.0179,  0, 0, 0],#J1547-5709 \citep{Cameron:2020pin}
[0.0136,  0, 0, 0],#J1618-4624 \citep{Cameron:2020pin}
[0.30,    0, 0, 0],#J1727-2951 \citep{Cameron:2020pin}

[315.19598238,    0, 0, 0],#J1755-25 \citep{Ng:2015zza}
[147.27431048,    0, 0, 0],#J1244-6359 \citep{Ng:2015zza}
[5.109272904279,  0, 0, 0],#J1101-6424 \citep{Ng:2015zza}
[3.1508076534271, 0, 0, 0],#J1614-2230 \citep{Ng:2015zza}

[13.94, 0, 0.90, 1.73E+08],#J1614-2230 \citep{Ng:2015zza}
[0.0, 0, 0, 1.36E+08],#J1614-2230 \citep{Ng:2015zza}
[0.0, 0, 0, 6.66E+08]#J1614-2230 \citep{Ng:2015zza}
]


#mollweid_vec(90,30,-60,-60,'black',0.5,6,6)
#mollweid_vec(10,45,120,-90,'orange',0.5,6,6)
#mollweid_vec(-70,-45,30,110,'red',0.5,6,6)

l = np.random.uniform(-180, 180, 20000)
b = np.random.uniform(-90, 90, 20000)

#hpx_map = cat2hpx(l, b, nside=32, radec=False)
hpx_map = cat2hpx(l, b, nside=64, radec=False)
hp.mollview(np.log10(hpx_map+1), cbar = True, coord='G')
hp.graticule(dpar = 20, dmer = 20, linestyle = ':', color = 'lightgrey')
plt.savefig('plots_v2/healpix_generated_random_test.png', transparent=True)
plt.savefig('plots_v2/healpix_generated_random_test.pdf', transparent=True)
plt.savefig('plots_v2/healpix_generated_random_test.eps', transparent=True)
plt.clf()




#cmap=plt.get_cmap('plasma')
#cmap=plt.get_cmap('inferno')
#cmap=plt.get_cmap('YlOrBr')
#cmap=plt.get_cmap('YlGn')
#cmap=plt.get_cmap('tab20b')
#cmap=plt.get_cmap('gist_ncar')
#cmap=plt.get_cmap('nipy_spectral')
cmap=plt.get_cmap('CMRmap')


#NSIDE = 64
NSIDE = 512
npix = np.arange(hp.nside2npix(NSIDE))*1

#You can find nside from the length of the image array by calling
#nside = hp.npix2nside(npix)
_path      = '/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/healpix_fit_maps/'
file15     = 'wmap_band_iqumap_r9_7yr_W_v4.fits'
#file2      = 'exposure_healpix_to_w637_evclass_128_evtype_48_PSF23.fits.gz'

#file2      = 'new_healpix_exposure_rebinned_using_exposure_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'new_healpix_intensity_Alm_using_counts_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'new_healpix_intensity_integrated_down_smoothed_using_counts_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'new_healpix_intensity_integrated_down_using_counts_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'new_healpix_intensity_integrated_up_smoothed_using_counts_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'new_healpix_intensity_integrated_up_using_counts_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'new_healpix_intensity_smoothed_using_counts_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'new_healpix_intensity_using_counts_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'counts_healpix_full_clean_gt20MeV.fits'
#file2      = 'counts_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'exposure_healpix_full_clean_gt20MeV.fits'
#file2      = 'exposure_healpix_to_w714_evclass_128_evtype_48_PSF23_zmax_100_order9.fits'
#file2      = 'wmap_band_iqumap_r9_7yr_W_v4.fits'
#file2      = 'wmap_temperature_analysis_mask_r9_7yr_v4.fits'
#file2      = '12yr_healpix_exposure_rebinned_using_exposure_healpix_to_w637_evclass_128_evtype_48_PSF23.fits.gz'
#file2      = '12yr_healpix_intensity_integrated_down_smoothed_using_counts_healpix_to_w637_evclass_128_evtype_48_PSF23_order9.fits.gz'
#file2      = '12yr_healpix_intensity_integrated_down_using_counts_healpix_to_w637_evclass_128_evtype_48_PSF23_order9.fits.gz'

#file2      = '12yr_healpix_intensity_integrated_up_smoothed_using_counts_healpix_to_w637_evclass_128_evtype_48_PSF23_order9.fits.gz'
#file2      = '12yr_healpix_intensity_integrated_up_using_counts_healpix_to_w637_evclass_128_evtype_48_PSF23_order9.fits.gz'
#file2      = '12yr_healpix_intensity_smoothed_using_counts_healpix_to_w637_evclass_128_evtype_48_PSF23_order9.fits.gz'
#file2      = '12yr_healpix_intensity_using_counts_healpix_to_w637_evclass_128_evtype_48_PSF23_order9.fits.gz'

#file2      = 'LFI_SkyMap_070_1024_R3.00_survey-1-3-5-6-7-8.fits'
#file2      = 'HFI_SkyMap_857_2048_R3.01_full.fits'
#file2      = 'DIRBE_4_256.fits'
#file2      = 'DIRBE_8_256.fits'
#file2      = 'DIRBE_ZSMA_3_256.fits'
#file2      = 'IRIS_NOHOLE_1_2048.fits'
#file2      = 'IRIS_NOHOLE_2_2048.fits'
#file2      = 'DIRBE_ZSMA_6_256.fits'
#file2      = 'IRIS_NOHOLE_3_2048.fits'
file2      = 'DIRBE_ZSMA_8_256.fits'
#file2      = 'irsa_lambda_sfd_ebv.fits'


#extMapFile = 'test_healpix_exposure_rebinned_using_exposure_healpix_full_clean_gt20MeV.fits'
extMapFile = '12yr_healpix_exposure_rebinned_using_exposure_healpix_to_w637_evclass_128_evtype_48_PSF23.fits.gz'
_path_clumpy = '/mnt/dampe_data/CLUMPY/output/'
#file_clumpy  = 'annihil_seed1332_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.fits'
#file_clumpy  = 'annihil_seed1332_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.fits'
#file_clumpy  = 'annihil_seed4444_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.fits'
file_clumpy  = 'annihil_seed3120_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.fits'
#_path_clumpy = '/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/DMASS-publish_v1.0/'
#file_clumpy = 'cmass-dr12v4-S-Reid-full.dat.fits'
#file_clumpy = 'Y1LSSmask_v2_redlimcut_il22_seeil4.0_4096ring.fits'

exposure_hpx = Map.read('/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/healpix_fit_maps/exposure_healpix_full_clean_gt20MeV.fits')
#print(events)
print(exposure_hpx.geom)
print(exposure_hpx.geom.axes[0])

#map15  = hp.read_map(os.path.join(_path,file15), h=True)

#map2   = hp.read_map(os.path.join(_path,file2))
map2   = hp.read_map(os.path.join(_path_clumpy,file_clumpy))


#extMap = hp.read_map(os.path.join(_path,extMapFile), h=True)

#print('extMap',extMap)
print('extMap',map2)
print('extMap size', len(map2))
#fig=plt.figure()

#hp.mollview(title="Galactic map of NS", cmap = cmap, unit = '$M_{T}$[$M_{\odot}$]', coord=['G','E'], cbar = True)
#hp.mollview(npix, title="Galactic map of NS", cmap = cmap, unit = '$M_{T}$[$M_{\odot}$]', cbar = True, coord=['G','E'])
#tick_locator = ticker.MaxNLocator(nbins=10)
hp.mollview(map2,
            #title="Galactic map of NS",
            #title="Galactic map of DM clumps + NSs. Einasto DM-profile.",
            title="",
            cmap = cmap,
            unit = 'spectral flux density [$MJy/sr$]',
            #unit = '$J_{tot}$ [GeV$^{2}$ cm$^{-5}$]',
            #cbar = True, coord=['G','E'], norm='log',
            cbar = True, coord='G', norm='log',
            #cbar = False, coord='G', norm='log',
            notext=False,
            #min=1e+15, max=1e18)
            min=1e+15, max=1e19)
            #min=1, max=6.3)
#fig = plt.gcf()
#if not ax:
#ax = plt.gca()
##fig.colorbar()
#cbar = plt.colorbar(ax, ticks=[-1, 0, 1], orientation='horizontal')  # set some values to ticks
#cbar.ax.set_xticklabels(['negative', 'zero', 'positive'])

#cbar =plt.colorbar(hp,
    #orientation='horizontal',
    #spacing='uniform',
    ##extend='max',
    ##label=r'NS mass $M_{T}[M_{\odot}]$ '#,
    #label='spectral flux density [$MJy/sr$]'
    ##boundaries=[0.2,2.8]
    #)
#cbar.locator = tick_locator
#cbar.update_ticks()

#hp.graticule(dpar = 30)
#hp.graticule(dpar = 90)
#hp.graticule(dpar = 20, dmer = 20, linestyle = ':', color = 'lightgrey')
hp.graticule(dpar=15, dmer=20, linestyle=':', color='black')#, local=False)
#plt.savefig('plots_v2/healpix_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.png', transparent=True)
#plt.savefig('plots_v2/healpix_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_gal2D_LOS0_90_FOVdiameter360.eps', transparent=True)
#plt.savefig('plots_v2/healpix_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_axial_gal2D_LOS0_90_FOVdiameter360.png', transparent=True)
#plt.savefig('plots_v2/healpix_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_axial_gal2D_LOS0_90_FOVdiameter360.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM500GeV_axial_gal2D_LOS0_90_FOVdiameter360.eps', transparent=True)
plt.savefig('plots_v2/healpix_seed3120_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_axial_gal2D_LOS0_90_FOVdiameter360.png', transparent=True)
plt.savefig('plots_v2/healpix_seed3120_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_axial_gal2D_LOS0_90_FOVdiameter360.pdf', transparent=True)
plt.savefig('plots_v2/healpix_seed3120_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_axial_gal2D_LOS0_90_FOVdiameter360.eps', transparent=True)

plt.clf()


print('number of NSs: ', len(column(l_GB_GL_dist,1)))



MAX = 900
MIN = 0


data=np.arange( len(column(l_GB_GL_dist,1)) )
theta_ns=np.arange( len(column(l_GB_GL_dist,1)) )
phi_ns=np.arange( len(column(l_GB_GL_dist,1)) )


lats = np.array(column(l_GB_GL_dist,0))
lons = np.array(column(l_GB_GL_dist,1))


theta_ns_=np.arange(len(column(l_GB_GL_dist,0)))
phi_ns_=np.arange(len(column(l_GB_GL_dist,0)))

for i in range(len(column(RA_DEC_deg_dist_kpc,0))):
    ##theta is the co-latitude, i.e. at the North Pole, at the Equator,at the South Pole
    ##phi is the longitude
    theta_ = 0.5*np.pi - np.deg2rad(column(RA_DEC_deg_dist_kpc,1)[i])
    phi_   = np.deg2rad(column(RA_DEC_deg_dist_kpc,0)[i])
    np.append(theta_ns_, theta_)
    np.append(phi_ns_, phi_)

#lons, lats = np.meshgrid(lons, lats)
vecs = coord.ang2vec(lons.flatten(), lats.flatten())   # or better directly: coord.rand_vec(ncrs)
#vecs = coord.ang2vec(phi_ns_.flatten(), theta_ns_.flatten())   # or better directly: coord.rand_vec(ncrs)
#vecs = coord.sph_unit_vectors(phi_ns_.flatten(), theta_ns_.flatten())   # or better directly: coord.rand_vec(ncrs)
#vecs = coord.ang2vec(phi_ns_, theta_ns_)   # or better directly: coord.rand_vec(ncrs)
#vecs = coord.sph_unit_vectors(phi_ns_, theta_ns_)   # or better directly: coord.rand_vec(ncrs)

#vecs = coord.ang2vec(column(l_GB_GL_dist,1), column(l_GB_GL_dist,0))   # or better directly: coord.rand_vec(ncrs)
#vecs = coord.ang2vec(column(l_GB_GL_dist,0), column(l_GB_GL_dist,1))   # or better directly: coord.rand_vec(ncrs)

print('vecs:',vecs)

#cmap=plt.get_cmap('plasma')
#cmap=plt.get_cmap('inferno')
#cmap=plt.get_cmap('YlOrBr')
#cmap=plt.get_cmap('YlGn')
#cmap=plt.get_cmap('tab20b')
#cmap=plt.get_cmap('gist_ncar')
#cmap=plt.get_cmap('nipy_spectral')
#cmap=plt.get_cmap('CMRmap')
cmap=plt.get_cmap('YlGnBu')
gridcolor=plt.get_cmap('Greys')

#vecs=coord.eq2gal(vecs)
#vecs=coord.gal2eq(vecs)
#vecs=coord.gal2sgal(vecs)
#vecs=coord.sgal2gal(vecs)
#vecs=coord.sgal2eq(vecs)
#vecs=coord.eq2sgal(vecs)
#vecs=coord.ecl2eq(vecs)
#vecs=coord.eq2ecl(vecs)


"""


# Plot an example map with sampled NS masses. If you specify the opath keyword in
# the skymap function, the plot will be automatically saved and closed
log10e = auger.rand_energy_from_auger(n=len(column(l_GB_GL_dist,1)), log10e_min=min(column(NS_and_compStar_mass_and_unc,0)))
#skymap.scatter(vecs, c=log10e, opath='plots_v2/healpix_NS_map_with_astrotools.pdf')
#fig, ax = skymap.scatter(vecs, s=100, cmap=cmap, tickalpha=0.2, plane='GP', coord_system='gal', c=log10e, cblabel=r'$M_{T}[M_{\odot}]$ ')
fig, ax = skymap.scatter(vecs, s=100, cmap=cmap, fontsize=16,
#fig, ax = skymap.scatter(vecs, cmap=cmap, fontsize=16,
                         dark_grid=True,
                         tickalpha=0.6,
                         gridcolor='lightgrey',
                         gridalpha=0.2,
                         #plane={'SGP','GP'},
                         #plane='SGP',
                         #plane='GP',
                         #coord_system='gal',
                         coord_system='eq',
                         c=log10e, cblabel=r'$M_{T}[M_{\odot}]$ ')

#plt.scatter()

#nside = 64      # resolution of the HEALPix map (default: 64)
#npix = hp.nside2npix(nside)
#index = hp.pix2ang(np.array(l_glon), np.array(l_glat), lonlat=True)
#count_map = np.histogram(index, bins=np.arange(npix + 1))[0]

#lats = np.array(l_glat)
#lons = np.array(l_glon)
##lons, lats = np.meshgrid(lons, lats)
#PHI, THETA = np.meshgrid(lons, lats)
#print('THETA:',  THETA)
#print('PHI:',    PHI)
#grid_pix = hp.ang2pix(nside, THETA, PHI)
#grid_map = count_map[grid_pix]


#image = ax.pcolormesh(l_glon, l_glat, grid_map,
                      #rasterized=True,
                      #cmap=plt.get_cmap('YlOrBr'), shading='auto'
                      #)
#cb = fig.colorbar(image, orientation='horizontal', shrink=.6, pad=0.05)
#cb.set_label(r'$M_{T}[M_{\odot}]$', size=22)

#ax = plt.axes(projection='geo aitoff')
ax.set_xticklabels([r"150$\degree$", r"120$\degree$", r"90$\degree$", r"60$\degree$", r"30$\degree$", r"GC", r"330$\degree$", r"300$\degree$", r"270$\degree$", r"240$\degree$", r"210$\degree$"])
ax.set_xticklabels([])
ax.set_title("Galactic")
width = 17# width of the figure
labelsize = 16
ax.set_xlabel('g.lon', size=labelsize)
ax.set_ylabel('g.lat', size=labelsize)
ax.grid(True, alpha=0.25)

plt.savefig('plots_v2/healpix_NS_map_with_astrotools.png')#, transparent=True)
plt.savefig('plots_v2/healpix_NS_map_with_astrotools.pdf')#, transparent=True)
plt.savefig('plots_v2/healpix_NS_map_with_astrotools.eps')#, transparent=True)
plt.clf()
"""



xNS_list = np.array(column(l_GB_GL_dist,1))
yNS_list = np.array(column(l_GB_GL_dist,0))
zNS_list = np.array(column(l_GB_GL_dist,2))
coordsNS = SkyCoord(l=xNS_list, b=yNS_list, unit='deg', frame="galactic")
scale_NS=[]
for i in range(len(column(RA_DEC_deg_dist_kpc,0))):
    scale_NS.append(10* column(NS_and_compStar_mass_and_unc,0)[i])


x_list = np.array(column(l_GB_GL_dist,1))
y_list = np.array(column(l_GB_GL_dist,0))
z_list = np.array(column(l_GB_GL_dist,2))
coords = SkyCoord(l=x_list, b=y_list, unit='deg', frame="galactic")
#scale_NSs_mass= [10* column(NS_and_compStar_mass_and_unc,0)]
scale_NSs_mass=[]
for i in range(len(column(RA_DEC_deg_dist_kpc,0))):
    scale_NSs_mass.append(20* column(NS_and_compStar_mass_and_unc,0)[i])




#file_name= '/mnt/dampe_data/CLUMPY/root_ntuples/ntuple_annihil_seed7061_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root'
#file_name= '/mnt/dampe_data/CLUMPY/root_ntuples/ntuple_annihil_seed9111_maxNSubclumps300_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root'
#file_name= '/mnt/dampe_data/CLUMPY/root_ntuples/ntuple_annihil_seed4444_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root'
#file_name= '/mnt/dampe_data/CLUMPY/root_ntuples/ntuple_annihil_seed3120_minSubClumpM_m07_maxSubClumpM_01_minNClumps10_mDM100GeV_nClumps300_massive_axial_gal2D_LOS0_90_FOVdiameter360.0deg_nside2048.root'
file_name= '/mnt/dampe_data/CLUMPY/root_ntuples/ntuple_annihil_seed1332_maxNSubclumps200_minSubClumpM_m07_maxSubClumpM01_mDM100GeV_gal2D_LOS180_0_FOVdiameter360.0deg_nside2048.root'
#file_name= '/mnt/dampe_data/CLUMPY/analysis/mergedALL_NS-DMclump_mDM100GeV_nside2048_r005_v5.root'
#file_name= '/mnt/dampe_data/CLUMPY/analysis/mergedALL_NS-DMclump_mDM500GeV_nside2048_r005_v5.root'
#tree_1   = uproot.open(file_name)["clumps_1kpc"]
#tree_075 = uproot.open(file_name)["clumps_075kpc"]
#tree_05  = uproot.open(file_name)["clumps_05kpc"]
#tree_025 = uproot.open(file_name)["clumps_025kpc"]
#tree_01  = uproot.open(file_name)["clumps_01kpc"]
#tree_005 = uproot.open(file_name)["clumps_005kpc"]
#tree_001 = uproot.open(file_name)["clumps_001kpc"]
tree = uproot.open(file_name)["CLUMPY_output"]

branches   = tree.arrays()

#branches = tree_005.arrays(namedecode='utf-8')
#branches_1   = tree_1.arrays()
#branches_075 = tree_075.arrays()
#branches_05  = tree_05.arrays()
#branches_025 = tree_025.arrays()
#branches_01  = tree_01.arrays()
#branches_005 = tree_005.arrays()
#branches_001 = tree_001.arrays()

branches = tree.arrays()
#branches.keys()
print("branches\n:  ",branches)
print("branches[GalacLongitude_deg]\n:  ",branches["GalacLongitude_deg"])
print("numentries\n:  ",len (branches['Distance_kpc']))


array_GLong     = np.array(branches['GalacLongitude_deg'])
array_GLat      = np.array(branches['GalacLatitude_deg'])
array_distance  = np.array(branches['Distance_kpc'])
array_distgal   = np.array(branches['DM_Dgal_kpc'])
array_Mdelta    = np.array(branches['DM_Mdelta_Msol'])
array_Mtidal    = np.array(branches['DM_Mtidal_Msol'])
array_Mequdens  = np.array(branches['DM_Mequdens_Msol'])
array_Rdelta    = np.array(branches['Rdelta_kpc'])
array_Rtidal    = np.array(branches['DM_Rtidal_kpc'])
array_Requdens  = np.array(branches['DM_Requdens_kpc'])


clumps_mass_1e2Msun = array_Mdelta >= 1e2
clumps_mass_5e2Msun = array_Mdelta >= 5e2
clumps_mass_1e3Msun = array_Mdelta >= 1e3
clumps_mass_1e4Msun = array_Mdelta >= 1e4
clumps_mass_1e5Msun = array_Mdelta >= 1e5
clumps_mass_1e6Msun = array_Mdelta >= 1e6
clumps_mass_1e7Msun = array_Mdelta >= 1e7
clumps_mass_1e8Msun = array_Mdelta >= 1e8
clumps_mass_1e9Msun = array_Mdelta >= 1e9

print("#DMclumps array_Mdelta >= 1e2: ",clumps_mass_1e2Msun.sum())
print("#DMclumps array_Mdelta >= 5e2: ",clumps_mass_5e2Msun.sum())
print("#DMclumps array_Mdelta >= 1e3: ",clumps_mass_1e3Msun.sum())
print("#DMclumps array_Mdelta >= 1e4: ",clumps_mass_1e4Msun.sum())
print("#DMclumps array_Mdelta >= 1e5: ",clumps_mass_1e5Msun.sum())
print("#DMclumps array_Mdelta >= 1e6: ",clumps_mass_1e6Msun.sum())
print("#DMclumps array_Mdelta >= 1e7: ",clumps_mass_1e7Msun.sum())
print("#DMclumps array_Mdelta >= 1e8: ",clumps_mass_1e8Msun.sum())
print("#DMclumps array_Mdelta >= 1e9: ",clumps_mass_1e9Msun.sum())


lats_list = branches['GalacLatitude_deg'][clumps_mass_1e7Msun]
lons_list = branches['GalacLongitude_deg'][clumps_mass_1e7Msun]
mass_list = branches['DM_Mdelta_Msol'][clumps_mass_1e7Msun]

data_pix = np.arange(clumps_mass_1e7Msun.sum())
#data_pix = branches['DM_Mdelta_Msol'][clumps_mass_1e7Msun]
theta    = np.pi/180.*(90.-np.array(lats_list))
phi      = np.pi/180.*np.array(lons_list)

pixel_indices  = hp.ang2pix(2048, theta, phi)
map_clumps_pix = np.zeros(hp.nside2npix(2048))
map_clumps_pix[pixel_indices] = data_pix
#for i in range(clumps_mass_1e7Msun.sum()):
#    map_clumps_pix[pixel_indices[i]] += data_pix[i]


clumps_pixel_mask=np.arange(clumps_mass_1e7Msun.sum())
clumps_pixel_mask=map_clumps_pix.astype(bool)
ns_mask = cm.make_mask_total(nside=2048, custom_mask=clumps_pixel_mask) #why it works only with nside=128?
apodized_mask5rad  = np.clip(hp.smoothing(ns_mask, fwhm=np.radians(5)), 0, None)


cmap=plt.get_cmap('inferno')
#cmap=plt.get_cmap('plasma')
#cmap=plt.get_cmap('YlOrBr')
#cmap=plt.get_cmap('YlGn')
#cmap=plt.get_cmap('tab20b')
#cmap=plt.get_cmap('gist_ncar')
#cmap=plt.get_cmap('nipy_spectral')
#cmap=plt.get_cmap('CMRmap')

hp.mollview(apodized_mask5rad,
            #title="DM clumps with $M > 10^{7} [M_{\odot}]$ with Gaussian smearing level: $2^{\circ}$",
            title="DM clumps with $M > 10^{7} [M_{\odot}]$ with Gaussian smearing level: $5^{\circ}$",
            #title="DM clumps with $M > 10^{7} [M_{\odot}]$ with Gaussian smearing level: $1^{\circ}$",
            #title="DM clumps with $M > 10^{3} [M_{\odot}]$ with Gaussian smearing level: $1^{\circ}$",
            #title="DM clumps with $M > 10^{3} [M_{\odot}]$ with Gaussian smearing level: $5^{\circ}$",
            #title="DM clumps with $M > 10^{4} [M_{\odot}]$ with Gaussian smearing level: $5^{\circ}$",
            #title="DM clumps with $M > 10^{5} [M_{\odot}]$ with Gaussian smearing level: $5^{\circ}$",
            #title="DM clumps with $M > 10^{7} [M_{\odot}]$ with Gaussian smearing level: $5^{\circ}$",
            cmap=cmap,
            #unit = 'DM clumps with $M[M_{\odot}]$',
            #cbar = True, coord=['G','E'], #norm='log',
            cbar = False,
            #cbar = True,
            coord='G', #norm='log',
            notext=False)#,#
            #min=1e0, max=1e4)
"""
hp.projview(
    apodized_mask5rad,
    cmap=cmap,
    coord=["G"],
    graticule=True,
    graticule_labels=True,
    unit="DM clumps mass $M[M_{\odot}]$",
    xlabel="g.lon",
    ylabel="g.lat",
    latitude_grid_spacing=15,
    longitude_grid_spacing=30,
    cb_orientation="horizontal",
    min=1e+2,
    max=1e+8,
    projection_type="mollweide"
)

hp.projview(
    apodized_mask5rad,
    cmap=cmap,
    coord=["G"],
    graticule=True,
    graticule_labels=True,
    unit=r"DM clumps mass $M[M_{\odot}]$",
    xlabel="g.lon",
    ylabel="g.lat",
    cb_orientation="horizontal",
    min=1e0,
    max=1e3,
    latitude_grid_spacing=15,
    longitude_grid_spacing=30,
    projection_type="mollweide"
    #title="Hammer projection",
    #fontsize={
        #"xlabel": 20,
        #"ylabel": 20,
        #"xtick_label": 20,
        #"ytick_label": 20,
        #"title": 20,
        #"cbar_label": 20,
        #"cbar_tick_label": 20,
    #},
    #xtick_label_color="r",
    #ytick_label_color="g",
    #graticule_color="black",
)
"""

hp.graticule(dpar=15, dmer=30, linestyle=':', color='black')#, local=False)
fontsize_ = 8
f = plt.gcf()# accessing the current figure...
HpxAx = f.get_children()[1]
coord_text_obj = HpxAx.get_children()[0]
#coord_text_obj.set_fontsize(fontsize_)

"""
tick_locator = ticker.MaxNLocator(nbins=8)
cbar =plt.colorbar(sc,
    orientation='horizontal',
    spacing='uniform',
    #extend='max',
    #label=r'NS mass $M_{T}[M_{\odot}]$ '#,
    label=r'$M[M_{\odot}]$ '#,
    #label=r'distance $[kpc]$ '#,
    #boundaries=[0.2,2.8]
    )
cbar.locator = tick_locator
cbar.update_ticks()
fig.tight_layout()
"""

#plt.savefig('plots_v2/healpix_seed1332_DMclumps1e3_mask_smeared1deg.png', transparent=True)
#plt.savefig('plots_v2/healpix_seed1332_DMclumps1e3_mask_smeared1deg.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_seed1332_DMclumps1e3_mask_smeared1deg.eps', transparent=True)
plt.savefig('plots_v2/healpix_seed1332_DMclumps1e7_mask_smeared5deg.png', transparent=True)
plt.savefig('plots_v2/healpix_seed1332_DMclumps1e7_mask_smeared5deg.pdf', transparent=True)
plt.savefig('plots_v2/healpix_seed1332_DMclumps1e7_mask_smeared5deg.eps', transparent=True)

#plt.savefig('plots_v2/healpix_seed4444_DMclumps_mask_smeared5deg.png', transparent=True)
#plt.savefig('plots_v2/healpix_seed4444_DMclumps_mask_smeared5deg.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_seed4444_DMclumps_mask_smeared5deg.eps', transparent=True)
#plt.savefig('plots_v2/healpix_seed4444_DMclumps1e3_mask_smeared5deg.png', transparent=True)
#plt.savefig('plots_v2/healpix_seed4444_DMclumps1e3_mask_smeared5deg.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_seed4444_DMclumps1e3_mask_smeared5deg.eps', transparent=True)
#plt.savefig('plots_v2/healpix_seed4444_DMclumps1e7_mask_smeared5deg.png', transparent=True)
#plt.savefig('plots_v2/healpix_seed4444_DMclumps1e7_mask_smeared5deg.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_seed4444_DMclumps1e7_mask_smeared5deg.eps', transparent=True)

#plt.savefig('plots_v2/healpix_seed9111_DMclumps1e3_mask_smeared5deg.png', transparent=True)
#plt.savefig('plots_v2/healpix_seed9111_DMclumps1e3_mask_smeared5deg.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_seed9111_DMclumps1e3_mask_smeared5deg.eps', transparent=True)
#plt.savefig('plots_v2/healpix_seed3120_DMclumps1e4_mask_smeared5deg.png', transparent=True)
#plt.savefig('plots_v2/healpix_seed3120_DMclumps1e4_mask_smeared5deg.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_seed3120_DMclumps1e4_mask_smeared5deg.eps', transparent=True)
#plt.savefig('plots_v2/healpix_seed3120_DMclumps1e5_mask_smeared5deg.png', transparent=True)
#plt.savefig('plots_v2/healpix_seed3120_DMclumps1e5_mask_smeared5deg.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_seed3120_DMclumps1e5_mask_smeared5deg.eps', transparent=True)
#plt.savefig('plots_v2/healpix_seed3120_DMclumps1e7_mask_smeared5deg.png', transparent=True)
#plt.savefig('plots_v2/healpix_seed3120_DMclumps1e7_mask_smeared5deg.pdf', transparent=True)
#plt.savefig('plots_v2/healpix_seed3120_DMclumps1e7_mask_smeared5deg.eps', transparent=True)
plt.clf()








file_name= '/mnt/dampe_data/CLUMPY/analysis/mergedALL_NS-DMclump_mDM100GeV_nside2048_r005_v5.root'
#file_name= '/mnt/dampe_data/CLUMPY/analysis/mergedALL_NS-DMclump_mDM500GeV_nside2048_r005_v5.root'
tree_1   = uproot.open(file_name)["clumps_1kpc"]
tree_075 = uproot.open(file_name)["clumps_075kpc"]
tree_05  = uproot.open(file_name)["clumps_05kpc"]
tree_025 = uproot.open(file_name)["clumps_025kpc"]
tree_01  = uproot.open(file_name)["clumps_01kpc"]
tree_005 = uproot.open(file_name)["clumps_005kpc"]
tree_001 = uproot.open(file_name)["clumps_001kpc"]
#tree = uproot.open(file_name)["CLUMPY_output"]
#branches   = tree.arrays()

#branches = tree_005.arrays(namedecode='utf-8')
branches_1   = tree_1.arrays()
branches_075 = tree_075.arrays()
branches_05  = tree_05.arrays()
branches_025 = tree_025.arrays()
branches_01  = tree_01.arrays()
branches_005 = tree_005.arrays()
branches_001 = tree_001.arrays()

#branches = tree.arrays()
#branches.keys()
print("branches\n:  ",branches)
print("branches[GalacLongitude_deg]\n:  ",branches["GalacLongitude_deg"])
print("numentries\n:  ",len (branches['Distance_kpc']))





#branches.keys()
print("branches_005\n:  ",branches_005)
#print("branches_005[GalacLongitude_deg]\n:  ",branches_005["GalacLongitude_deg"])
#print("numentries\n:  ",len (branches_005['Distance_kpc']))

print("Ttree_005 show:\n",tree_005.show())
#print ("Ttree_005 arrays:\n", tree_005.arrays() )
print("")
print("Ttree_005 name:  ",tree_005.name)
print("Ttree_005 title: ",tree_005.title)
#print("Ttree_005 numentries:\n",tree_005.numentries)

print('tree_005.get(Distance_kpc):\n',tree_005.get('clumps_005kpc') )
print('tree_005.iterkeys(   ):\n',tree_005.iterkeys(recursive=False) )
print('tree_005.itervalues(   ):\n',tree_005.itervalues(recursive=False) )

################################################
#Test read some ROOT Ttree_005 branches_005
################################################
print ("Distance_kpc:               ", branches_005['distance_kpc'] )
print ("nEvents in Distance_kpc:    ", len (branches_005['distance_kpc']) )
print ("Distance_kpc.interpretation:", tree_005['distance_kpc'].interpretation )
print("")
print ("Distance_kpc:               ", branches_001['distance_kpc'] )
print ("nEvents in Distance_kpc:    ", len (branches_001['distance_kpc']) )
print ("Distance_kpc.interpretation:", tree_001['distance_kpc'].interpretation )
print("")
print ("Distance_kpc:               ", branches_01['distance_kpc'] )
print ("nEvents in Distance_kpc:    ", len (branches_01['distance_kpc']) )
print ("Distance_kpc.interpretation:", tree_01['distance_kpc'].interpretation )
print("")
print ("Distance_kpc:               ", branches_1['distance_kpc'] )
print ("nEvents in Distance_kpc:    ", len (branches_1['distance_kpc']) )
print ("Distance_kpc.interpretation:", tree_1['distance_kpc'].interpretation )
print("")
print ("clump GalacLongitude_deg[0]: ", branches_005['GalacLongitude_deg'][0])
print ("clump GalacLatitude_deg[0]:  ", branches_005['GalacLatitude_deg'][0])
print ("clump distance_kpc[0]:       ", branches_005['distance_kpc'][0])
print ("clump Mdelta_Msol[0]:        ", branches_005['DM_Mdelta_Msol'][0])

DMclumps001_GalacLong  = np.array(branches_001['GalacLongitude_deg'])
DMclumps001_GalacLatid = np.array(branches_001['GalacLatitude_deg'])
DMclumps001_dist       = np.array(branches_001['distance_kpc'])
DMclumps001_distSmooth = np.array(branches_001['distanceSmooth_kpc'])
DMclumps001_Mdelta     = np.array(branches_001['DM_Mdelta_Msol'])
DMclumps001_Rdelta     = np.array(branches_001['DM_Rdelta_kpc'])
DMclumps001_Mtilda     = np.array(branches_001['DM_Mtidal_Msol'])
DMclumps001_Rtilda     = np.array(branches_001['DM_Rtidal_kpc'])

DMclumps005_GalacLong  = np.array(branches_005['GalacLongitude_deg'])
DMclumps005_GalacLatid = np.array(branches_005['GalacLatitude_deg'])
DMclumps005_dist       = np.array(branches_005['distance_kpc'])
DMclumps005_distSmooth = np.array(branches_005['distanceSmooth_kpc'])
DMclumps005_Mdelta     = np.array(branches_005['DM_Mdelta_Msol'])
DMclumps005_Rdelta     = np.array(branches_005['DM_Rdelta_kpc'])
DMclumps005_Mtilda     = np.array(branches_005['DM_Mtidal_Msol'])
DMclumps005_Rtilda     = np.array(branches_005['DM_Rtidal_kpc'])

DMclumps01_GalacLong  = np.array(branches_01['GalacLongitude_deg'])
DMclumps01_GalacLatid = np.array(branches_01['GalacLatitude_deg'])
DMclumps01_dist       = np.array(branches_01['distance_kpc'])
DMclumps01_distSmooth = np.array(branches_01['distanceSmooth_kpc'])
DMclumps01_Mdelta     = np.array(branches_01['DM_Mdelta_Msol'])
DMclumps01_Rdelta     = np.array(branches_01['DM_Rdelta_kpc'])
DMclumps01_Mtilda     = np.array(branches_01['DM_Mtidal_Msol'])
DMclumps01_Rtilda     = np.array(branches_01['DM_Rtidal_kpc'])

DMclumps025_GalacLong  = np.array(branches_025['GalacLongitude_deg'])
DMclumps025_GalacLatid = np.array(branches_025['GalacLatitude_deg'])
DMclumps025_dist       = np.array(branches_025['distance_kpc'])
DMclumps025_distSmooth = np.array(branches_025['distanceSmooth_kpc'])
DMclumps025_Mdelta     = np.array(branches_025['DM_Mdelta_Msol'])
DMclumps025_Rdelta     = np.array(branches_025['DM_Rdelta_kpc'])
DMclumps025_Mtilda     = np.array(branches_025['DM_Mtidal_Msol'])
DMclumps025_Rtilda     = np.array(branches_025['DM_Rtidal_kpc'])

DMclumps05_GalacLong  = np.array(branches_05['GalacLongitude_deg'])
DMclumps05_GalacLatid = np.array(branches_05['GalacLatitude_deg'])
DMclumps05_dist       = np.array(branches_05['distance_kpc'])
DMclumps05_distSmooth = np.array(branches_05['distanceSmooth_kpc'])
DMclumps05_Mdelta     = np.array(branches_05['DM_Mdelta_Msol'])
DMclumps05_Rdelta     = np.array(branches_05['DM_Rdelta_kpc'])
DMclumps05_Mtilda     = np.array(branches_05['DM_Mtidal_Msol'])
DMclumps05_Rtilda     = np.array(branches_05['DM_Rtidal_kpc'])

DMclumps075_GalacLong  = np.array(branches_075['GalacLongitude_deg'])
DMclumps075_GalacLatid = np.array(branches_075['GalacLatitude_deg'])
DMclumps075_dist       = np.array(branches_075['distance_kpc'])
DMclumps075_distSmooth = np.array(branches_075['distanceSmooth_kpc'])
DMclumps075_Mdelta     = np.array(branches_075['DM_Mdelta_Msol'])
DMclumps075_Rdelta     = np.array(branches_075['DM_Rdelta_kpc'])
DMclumps075_Mtilda     = np.array(branches_075['DM_Mtidal_Msol'])
DMclumps075_Rtilda     = np.array(branches_075['DM_Rtidal_kpc'])

DMclumps1_GalacLong  = np.array(branches_1['GalacLongitude_deg'])
DMclumps1_GalacLatid = np.array(branches_1['GalacLatitude_deg'])
DMclumps1_dist       = np.array(branches_1['distance_kpc'])
DMclumps1_distSmooth = np.array(branches_1['distanceSmooth_kpc'])
DMclumps1_Mdelta     = np.array(branches_1['DM_Mdelta_Msol'])
DMclumps1_Rdelta     = np.array(branches_1['DM_Rdelta_kpc'])
DMclumps1_Mtilda     = np.array(branches_1['DM_Mtidal_Msol'])
DMclumps1_Rtilda     = np.array(branches_1['DM_Rtidal_kpc'])

x_list = np.concatenate((DMclumps001_GalacLatid, DMclumps005_GalacLatid, DMclumps01_GalacLatid, DMclumps025_GalacLatid, DMclumps05_GalacLatid, DMclumps075_GalacLatid, DMclumps1_GalacLatid))
y_list = np.concatenate((DMclumps001_GalacLong,  DMclumps005_GalacLong,  DMclumps01_GalacLong,  DMclumps025_GalacLong,  DMclumps05_GalacLong,  DMclumps075_GalacLong,  DMclumps1_GalacLong))
z_list = np.concatenate((DMclumps001_distSmooth, DMclumps005_distSmooth, DMclumps01_distSmooth, DMclumps025_distSmooth, DMclumps05_distSmooth, DMclumps075_distSmooth, DMclumps1_distSmooth))
#x_list = DMclumps005_GalacLatid
#y_list = DMclumps005_GalacLong
#z_list = DMclumps005_distSmooth
#scale_NSs_mass = DMclumps005_Mdelta
scale_NSs_mass = np.concatenate((DMclumps001_Mdelta, DMclumps005_Mdelta, DMclumps01_Mdelta, DMclumps025_Mdelta, DMclumps05_Mdelta, DMclumps075_Mdelta, DMclumps1_Mdelta))
#scale_NSs_mass = np.concatenate((DMclumps001_Rdelta, DMclumps005_Rdelta, DMclumps01_Rdelta, DMclumps025_Rdelta, DMclumps05_Rdelta, DMclumps075_Rdelta, DMclumps1_Rdelta))


#lats_list = branches['GalacLatitude_deg'][clumps_mass_1e7Msun]
#lons_list = branches['GalacLongitude_deg'][clumps_mass_1e7Msun]
#mass_list = branches['DM_Mdelta_Msol'][clumps_mass_1e7Msun]
lats_list = x_list
lons_list = y_list
mass_list = scale_NSs_mass

#data_pix = np.arange(clumps_mass_1e7Msun.sum())
#data_pix = branches['DM_Mdelta_Msol'][clumps_mass_1e7Msun]
#data_pix = np.arange(scale_NSs_mass)
data_pix = scale_NSs_mass
theta    = np.pi/180.*(90.-np.array(lats_list))
phi      = np.pi/180.*np.array(lons_list)

pixel_indices  = hp.ang2pix(2048, theta, phi)
map_clumps_pix = np.zeros(hp.nside2npix(2048))
map_clumps_pix[pixel_indices] = data_pix
#for i in range(clumps_mass_1e7Msun.sum()):
#    map_clumps_pix[pixel_indices[i]] += data_pix[i]


#clumps_pixel_mask=np.arange(clumps_mass_1e7Msun.sum())
clumps_pixel_mask=scale_NSs_mass
clumps_pixel_mask=map_clumps_pix.astype(bool)
ns_mask = cm.make_mask_total(nside=2048, custom_mask=clumps_pixel_mask) #why it works only with nside=128?
apodized_mask5rad  = np.clip(hp.smoothing(ns_mask, fwhm=np.radians(2)), 0, None)

cmap=plt.get_cmap('inferno')
#cmap=plt.get_cmap('plasma')
#cmap=plt.get_cmap('YlOrBr')
#cmap=plt.get_cmap('YlGn')
#cmap=plt.get_cmap('tab20b')
#cmap=plt.get_cmap('gist_ncar')
#cmap=plt.get_cmap('nipy_spectral')
#cmap=plt.get_cmap('CMRmap')

hp.mollview(apodized_mask5rad,
            #title="DM clumps with $M > 10^{7} [M_{\odot}]$ with Gaussian smearing level: $1^{\circ}$",
            #title="DM clumps that wraps the NSs",
            title="DM clumps with $M > 10^{3} [M_{\odot}]$ with Gaussian smearing level: $2^{\circ}$",
            #title="DM clumps with $M > 10^{3} [M_{\odot}]$ with Gaussian smearing level: $5^{\circ}$",
            #title="DM clumps with $M > 10^{4} [M_{\odot}]$ with Gaussian smearing level: $5^{\circ}$",
            #title="DM clumps with $M > 10^{5} [M_{\odot}]$ with Gaussian smearing level: $5^{\circ}$",
            #title="DM clumps with $M > 10^{7} [M_{\odot}]$ with Gaussian smearing level: $5^{\circ}$",
            cmap=cmap,
            #unit = 'DM clumps with $M[M_{\odot}]$',
            #cbar = True, coord=['G','E'], #norm='log',
            cbar = False,
            #cbar = True,
            coord='G', #norm='log',
            notext=False)#,#
            #min=1e0, max=1e4)

plt.savefig('plots_v2/healpix_100GeV_DMclumps_map_Mdelta_GalFrame_Gauss2deg_v1.png', transparent=True)
plt.savefig('plots_v2/healpix_100GeV_DMclumps_map_Mdelta_GalFrame_Gauss2deg_v1.pdf', transparent=True)
plt.savefig('plots_v2/healpix_100GeV_DMclumps_map_Mdelta_GalFrame_Gauss2deg_v1.eps', transparent=True)
plt.clf()



from astropy.io import fits

c1 = fits.Column(name='GalacLatid',   array=x_list,         format='K')
c2 = fits.Column(name='GalacLong,',   array=y_list,         format='K')
c3 = fits.Column(name='DMclump mass', array=scale_NSs_mass, format='K')
t = fits.BinTableHDU.from_columns([c1, c2, c3])
t.writeto('plots_v2/CLUMPY_test_table1.fits', overwrite=True)


"""
_path_clumpy = '/home/blessed/Documents/_Papers_ToSubmit/_DMStarsClose2GC_Minihalo/python_code/plots_v2/'
file_clumpy  = 'CLUMPY_test_table1.fits'
map2   = hp.read_map(os.path.join(_path_clumpy,file_clumpy))

#print(map2.geom)
#print(map2.geom.axes[0])
print('extMap',map2)
print('extMap size', len(map2))
hp.mollview(map2,
            #title="Galactic map of NS",
            #title="Galactic map of DM clumps + NSs. Einasto DM-profile.",
            title="",
            cmap = cmap,
            #unit = 'spectral flux density [$MJy/sr$]',
            #unit = '$J_{tot}$ [GeV$^{2}$ cm$^{-5}$]',
            #cbar = True, coord=['G','E'], norm='log',
            cbar = True, coord='G', norm='log',
            #cbar = False, coord='G', norm='log',
            notext=False)#,
            #min=1e+15, max=1e18)
            #min=1e+15, max=1e19)

hp.graticule(dpar=15, dmer=20, linestyle=':', color='black')#, local=False)
plt.savefig('plots_v2/healpix_100GeV_DMclumps_map_Mdelta_GalFrame_Gauss2deg_v1_.png', transparent=True)
plt.savefig('plots_v2/healpix_100GeV_DMclumps_map_Mdelta_GalFrame_Gauss2deg_v1_.pdf', transparent=True)
plt.savefig('plots_v2/healpix_100GeV_DMclumps_map_Mdelta_GalFrame_Gauss2deg_v1_.eps', transparent=True)
plt.clf()
"""
