#from root_numpy import root2array, tree2array
#from root_numpy import testdata
import uproot 
import ROOT

import numpy as np
from numpy import asarray
from numpy import savetxt
from numpy import save

from rfpimp import *
import pandas as pd
import seaborn as sns
import sys
import os
import h5py
import pickle
import math

# Remove Units from Value List
# Using replace() + strip() + list comprehension
import re

from healpy.newvisufunc import projview, newprojplot
# classic healpy mollweide projections plot with graticule and axis labels
from NPTFit import create_mask as cm # Module for creating masks

import astropy.units as u
import astropy.coordinates as apycoords
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
#from astropy.utils.compat import argparse

#If you didn't install tk in your python then you can not use the tk backend.
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from plotnine import *
from plotnine.data import mpg

M13_B_c = SkyCoord("16h41m40.39144s +36d25m58.4880s", frame=Galactic)  # Units from string
M13_D_c = SkyCoord("16h41m42.395232 +36d27m28.2021", frame=Galactic)  # Units from string
M13_E_c = SkyCoord("16h41m42.0221 +36d27m34.9676", frame=Galactic)  # Units from string
M13_F_c = SkyCoord("16h41m44.6058 +36d28m16.0034", frame=Galactic)  # Units from string

M13_A_c = SkyCoord("16h41m40.87019 +36d27m14.9788", frame=Galactic)  # Units from string
M13_C_c = SkyCoord("16h41m41.00748 +36d27m02.7438", frame=Galactic)  # Units from string

J2043_1711_c = SkyCoord("20h43m20.881730 17d11m28.91265", frame=Galactic)  # Units from string
J1807_2459B_c = SkyCoord("18h07m20.355604 -24d59m52.9015", frame=Galactic)  # Units from string
J1740_5340_c = SkyCoord("17h40m44.611 -53d40m41.57", frame=Galactic)  # Units from string
J1753_2240_c = SkyCoord("17h53m39.847 -22d40m42.", frame=Galactic)  # Units from string

RX_J1856_3754_c = SkyCoord("18h56m35.85 -37d54m36.7", frame=Galactic)  # Units from string



J1807_2500B_c = SkyCoord("18h07m20.871209 -25d00m1.915", frame=Galactic)  # Units from string

KS1731_260_c= SkyCoord("17h34m13.46 -26d05m18.6", frame=Galactic)  # Units from string

OAO1657_41_c= SkyCoord("17h00m48.884 -41d39m21.46", frame=Galactic)  # Units from string


XTE_J1855_026_c= SkyCoord("18h55m30.4058947392 -2d36m16.737653700", frame=Galactic)  # Units from string

Cyg_X7_c= SkyCoord("20h22m00.0 40d42m00.0", frame=Galactic)  # Units from string

V395Car_c= SkyCoord("9h22m34.6761959832 -63d17m41.355773520", frame=Galactic)  # Units from string
V1341Cyg_c= SkyCoord("21h44m41.1544345272 38d19m17.066570988", frame=Galactic)  # Units from string

LZAqr_c= SkyCoord("21h23m14.54 -5d47m52.9", frame=Galactic)  # Units from string

GPVel_c= SkyCoord("9h02m06.8608812864 -40d33m16.899168060", frame=Galactic)  # Units from string


X_LMC_X4_c= SkyCoord("5h32m49.5556042296 -66d22m13.202721768", frame=Galactic)  # Units from string

Sk160_c= SkyCoord("1h17m05.1457288392 -73d26m36.014808156", frame=Galactic)  # Units from string
QVnor_c= SkyCoord("15h42m23.3633146128 -52d23m09.577395960", frame=Galactic)  # Units from string
HD153919_c= SkyCoord("17h03m56.7725629224 -37d50m38.913331452", frame=Galactic)  # Units from string

V691CrA_c= SkyCoord("18h25m46.8185061768 -37d06m18.529336656", frame=Galactic)  # Units from

V779_Cen_c= SkyCoord("11h21m15.0920532528 -60d37m25.630264596", frame=Galactic)  # Units from

HZ_Her_c= SkyCoord("16h57m49.8110126616 35d20m32.486555472", frame=Galactic)  # Units from

NGC675_c= SkyCoord("1h49m08.786 13d03m32.45", frame=Galactic)  # Units from
NGC6440_c= SkyCoord("17h48m52.67 -20d21m34.5", frame=Galactic)  # Units from
M28_c= SkyCoord("18h24m32.89 -24d52m11.4", frame=Galactic)  # Units from


wCen_c = SkyCoord("11h55m01.3260078696 -59d15m13.535219196", frame=Galactic)  # Units from string
wCen_J2140_2311B = SkyCoord("21h40m25.2 -23d11m45", frame=Galactic)  # Units from string


J1640_2224c = SkyCoord("16h40m16.745069  +22d24m08.82406", frame=Galactic)  # Units from
J2017_0603c = SkyCoord("20h17m22.705223  +06d03m05.56393", frame=Galactic)  # Units from
J2043_1711c = SkyCoord("20h43m20.8813524 +17d11m28.89999 ", frame=Galactic)  # Units from
J2302_4442c = SkyCoord("23h02m46.978387  +44d42m22.08051", frame=Galactic)  # Units from
J2317_1439c = SkyCoord("23h17m09.236381  +14d39m31.26102", frame=Galactic)  # Units from


J1750_37A_NGC6441A_c = SkyCoord("17h50m13.8016 -37d03m10.95", frame=Galactic)  # Units from
print("J1750_37A_NGC6441A_c:",J1750_37A_NGC6441A_c)


M13B_c= SkyCoord(ra=250.41829767*u.degree, dec=36.43291333*u.degree, frame='icrs')  # Units from
M13D_c= SkyCoord(ra=250.4266468*u.degree,  dec=36.45783392*u.degree, frame='icrs')  # Units from
M13E_c= SkyCoord(ra=250.42509208*u.degree, dec=36.45971322*u.degree, frame='icrs')  # Units from
M13F_c= SkyCoord(ra=250.4358575*u.degree,  dec=36.47111206*u.degree, frame='icrs')  # Units from
M13A_c= SkyCoord(ra=250.42029246*u.degree, dec=36.45416078*u.degree, frame='icrs')  # Units from
M13C_c= SkyCoord(ra=250.4208645*u.degree,  dec=36.45076217*u.degree, frame='icrs')  # Units from

J1750_37A_NGC6441A_cc= SkyCoord(ra=267.5575067	*u.degree,  dec=-37.05304167*u.degree, frame='icrs')  # Units from



J1537_5312 = SkyCoord("15h37m37.69466  -53d12m25.057", frame=Galactic)  # Units from
J1547_5709 = SkyCoord("15h47m24.1248   -57d09m17.5699", frame=Galactic)  # Units from
J1618_4624 = SkyCoord("16h18m52.77579  -46d24m34.950", frame=Galactic)  # Units from
J1727_2951 = SkyCoord("17h27m00.402    -29d51m40.8", frame=Galactic)  # Units from


J1755_25   = SkyCoord("17h55m6.0      -25d53m0",      frame=Galactic)  # Units from
J1244_6359 = SkyCoord("12h44m47.693   -63d59m47.4",   frame=Galactic)  # Units from
J1101_6424 = SkyCoord("11h01m37.1923  -64d24m39.332", frame=Galactic)  # Units from
J1614_2230 = SkyCoord("16h14m36.5051  -22d30m31.081", frame=Galactic)  # Units from


print("J1750_37A_NGC6441A_cc:",J1750_37A_NGC6441A_cc.galactic)
print("------------------------")
print("------------------------")
print("M13_B_c:",M13_B_c)
print("M13_D_c:",M13_D_c)
print("M13_E_c:",M13_E_c)
print("M13_F_c:",M13_F_c)
print("------------------------")
print("M13_A_c:",M13_A_c)
print("M13_C_c:",M13_C_c)

print("------------------------")
print("J2043_1711_c:",J2043_1711_c)
print("J1807_2459B_c:",J1807_2459B_c)
print("J1740_5340_c:",J1740_5340_c)
print("J1753_2240_c:",J1753_2240_c)
print("------------------------")
print("RX_J1856_3754_c:",RX_J1856_3754_c)

print("J1807_2500B_c:",J1807_2500B_c)
print("KS1731_260_c:",KS1731_260_c)
print("OAO1657_41_c:",OAO1657_41_c)
print("XTE_J1855_026_c:",XTE_J1855_026_c)
print("wCen_c:          ",wCen_c)
print("wCen_J2140_2311B:",wCen_J2140_2311B)


print("Cyg_X7_c:",Cyg_X7_c)

print("V395Car_c:",V395Car_c)
print("V1341Cyg_c:",V1341Cyg_c)
print("LZAqr_c:",LZAqr_c)
print("GPVel_c:",GPVel_c)
print("------------------------")
print("X_LMC_X4_c:",X_LMC_X4_c)
print("Sk160_c:",Sk160_c)
print("QVnor_c:",QVnor_c)
print("HD153919_c:",HD153919_c)

print("V691CrA_c:",V691CrA_c)
print("V779_Cen_c:",V779_Cen_c)
print("HZ_Her_c:",HZ_Her_c)

print("------------------------")
print("NGC675_c:",NGC675_c)
print("NGC6440_c:",NGC6440_c)
print("M28_c:",M28_c)


print("------------------------")
print("J1640_2224c:",J1640_2224c)
print("J2017_0603c:",J2017_0603c)
print("J2043_1711c:",J2043_1711c)
print("J2302_4442c:",J2302_4442c)
print("J2317_1439c:",J2317_1439c)
print("........................")
print("------------------------")
print("------------------------")

print("J1537_5312:",J1537_5312)
print("J1547_5709:",J1547_5709)
print("J1618_4624:",J1618_4624)
print("J1727_2951:",J1727_2951)
print("J1537_5312:",J1537_5312.galactic)
print("J1547_5709:",J1547_5709.galactic)
print("J1618_4624:",J1618_4624.galactic)
print("J1727_2951:",J1727_2951.galactic)
print("........................")
print("J1755_25:",J1755_25)
print("J1244_6359:",J1244_6359)
print("J1101_6424:",J1101_6424)
print("J1614_2230:",J1614_2230)
print("------------------------")



RX_J1856_3754_c= SkyCoord(ra=284.149375*u.degree,  dec=-37.91019444*u.degree, frame='icrs')  # Units from

J1740_5340_c= SkyCoord(ra=265.1858792*u.degree,  dec=-53.67821389*u.degree, frame='icrs')  # Units



print("M13A_c:",M13A_c.galactic)
print("M13B_c:",M13B_c.galactic)
print("M13C_c:",M13C_c.galactic)
print("M13D_c:",M13D_c.galactic)
print("M13E_c:",M13E_c.galactic)
print("M13F_c:",M13F_c.galactic)

print("RX_J1856_3754_c:",RX_J1856_3754_c.galactic)
print("J1740_5340_c:",J1740_5340_c.galactic)


NS_massErr = [
[1.52, 0.22], [1.58, 0.06], [1.96, 0.36], [1.96, 0.19], [1.02, 0.17], [1.53, 0.42], [1.41, 0.24], [1.71, 0.21], [1.21, 0.12], [1.57, 0.16], [1.74, 0.3],
[2.12, 0.16], [1.57, 0.11], [1.073, 0.36], [1.44, 0.1], [1.41, 0.1], [1.4, 0.1], [2.08, 0.19], [1.3332, 0.001], [1.35, 0.1], [1.26, 0.08], [1.4, 0.1],
[1.37, 0.13], [1.4398, 0.002], [1.4, 0.1], [2.4, 0.12], [1.358, 0.01], [1.38, 0.06], [1.77, 0.45], [1.46, 0.38], [1.41, 0.04], [1.4378, 0.0013],
[2.01, 0.04], [1.76, 0.15], [1.559, 0.005], [1.34, 0.08], [1.25, 0.05], [1.7, 0.1], [1.3381, 0.0007], [1.2489, 0.0007], [2.08, 0.1], [1.26, 0.14],
[1.71, 0.02], [1.83, 0.11], [1.7, 0.3], [1.71, 0.16], [1.35, 0.1], [1.27, 0.01], [0.86, 1.52], [1.56, 0.13], [2.5, 0.9], [1.908, 0.016],
[1.35, 0.07], [1.47, 0.07], [1.14, 0.43], [1.4, 0.56], [1.4, 0.56], [1.4, 0.56], [1.4, 0.56], [1.91, 0.02], [1.312, 0.017], [1.3384, 0.09],
[1.4, 0.56], [1.24, 0.11], [1.34, 0.1], [1.3655, 0.0021], [2.3, 1.3], [1.56, 0.24], [1.4, 0.56], [1.84, 0.11], [1.4, 0.56], [1.338, 0.002],
[1.4, 0.7], [1.666, 0.01], [1.291, 0.011], [1.48, 0.03], [1.6, 0.6], [1.33, 0.11], [1.65, 0.05], [1.29, 0.01], [1.29, 0.1], [1.832, 0.0029],
[1.34, 0.17], [1.496, 0.023], [1, 0.5], [1.51, 0.1], [1.38, 0.12], [1.33, 0.3], [1.4, 0.21], [1.8, 0.4], [1.831, 0.01], [1.393, 0.013],
[1.32, 0.1], [1.874, 0.32], [1.728, 0.066], [1.73, 0.1], [1.69, 0.1], [1.6, 0.1], [1.48, 0.1]
]

NS_mass_and_massErr = [
[1.52, 0.22, 0.18],  [1.58, 0.06, 0.06],  [1.96, 0.36, 0.36],  [1.96, 0.19, 0.19],  [1.02, 0.17, 0.17],  [1.53, 0.42, 0.42],  [1.41, 0.24, 0.24],
[1.71, 0.21, 0.21],  [1.21, 0.12, 0.12],  [1.57, 0.16, 0.16],  [1.74, 0.3, 0.3],  [2.12, 0.16, 0.16],  [1.57, 0.11, 0.11],  [1.073, 0.36, 0.36],
[1.44, 0.1, 0.1],  [1.41, 0.1, 0.1],  [1.4, 0.1, 0.1],  [2.08, 0.19, 0.19],  [1.3332, 0.001, 0.001],  [1.35, 0.1, 0.1],  [1.26, 0.08, 0.17],
[1.4, 0.1, 0.1],  [1.37, 0.13, 0.1],  [1.4398, 0.002, 0.002],  [1.4, 0.1, 0.1],  [2.4, 0.12, 0.12],  [1.358, 0.01, 0.01],  [1.38, 0.06, 0.1],
[1.77, 0.45, 0.45],  [1.46, 0.38, 0.38],  [1.41, 0.04, 0.08],  [1.4378, 0.0013, 0.0013],  [2.01, 0.04, 0.04],  [1.76, 0.15, 0.15],
[1.559, 0.005, 0.005],  [1.34, 0.08, 0.08],  [1.25, 0.05, 0.06],  [1.7, 0.1, 0.17],  [1.3381, 0.0007, 0.0007],  [1.2489, 0.0007, 0.0007],
[2.08, 0.1, 0.1],  [1.26, 0.14, 0.14],  [1.71, 0.02, 0.02],  [1.83, 0.11, 0.11],  [1.7, 0.3, 0.3],  [1.71, 0.16, 0.16],  [1.35, 0.1, 0.1],
[1.27, 0.01, 0.13],  [0.86, 1.52, 1.52],  [1.56, 0.13, 0.13],  [2.5, 0.9, 0.7],  [1.908, 0.016, 0.016],  [1.35, 0.07, 0.07],
[1.47, 0.07, 0.06],  [1.14, 0.43, 0.25],  [1.4, 0.56, 0.56],  [1.4, 0.56, 0.56],  [1.4, 0.56, 0.56],  [1.4, 0.56, 0.56],  [1.91, 0.02, 0.1],
[1.312, 0.017, 0.017],  [1.3384, 0.09, 0.09],  [1.4, 0.56, 0.56],  [1.24, 0.11, 0.11],  [1.34, 0.1, 0.1],  [1.3655, 0.0021, 0.0021],
[2.3, 1.3, 1.3],  [1.56, 0.24, 0.24],  [1.4, 0.56, 0.56],  [1.84, 0.11, 0.11],  [1.4, 0.56, 0.56],  [1.338, 0.002, 0.338],  [1.4, 0.7, 0.7],
[1.666, 0.01, 0.012],  [1.291, 0.011, 0.011],  [1.48, 0.03, 0.03],  [1.6, 0.6, 0.6],  [1.33, 0.11, 0.11],  [1.65, 0.05, 0.05],
[1.29, 0.01, 0.09],  [1.29, 0.1, 0.1],  [1.832, 0.0029, 0.0029],  [1.34, 0.17, 0.15],  [1.496, 0.023, 0.023],  [1, 0.5, 0.5],
[1.51, 0.1, 0.1],  [1.38, 0.12, 0.13],  [1.33, 0.3, 0.28],  [1.4, 0.21, 0.18],  [1.8, 0.4, 0.4],  [1.831, 0.01, 0.01],
[1.393, 0.013, 0.013],  [1.32, 0.1, 0.1],  [1.874, 0.32, 0.068],  [1.728, 0.066, 0.136],  [1.73, 0.1, 0.1],  [1.69, 0.1, 0.1],
[1.6, 0.1, 0.1],  [1.48, 0.1, 0.1]
]



l_GB_GL_dist = [
[-7.9138, 2.7882, 8.400], [48.3400, 19.8500, 1.051], [2.1220, 49.9680, 4.167], [-27.3200, 65.0300, 10.400], [-12.0210, 104.9310, 3.159], [-17.1369, 184.1245, 0.522],
[-1.1870, 168.2747, 1.900], [-35.0400, 244.5100, 12.100], [-4.5000, 245.2400, 1.100], [-4.5000, 245.2400, 1.100],
[-3.8600, 295.7900, 10.000], [54.2800, 80.8100, 0.964], [0.9500, 6.5000, 0.730],
[2.877, 9.966, 7.400], [-2.20605091, 5.83588493, 3.000], [0.4400, 12.8200, 4.419], [24.7350, 72.8300, 4.356],
[15.6120, 53.3430, 0.915], [-1.0100, 37.3400, 7.000], [0.1500, 41.6000, 7.400],
[0.1900, 45.2500, 5.600], [-16.8700, 19.9800, 2.004], [-46.0753, 62.0185, 0.268],
[-32.52909252, 276.3349498, 7.900], [-43.55935025, 300.4149551, 9.000], [-11.29080163, 356.8502712, 2.500],
[2.163693898, 327.4196998, 6.400], [37.52303328, 58.14903389, 6.600], [2.173484363, 347.754429, 1.800],
[0.335517538, 292.0903507, 5.700], [0.335517538, 292.0903507, 5.700], [-0.3539, 351.4972, 7.500], [-43.804, 303.514, 59.700], [-32.52909252, 276.3349498, 7.900],
[-43.55935025, 300.4149551, 9.000], [3.929852619, 263.0582983, 1.900],
[-0.8505, 330.9263, 3.600], [37.52303328, 58.14903389, 6.600], [2.173484363, 347.754429, 1.800], [-11.29080163, 356.8502712, 2.000], [75.4140, 311.3100, 0.709],
[15.9600, 350.9760, 1.800], [-0.5970, 21.5870, 4.685], [3.5010, 55.7770, 0.300], [2.1554, 78.4883, 6.100],
[-19.8109, 279.9777, 7.100], [-0.0470, 359.9440, 8.300], [-0.1750, 359.7880, 8.109],
[-0.0200, 0.1500, 8.217], [-0.2330, 0.1260, 8.149], [-0.0700, 6.8400, 4.000], [16.8060, 44.6420, 2.361], [0.3870, 12.9040, 4.492], [-1.0010, 16.8050, 4.440],
[3.65256, 1.07299, 8.000], [40.91351866, 58.99968688, 7.100], [40.91268749, 58.99528387, 7.100], [0.319238, 344.369163, 6.400], [-17.21600527, 358.5988662, 0.167],
[2.830882823, 295.7666832, 8.178], [-2.092058679, 31.07638514, 10.000],
[-9.33975534, 281.8354681, 10.000], [-9.33975534, 281.8354681, 10.000], [6.7700, 20.7900, 7.800], [3.0600, 42.2900, 1.200], [-4.697, 59.197, 1.400], [-11.31636463, 87.32819642, 8.000], [-11.31636463, 87.32819642, 8.000], [-30.0400, 169.9900, 1.300], [-36.7736, 183.3368, 2.100], [-41.9600, 253.3900, 0.157],
[-2.0100, 200.5700, 0.420], [29.5990, 149.7300, 1.136], [21.0900, 202.7300, 1.110],
[-5.7440, 283.6680, 4.04], [50.8600, 160.3500, 0.700], [51.1000, 231.7900, 0.645], [45.7800, 243.4900, 1.370], [12.2540, 280.8510, 0.340], [13.8000, 298.9700, 1.613],
[16.451, 344.09, 1.887], [20.1900, 352.6400, 0.700], [25.2200, 28.7500, 1.311], [17.7400, 27.7200, 1.471], [-11.9669, 338.1647, 2.400], [21.6410, 37.8850, 1.667],
[1.6900, 3.8400, 6.900], [1.663, 6.299, 3.232], [-2.2060, 5.8360, 3.000], [5.367, 44.875, 2.083], [-19.6000, 359.7300, 1.140],
[1.795, 46.564, 1.496], [-25.7300, 336.5250, 4.000], [-9.1200, 30.0300, 1.111], [4.7100, 69.2900, 6.940], [2.5536, 66.8583, 6.500], [-1.1687, 61.0975, 7.269], [-8.675, 60.522, 2.162],
[-6.6200, 64.7500, 1.163], [-15.3100, 61.9200, 1.600], [-3.926, 77.832, 5.100], [1.302, 86.861, 4.120], [-42.0800, 47.7800, 0.714], [-43.006, 72.991, 0.971], [40.9127738, 58.97150791, 7.100], [40.90884273, 59.00526155, 7.100], [40.91029222, 59.00755123, 7.100], [40.90293961, 59.02380103, 7.100], [3.929852619, 263.0582983, 1.900], [-36.19939, 46.48299, 8.500], [46.8060, 3.8590, 7.500],
[-44.9020, 305.8960, 4.690], [-44.9030, 305.9000, 4.900], [3.80168, 7.72912, 8.500], [0.6110, 8.3820, 0.760], [-5.58068, 7.79821, 5.600], [-47.44257, 143.961909, 3.400],
[1.6700, 3.8100, 6.900], [1.6700, 3.8100, 6.900], [1.6700, 3.8100, 6.900], [1.6700, 3.8100, 6.900], [1.6700, 3.8100, 6.900], [1.6700, 3.8100, 6.900], [1.6700, 3.8100, 6.900]
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
[272.979308, -17.61047, 4.419],
[274.1497265, 45.1760727, 4.356],
[277.394445, 24.9383869, 0.915],
[285.7741384, 3.455335864, 7.000],
[286.70358, 7.77386, 7.400],
[288.3710592, 11.034928, 5.600],
[292.623815, -18.862853, 2.004],
[335.524871, -1.62103475, 0.268],
[83.20648168, -66.37033409, 7.900],
[19.27144054, -73.44333745, 9.000],
[276.4450771, -37.10514704, 2.500],
[235.5973471, -52.38599372, 6.400],
[254.4575459, 35.34235738, 6.600],
[255.9865524, -37.84414259, 1.800],
[170.3128836, -60.62378618, 5.700],
[170.3128836, -60.62378618, 5.700],
[261.297496, -36.282665, 7.500],
[11.3965, -73.3175, 59.700],
[83.20648168, -66.37033409, 7.900],
[19.27144054, -73.44333745, 9.000],
[135.528587, -40.55469421, 1.900],
[243.1793333, -52.4232222, 3.600],
[254.4575459, 35.34235738, 6.600],
[255.9865524, -37.84414259, 1.800],
[276.4450771, -37.10514704, 2.000],
[195.0149029, 12.68235336, 0.709],
[245.9092575, -26.5316025, 1.800],
[278.170278, -10.3591, 4.685],
[290.436729, 21.883958, 0.300],
[305.5, 40.7, 6.100],
[117.14084, -67.75153, 7.100],
[266.417359, -29.008304, 8.300],
[266.4492935, -29.20855, 8.109],
[266.513989, 28.8370514, 8.217],
[266.7077283, 28.9497194, 8.149],
[270.332562, -23.079066, 4.000],
[272.65533, 17.743717, 2.361],
[273.0660802, -17.5605197, 4.492],
[276.2623, -14.7816, 4.440],
[263.5560833, -26.0885, 8.000],
[250.4202925, 36.45416078, 7.100],
[250.4208645, 36.45076217, 7.100],
[255.2036833, -41.65596111, 6.400],
[284.149375, -37.91019444, 0.167],
[178.755525, -59.25375978, 8.178],
[283.8766912, -2.60464935, 10.000],
[140.6444842, -63.29482105, 10.000],
[140.6444842, -63.29482105, 10.000],
[271.207898, -7.59019, 7.800],
[284.4016262, 9.721442158, 1.200],
[299.9032078, 20.80420061, 1.400],
[326.1714768, 38.32140738, 8.000],
[326.1714768, 38.32140738, 8.000],
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
[310.8370056, 17.19136111, 1.600],
[311.2562713, 36.55039017, 5.100],
[313.4692835, 46.84769947, 4.120],
[326.460248, -7.83847452, 0.714],
[338.5961379, 6.191301758, 0.971],
[250.4182977, 36.43291333, 7.100],
[250.4266468, 36.45783392, 7.100],
[250.4250921, 36.45971322, 7.100],
[250.4358575, 36.47111206, 7.100],
[135.528587, -40.55469421, 1.900],
[320.8105833, -5.79802778, 8.500],
[229.631075, 2.087631, 7.500],
[6.02793, -72.06855742, 4.690],
[14.4333, -72.0219, 4.900],
[267.2194583, -20.35958333, 8.500],
[270.5222316, -21.4010136, 0.760],
[276.1370417, -24.86983333, 5.600],
[27.28660833, 13.05901389, 3.400],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900],
[267.02, -24.77917, 6.900]
]

NS_J2000_name = [
['NS-NS binaries' , '4U 1820-30'],
['NS-NS binaries' , 'B1534+12'],
['NS-NS binaries' , 'B1913+16'],
['NS-NS binaries' , 'B2127+11C'],
['NS-NS binaries' , 'B2303+46'],
['NS-NS binaries' , 'J0453+1559'],
['NS-NS binaries' , 'J0509+3801'],
['NS-NS binaries' , 'J0514-4002A / NGC 1851'],
['NS-NS binaries' , 'J0737-3039A'],
['NS-NS binaries' , 'J0737-3039B'],
['NS-NS binaries' , 'J1141-6545'],
['NS-NS binaries' , 'J1518+4904'],
['NS-NS binaries' , 'J1756-2251'],
['NS-NS binaries' , 'J1757-1854'],
['NS-NS binaries' , 'J1807-2500B'],
['NS-NS binaries' , 'J1811-1736'],
['NS-NS binaries' , 'J1816+4510'],
['NS-NS binaries' , 'J1829+2456'],
['NS-NS binaries' , 'J1903+0327'],
['NS-NS binaries' , 'J1906+0746'],
['NS-NS binaries' , 'J1913+1102'],
['NS-NS binaries' , 'J1930-1852'],
['NS-NS binaries' , 'J2222-0137'],
['NS-NS binaries' , 'LMC X-4'],
['NS-NS binaries' , 'Sk 160/SMC X-1'],
['NS in X-ray binaries' , '2A 1822-371/V* V691 CrA'],
['NS in X-ray binaries' , '4U 1538-522/V* QV Nor'],
['NS in X-ray binaries' , '4U 1656+35/Her X-1'],
['NS in X-ray binaries' , '4U 1700-377/HD 153919'],
['NS in X-ray binaries' , 'Cen X-3/V* V779 Cen'],
['NS in X-ray binaries' , 'Cen X-3/V* V779 Cen'],
['NS in X-ray binaries' , 'EXO 1722-363'],
['NS in X-ray binaries' , 'J0045-7319'],
['NS in X-ray binaries' , 'LMC X-4'],
['NS in X-ray binaries' , 'Sk 160/SMC X-1'],
['NS in X-ray binaries' , 'Vela X-1/V* GP Vel'],
['radio millisecond pulsars' , '4U 1608-522'],
['radio millisecond pulsars' , '4U 1656+35/Her X-1'],
['radio millisecond pulsars' , '4U 1700-377/HD 153919'],
['radio millisecond pulsars' , '4U 1822-371/V* V691 CrA'],
['radio millisecond pulsars' , 'B1257+12'],
['radio millisecond pulsars' , 'B1620-26'],
['radio millisecond pulsars' , 'B1829-10'],
['radio millisecond pulsars' , 'B1919+21'],
['radio millisecond pulsars' , 'Cyg X-7'],
['radio millisecond pulsars' , 'EXO 0748-676'],
['radio millisecond pulsars' , 'J1745-2900'],
['radio millisecond pulsars' , 'J1745-2912'],
['radio millisecond pulsars' , 'J1746-2849'],
['radio millisecond pulsars' , 'J1746-2856'],
['radio millisecond pulsars' , 'J1801-2304'],
['radio millisecond pulsars' , 'J1810+1744'],
['radio millisecond pulsars' , 'J1812-1733'],
['radio millisecond pulsars' , 'J1825-1446'],
['radio millisecond pulsars' , 'KS 1731-260'],
['radio millisecond pulsars' , 'M13 A'],
['radio millisecond pulsars' , 'M13 C'],
['radio millisecond pulsars' , 'OAO 1657-415'],
['radio millisecond pulsars' , 'RX J1856–3754'],
['radio millisecond pulsars' , 'w Cen'],
['radio millisecond pulsars' , 'XTE J1855-026'],
['WD-NS binaries' , '2S 0921-630/V* V395 Car'],
['WD-NS binaries' , '2S 0921-630/V* V395 Car'],
['WD-NS binaries' , 'B1802-07'],
['WD-NS binaries' , 'B1855+09'],
['WD-NS binaries' , 'B1957+20'],
['WD-NS binaries' , 'Cyg X-2/V* V1341 Cyg'],
['WD-NS binaries' , 'Cyg X-2/V* V1341 Cyg'],
['WD-NS binaries' , 'J0337+1715'],
['WD-NS binaries' , 'J0348+0432'],
['WD-NS binaries' , 'J0437-4715'],
['WD-NS binaries' , 'J0621+1002'],
['WD-NS binaries' , 'J0740+6620'],
['WD-NS binaries' , 'J0751+1807'],
['WD-NS binaries' , 'J0955-6150'],
['WD-NS binaries' , 'J1012+5307'],
['WD-NS binaries' , 'J1022+1001'],
['WD-NS binaries' , 'J1023+0038'],
['WD-NS binaries' , 'J1045-4509'],
['WD-NS binaries' , 'J1227-4853'],
['WD-NS binaries' , 'J1600-3053'],
['WD-NS binaries' , 'J1614-2230'],
['WD-NS binaries' , 'J1713+0747'],
['WD-NS binaries' , 'J1738+0333'],
['WD-NS binaries' , 'J1740-5340'],
['WD-NS binaries' , 'J1741+1351'],
['WD-NS binaries' , 'J1748-2446I'],
['WD-NS binaries' , 'J1753-2240'],
['WD-NS binaries' , 'J1807-2459B'],
['WD-NS binaries' , 'J1853+1303'],
['WD-NS binaries' , 'J1909-3744'],
['WD-NS binaries' , 'J1910+1256'],
['WD-NS binaries' , 'J1910-5959A'],
['WD-NS binaries' , 'J1918-0642'],
['WD-NS binaries' , 'J1946+3417'],
['WD-NS binaries' , 'J1949+3106'],
['WD-NS binaries' , 'J1950+2414'],
['WD-NS binaries' , 'J2016+1948'],
['WD-NS binaries' , 'J2019+2425'],
['WD-NS binaries' , 'J2043+1711'],
['WD-NS binaries' , 'J2045+3633'],
['WD-NS binaries' , 'J2053+4650'],
['WD-NS binaries' , 'J2145-0750'],
['WD-NS binaries' , 'J2234+0611'],
['WD-NS binaries' , 'M13 B'],
['WD-NS binaries' , 'M13 D'],
['WD-NS binaries' , 'M13 E'],
['WD-NS binaries' , 'M13 F'],
['WD-NS binaries' , 'Vela X-1/V* GP Vel'],
['WD-NS binaries' , 'XTE J2123-058/V* LZ Aqr '],
['WD-NS binaries GalCluster pulsar' , 'B1516+02B/J1518+0204B (M5)'],
['WD-NS binaries GalCluster pulsar' , 'J0024-7204H/B0021-72H (47 Tucanae)'],
['WD-NS binaries GalCluster pulsar' , 'J0024-7204I/B0021-72I (47 Tucanae)'],
['WD-NS binaries GalCluster pulsar' , 'J1748-2021B / NGC 6440'],
['WD-NS binaries GalCluster pulsar' , 'J1802-2124'],
['WD-NS binaries GalCluster pulsar' , 'J1824-2452C / M28'],
['WD-NS binaries GalCluster pulsar' , 'J1911-5958A / NGC 6752'],
['WD-NS binaries GalCluster pulsar' , 'Ter 5ai / J1748-2446ai'],
['WD-NS binaries GalCluster pulsar' , 'Ter 5I / J1748-2446I'],
['WD-NS binaries GalCluster pulsar' , 'Ter 5J / J1748-2446J'],
['WD-NS binaries GalCluster pulsar' , 'Ter 5U / J1748-2446U'],
['WD-NS binaries GalCluster pulsar' , 'Ter 5W / J1748-2446W'],
['WD-NS binaries GalCluster pulsar' , 'Ter 5X / J1748-2446X'],
['WD-NS binaries GalCluster pulsar' , 'Ter 5Z / J1748-2446Z']
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
[1.27, 0.01, 0.13, 1.01, 0.001, 0.001],
[1.56, 0.13, 0.13, 1.05, 0.45, 0.45],
[1.312, 0.017, 0.017, 1.258, 0.017, 0.017],
[1.3384, 0.09, 0.09, 1.3946, 0.09, 0.09],
[1.3655, 0.0021, 0.0021, 1.2068, 0.0022, 0.0016],
[1.56, 0.24, 0.24, 1.18, 0.03, 0.03],
[1.84, 0.11, 0.11, 0.193, 0.012, 0.012],
[1.338, 0.002, 0.338, 1.256, 0.346, 0.003],
[1.666, 0.01, 0.012, 1.033, 0.011, 0.008],
[1.291, 0.011, 0.011, 1.322, 0.011, 0.011],
[1.65, 0.05, 0.05, 1.24, 0.05, 0.05],
[1.29, 0.1, 0.1, 1.3, 0.1, 0.1],
[1.831, 0.01, 0.01, 1.3194, 0.04, 0.04],
[1.57, 0.11, 0.11, 1.29, 0.05, 0.05],
[1.21, 0.12, 0.12, 1.04, 0.09, 0.09],
[0.97, 0.24, 0.24, 0.33, 0.05, 0.05],
[1.02, 0.17, 0.17, 16.4, 5.2, 4],
[1.5, 0.3, 0.3, 2.3, 0.3, 0.3],
[2.44, 0.27, 0.27, 58, 11, 11],
[1.57, 0.16, 0.16, 19.7, 4.3, 4.3],
[1.24, 0.24, 0.24, 19.7, 4.3, 4.3],
[1.46, 0.38, 0.38, 13.6, 1.6, 1.6],
[1.58, 0.34, 0.34, 8.8, 1.8, 1.8],
[1.31, 0.14, 0.14, 15.6, 1.8, 1.8],
[1.05, 0.09, 0.09, 15.5, 1.5, 1.5],
[1.88, 0.13, 0.13, 23.1, 0.2, 0.2],
[1.52, 0.22, 0.18, 0, 0, 0],
[1.073, 0.36, 0.36, 0, 0, 0],
[1.96, 0.19, 0.19, 0, 0, 0],
[1.96, 0.36, 0.36, 0, 0, 0],
[1.4, 0.1, 0.1, 0, 0, 0],
[1.35, 0.1, 0.1, 0, 0, 0],
[1.4, 0.1, 0.1, 0, 0, 0],
[1.4, 0.1, 0.1, 0, 0, 0],
[0.832, 1.19, 0.0551, 0, 0, 0],
[1.77, 0.45, 0.45, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[2.3, 1.3, 1.3, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.4, 0.56, 0.56, 0, 0, 0],
[1.3, 0.5, 0.5, 0, 0, 0],
[1.48, 0.21, 0.64, 0, 0, 0],
[1.48, 0.21, 0.64, 0, 0, 0],
[1.74, 0.3, 0.3, 16, 2, 2],
[1.7, 0.3, 0.3, 0, 0, 0],
[1.43, 0.26, 0.61, 0, 0, 0],
[1.41, 0.24, 0.24, 0, 0, 0],
[1.44, 0.1, 0.1, 0.35, 0.03, 0.03],
[1.44, 0.1, 0.1, 0.35, 0.03, 0.03],
[1.26, 0.08, 0.17, 0.36, 0.67, 0.15],
[1.37, 0.13, 0.1, 0.244, 0.014, 0.012],
[2.4, 0.12, 0.12, 0.035, 0.002, 0.002],
[1.78, 0.23, 0.23, 0.6, 0.13, 0.13],
[1.5, 0.3, 0.3, 0.63, 0.16, 0.16],
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
[1.29, 0.01, 0.09, 0.231, 0.01, 0.01],
[1.832, 0.0029, 0.0029, 0.2659, 0.003, 0.003],
[1.34, 0.17, 0.15, 0.81, 0.06, 0.05],
[1.496, 0.023, 0.023, 0.28, 0.005, 0.005],
[1, 0.5, 0.5, 0.43, 0.5, 0.5],
[1.51, 0.1, 0.1, 0.33, 0.1, 0.1],
[1.38, 0.12, 0.13, 0.173, 0.01, 0.01],
[1.33, 0.3, 0.28, 0.94, 0.14, 0.13],
[1.4, 0.21, 0.18, 0.86, 0.07, 0.06],
[1.8, 0.4, 0.4, 0.9, 0.05, 0.05],
[1.353, 0.014, 0.017, 0.298, 0.15, 0.12],
[1.48, 0.21, 0.64, 0.1876, 0.1, 0.1],
[1.48, 0.21, 0.64, 0.2085, 0.1, 0.1],
[1.48, 0.21, 0.64, 0.0225, 0.1, 0.1],
[1.48, 0.21, 0.64, 0.1517, 0.1, 0.1],
[1.77, 0.08, 0.08, 0.35, 0.1, 0.1],
[1.46, 0.3, 0.39, 0.53, 0.28, 0.39],
[2.08, 0.19, 0.19, 0.21, 0.1, 0.1],
[1.61, 0.04, 0.04, 0.18, 0.086, 0.016],
[1.41, 0.1, 0.1, 0.15, 0.1, 0.1],
[2.62, 0.223606798, 0.223606798, 0.12, 0.1, 0.1],
[1.24, 0.11, 0.11, 0.78, 0.04, 0.04],
[1.616, 0.007, 0.007, 0.27, 0.1, 0.1],
[1.4, 0.16, 0.1, 0.18, 0.02, 0.02],
[1.32, 0.1, 0.1, 0.49, 0.1, 0.1],
[1.874, 0.32, 0.068, 0.25, 0.1, 0.1],
[1.728, 0.066, 0.136, 0.35, 0.1, 0.1],
[1.73, 0.1, 0.1, 0.39, 0.1, 0.1],
[1.69, 0.1, 0.1, 0.25, 0.1, 0.1],
[1.6, 0.1, 0.1, 0.25, 0.1, 0.1],
[1.48, 0.1, 0.1, 0.22, 0.1, 0.1]
]

Baryonic_freq=[
26.3821, 16.9405, 32.7554, 0.9378, 21.8427, 13.0648, 200.3777, 44.0541, 0.3606, 2.5387, 24.4290, 35.1351, 46.5176, 238.8810, 9.5986, 313.1749, 24.3844, 465.1352, 6.9409, 36.6502, 5.3902, 30.4712, 1.079591939, 160.8097, 90.2873, 3.0271, 0.7478, 0.2657, 5.3368, 0.6764, 1.0579, 2.4048, 602.4096, 1.8576, 3.5817, 96.36223457, 268.6668662, 43.2884, 186.4941, 622.1220305, 365.9534, 25.5606, 173.6879, 34.6574, 346.5320, 287.4579, 500.2501, 190.2678, 60.7794, 592.4215, 133.7931, 592.9878, 277.9377, 317.3789, 218.8118, 170.9374, 266.8692, 104.4911, 10.51106825, 238.8814, 244.3913777, 339.3157, 200.6588053, 306.1674, 130.7895, 315.4436, 76.1140, 232.3002, 15.39873763, 254.1603, 420.1894, 31.5638, 79.4516, 62.2959, 279.5966, 283.4409432, 320.6886957, 402.0937023, 332.9448049, 125.8346, 311.4934, 59.66541822, 8.960506248, 79.0664, 240.442414, 47.1068, 104.4911, 12.4474, 304.0308, 237.8019, 333.4156, 406.0765
]

Baryonic_freq_NSmass=[
[26.3821, 1.3332], [16.9405, 1.4398], [32.7554, 1.358], [0.9378, 1.38], [21.8427, 1.559], [13.0648, 1.34], [200.3777, 1.25], [44.0541, 1.3381], [0.3606, 1.2489], [2.5387, 1.27], [24.4290, 1.56], [35.1351, 1.312], [46.5176, 1.3384], [238.8810, 1.3655], [9.5986, 1.56], [313.1749, 1.84], [24.3844, 1.338], [465.1352, 1.666], [6.9409, 1.291], [36.6502, 1.65], [5.3902, 1.29], [30.4712, 1.831], [1.079591939, 1.58], [160.8097, 1.4], [90.2873, 1.35], [3.0271, 1.4], [0.7478, 1.4], [0.2657, 1.4], [5.3368, 1.4], [0.6764, 1.4], [1.0579, 1.4], [2.4048, 1.4], [602.4096, 2.3], [1.8576, 1.4], [3.5817, 1.4], [96.36223457, 1.48], [268.6668662, 1.48], [43.2884, 1.26], [186.4941, 1.37], [622.1220305, 2.4], [365.9534, 1.4378], [25.5606, 2.01], [173.6879, 1.76], [34.6574, 1.7], [346.5320, 2.08], [287.4579, 1.26], [500.2501, 1.71], [190.2678, 1.83], [60.7794, 1.7], [592.4215, 1.761], [133.7931, 1.35], [592.9878, 0.86], [277.9377, 2.5], [317.3789, 1.908], [218.8118, 1.35],
[170.9374, 1.47], [266.8692, 1.14], [104.4911, 1.91], [10.51106825, 1.25], [238.8814, 1.34], [244.3913777, 1.4], [339.3157, 1.48], [200.6588053, 1.6], [306.1674, 1.33], [130.7895, 1.29], [315.4436, 1.832], [76.1140, 1.34], [232.3002, 1.496], [15.39873763, 1], [254.1603, 1.51], [420.1894, 1.38], [31.5638, 1.33], [79.4516, 1.4], [62.2959, 1.8], [279.5966, 1.353], [283.4409432, 1.48], [320.6886957, 1.48], [402.0937023, 1.48], [332.9448049, 1.48], [125.8346, 2.08], [311.4934, 1.61], [59.66541822, 2.548], [8.960506248, 1.97], [79.0664, 1.24], [240.442414, 1.616], [47.1068, 1.32], [104.4911, 1.874], [12.4474, 1.728], [304.0308, 1.73], [237.8019, 1.69], [333.4156, 1.6], [406.0765, 1.48]
]

def column(array, i):
    return [row[i] for row in array]

print("NS_massErr:",NS_massErr[1])
print("NS_massErr:",NS_massErr[0:1])
print("NS_massErr:",NS_massErr[0:2])
print("NS_massErr:",NS_massErr[1:1])
print("NS_massErr:", column(NS_massErr,0))
print("NS_massErr:", column(NS_massErr,1))
print("NS_massErr:", column(NS_massErr,0)[0],column(NS_massErr,1)[0])
print("NS_massErr:", column(NS_massErr,0)[1],column(NS_massErr,1)[1])

print("NS_massErr:", NS_mass_and_massErr[0])

print("NS_massErr:", column(NS_mass_and_massErr,0)[0],column(NS_mass_and_massErr,1)[0],column(NS_mass_and_massErr,2)[0])
print("NS_massErr:", column(NS_mass_and_massErr,0)[1],column(NS_mass_and_massErr,1)[1],column(NS_mass_and_massErr,2)[1])
print("NS_massErr:", column(NS_mass_and_massErr,0)[2],column(NS_mass_and_massErr,1)[2],column(NS_mass_and_massErr,2)[2])

print("NS_massErr:", len(l_GB_GL_dist))

print("NS_massErr:", column(l_GB_GL_dist,0)[0],column(l_GB_GL_dist,1)[0],column(l_GB_GL_dist,2)[0])
print("NS_massErr:", column(l_GB_GL_dist,0)[1],column(l_GB_GL_dist,1)[1],column(l_GB_GL_dist,2)[1])
print("NS_massErr:", column(l_GB_GL_dist,0)[2],column(l_GB_GL_dist,1)[2],column(l_GB_GL_dist,2)[2])

print("NS_J2000_name:", column(NS_J2000_name,0)[2],column(NS_J2000_name,1)[2])
print("NS_mass:",
      column(NS_and_compStar_mass_and_unc,0)[2],
      "+/-",
      column(NS_and_compStar_mass_and_unc,1)[2],
      column(NS_and_compStar_mass_and_unc,2)[2],
      "compStar_mass:",
      column(NS_and_compStar_mass_and_unc,3)[2],
      "+/-",
      column(NS_and_compStar_mass_and_unc,4)[2],
      column(NS_and_compStar_mass_and_unc,5)[2])

print("GB:", column(l_GB_GL_dist,0))
print("GL:", column(l_GB_GL_dist,1))


NSbinaries = []
for i in range(len(NS_J2000_name)):
    if column(NS_J2000_name,0)[i] == 'NS-NS binaries':
        NSbinaries.append(column(NS_and_compStar_mass_and_unc,3)[i])

WDcompanions = []
for i in range(len(NS_J2000_name)):
    if column(NS_J2000_name,0)[i] == 'WD-NS binaries' or column(NS_J2000_name,0)[i] == 'WD-NS binaries GalCluster pulsar':
        WDcompanions.append(column(NS_and_compStar_mass_and_unc,3)[i])

NS_isolated = []
for i in range(len(NS_J2000_name)):
    if column(NS_J2000_name,0)[i] == 'radio millisecond pulsars':
        NS_isolated.append(column(NS_and_compStar_mass_and_unc,0)[i])

NSNS_binaries = []
for i in range(len(NS_J2000_name)):
    if column(NS_J2000_name,0)[i] == 'NS-NS binaries':
        NSNS_binaries.append(column(NS_and_compStar_mass_and_unc,0)[i])

WDNS_binaries = []
for i in range(len(NS_J2000_name)):
    if column(NS_J2000_name,0)[i] == 'WD-NS binaries':
        WDNS_binaries.append(column(NS_and_compStar_mass_and_unc,0)[i])

Xray_binaries = []
for i in range(len(NS_J2000_name)):
    if column(NS_J2000_name,0)[i] == 'NS in X-ray binaries':
        Xray_binaries.append(column(NS_and_compStar_mass_and_unc,0)[i])


print("NS masses:", column(NS_and_compStar_mass_and_unc,0))
print("NS masses:", len(column(NS_and_compStar_mass_and_unc,0)))

print("NS binaries masses:", NSbinaries)
print("NS binaries masses:", len(NSbinaries))


NS_array1 = np.array(column(NS_and_compStar_mass_and_unc,0), float)
NS_array2 = np.array(NSbinaries, float)

print("NS_array1:", NS_array1)
print("NS_array2:", NS_array2)
print("length NS_array1:", len(NS_array1))
print("length NS_array2:", len(NS_array2))

ns_masses_extra = np.concatenate((NS_array1, NS_array2))
print("NS binaries masses extra:", ns_masses_extra)
print("NS binaries masses extra:", len(ns_masses_extra))

print("length NS_isolated:", len(NS_isolated))
print("length NSNS_binaries:", len(NSNS_binaries))
print("length WDNS_binaries:", len(WDNS_binaries))
print("length Xray_binaries:", len(Xray_binaries))


#ns_masses=column(NS_and_compStar_mass_and_unc,0)
ns_masses=ns_masses_extra


#n, bins, patches = plt.hist(ns_masses, 40, density=True, facecolor='g', alpha=0.75)
n, bins, patches = plt.hist(ns_masses, 40, density=False, facecolor='g', alpha=0.75)


#plt.set_size_inches(18.5, 10.5)
plt.xlabel('NS mass [$M_{sun}$]')
plt.ylabel('#Entries')
#plt.title('Histogram of IQ')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(0.5, 3)
plt.ylim(0, 26)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(5, 7)
plt.savefig('plots/NS_mass_dist.png', dpi=200, transparent=True)
plt.savefig('plots/NS_mass_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/NS_mass_dist.eps', dpi=200, transparent=True)
plt.clf()



#n, bins, patches = plt.hist(ns_masses, 40, density=True, facecolor='g', alpha=0.75)
n, bins, patches = plt.hist(WDcompanions, 40, density=False, facecolor='g', alpha=0.75)

#plt.set_size_inches(18.5, 10.5)
plt.xlabel('companion WD mass [$M_{sun}$]')
plt.ylabel('#Entries')
plt.title('WD-NS system')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(0, 1.5)
plt.ylim(0, 10)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(5, 7)
plt.savefig('plots/WDNS_mass_dist.png', dpi=200, transparent=True)
plt.savefig('plots/WDNS_mass_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/WDNS_mass_dist.eps', dpi=200, transparent=True)
plt.clf()



n, bins, patches = plt.hist(NS_isolated, 20, density=False, facecolor='g', alpha=0.75)

#plt.set_size_inches(18.5, 10.5)
plt.xlabel('NS mass [$M_{sun}$]')
plt.ylabel('#Entries')
plt.title('NS isolated')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(0.5, 3)
plt.ylim(0, 15)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(5, 7)
plt.savefig('plots/NSisolated_mass_dist.png', dpi=200, transparent=True)
plt.savefig('plots/NSisolated_mass_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/NSisolated_mass_dist.eps', dpi=200, transparent=True)
plt.clf()

n, bins, patches = plt.hist(NSNS_binaries, 20, density=False, facecolor='g', alpha=0.75)

#plt.set_size_inches(18.5, 10.5)
plt.xlabel('NS mass [$M_{sun}$]')
plt.ylabel('#Entries')
plt.title('NS-NS binaries')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(0.5, 3)
plt.ylim(0, 8)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(5, 7)
plt.savefig('plots/NSNSbinaries_mass_dist.png', dpi=200, transparent=True)
plt.savefig('plots/NSNSbinaries_mass_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/NSNSbinaries_mass_dist.eps', dpi=200, transparent=True)
plt.clf()


n, bins, patches = plt.hist(WDNS_binaries, 20, density=False, facecolor='g', alpha=0.75)

#plt.set_size_inches(18.5, 10.5)
plt.xlabel('NS mass [$M_{sun}$]')
plt.ylabel('#Entries')
plt.title('WD-NS binaries')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(0.5, 3)
plt.ylim(0, 10)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(5, 7)
plt.savefig('plots/WDNSbinaries_mass_dist.png', dpi=200, transparent=True)
plt.savefig('plots/WDNSbinaries_mass_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/WDNSbinaries_mass_dist.eps', dpi=200, transparent=True)
plt.clf()

n, bins, patches = plt.hist(Xray_binaries, 20, density=False, facecolor='g', alpha=0.75)

#plt.set_size_inches(18.5, 10.5)
plt.xlabel('NS mass [$M_{sun}$]')
plt.ylabel('#Entries')
plt.title('X-ray binaries')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(0.5, 3)
plt.ylim(0, 5)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(5, 7)
plt.savefig('plots/XRaybinaries_mass_dist.png', dpi=200, transparent=True)
plt.savefig('plots/XRaybinaries_mass_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/XRaybinaries_mass_dist.eps', dpi=200, transparent=True)
plt.clf()



n, bins, patches = plt.hist(column(l_GB_GL_dist,0), 40, density=False, facecolor='g', alpha=0.75)

plt.xlabel('NS GB [deg]')
plt.ylabel('#Entries')
plt.title('Galactic coordinates of NS')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(-90, 90)
plt.ylim(0, 40)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(14, 7)
plt.savefig('plots/NS_GB_dist.png', dpi=200, transparent=True)
plt.savefig('plots/NS_GB_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/NS_GB_dist.eps', dpi=200, transparent=True)
plt.clf()


n, bins, patches = plt.hist(column(l_GB_GL_dist,1), 40, density=False, facecolor='g', alpha=0.75)

plt.xlabel('NS GL [deg]')
plt.ylabel('#Entries')
plt.title('Galactic coordinates of NS')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(0, 360)
plt.ylim(0, 30)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(14, 7)
plt.savefig('plots/NS_GL_dist.png', dpi=200, transparent=True)
plt.savefig('plots/NS_GL_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/NS_GL_dist.eps', dpi=200, transparent=True)
plt.clf()



n, bins, patches = plt.hist(column(RA_DEC_deg_dist_kpc,0), 40, density=False, facecolor='g', alpha=0.75)

plt.xlabel('NS ra [deg]')
plt.ylabel('#Entries')
plt.title('Galactic coordinates of NS')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(0, 360)
plt.ylim(0, 40)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(14, 7)
plt.savefig('plots/NS_RA_dist.png', dpi=200, transparent=True)
plt.savefig('plots/NS_RA_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/NS_RA_dist.eps', dpi=200, transparent=True)
plt.clf()


n, bins, patches = plt.hist(column(RA_DEC_deg_dist_kpc,1), 40, density=False, facecolor='g', alpha=0.75)

plt.xlabel('NS dec [deg]')
plt.ylabel('#Entries')
plt.title('Galactic coordinates of NS')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(-80, 80)
plt.ylim(0, 30)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(14, 7)
plt.savefig('plots/NS_DEC_dist.png', dpi=200, transparent=True)
plt.savefig('plots/NS_DEC_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/NS_DEC_dist.eps', dpi=200, transparent=True)
plt.clf()


n, bins, patches = plt.hist(Baryonic_freq, 50, density=False, facecolor='g', alpha=0.75)

plt.xlabel('NS barycentric rotation frequency, $F_{0}$ [Hz]')
plt.ylabel('#Entries')
#plt.title('Galactic coordinates of NS')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.xlim(0, 500)
plt.ylim(0, 30)
plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
fig.set_size_inches(14, 7)
plt.savefig('plots/NS_baryonic_freq_dist.png', dpi=200, transparent=True)
plt.savefig('plots/NS_baryonic_freq_dist.pdf', dpi=200, transparent=True)
plt.savefig('plots/NS_baryonic_freq_dist.eps', dpi=200, transparent=True)
plt.clf()





#cmap=plt.get_cmap('inferno')
#cmap=plt.get_cmap('plasma')
#cmap=plt.get_cmap('YlOrBr')
#cmap=plt.get_cmap('YlGn')
#cmap=plt.get_cmap('tab20b')
#cmap=plt.get_cmap('gist_ncar')
#cmap=plt.get_cmap('nipy_spectral')
#cmap=plt.get_cmap('CMRmap')

#x_bins = np.linspace(0, 500, 100)
#y_bins = np.linspace(0, 3,  40)

x_min = np.min(column(Baryonic_freq_NSmass,0))
x_max = np.max(column(Baryonic_freq_NSmass,0))

y_min = np.min(column(Baryonic_freq_NSmass,1))
y_max = np.max(column(Baryonic_freq_NSmass,1))

x_bins = np.linspace(x_min, x_max, 25)
y_bins = np.linspace(y_min, y_max, 20)

print("length NS_isolated:", len(Baryonic_freq_NSmass))
print("length Baryonic_freq_NSmass:", column(Baryonic_freq_NSmass,0))
print("length Baryonic_freq_NSmass:", column(Baryonic_freq_NSmass,1))

print("length Baryonic_freq_NSmass:", x_bins)
print("length Baryonic_freq_NSmass:", y_bins)

fig, ax = plt.subplots(figsize =(14, 7))
plt.hist2d(column(Baryonic_freq_NSmass,0),
                              column(Baryonic_freq_NSmass,1),
                              bins = [x_bins, y_bins],
                              cmap = plt.get_cmap('YlOrBr'),
                              cmin=1,
                              cmax=20)

#plt.xlabel('NS barycentric rotation frequency, $F_{0}$ [Hz]')
#plt.ylabel('#Entries')
#plt.title('Galactic coordinates of NS')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.xlim(0, 500)
#plt.ylim(0.5, 3)
#plt.grid(False)

#figsize_inch = 8, 6
#fig = plt.figure(figsize=figsize_inch, dpi=200)
#plt.rcParams["figure.figsize"] = (10,3)
#fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)
#fig.set_size_inches(14, 7)

# Adding color bar
plt.colorbar()
ax.set_xlabel('NS barycentric rotation frequency, $F_{0}$ [Hz]')
ax.set_ylabel('NS mass [$M_{sun}$]')
plt.savefig('plots/2d_NS_baryonic_freq_dist_vs_mass.png', dpi=200, transparent=True)
plt.savefig('plots/2d_NS_baryonic_freq_dist_vs_mass.pdf', dpi=200, transparent=True)
plt.savefig('plots/2d_NS_baryonic_freq_dist_vs_mass.eps', dpi=200, transparent=True)
plt.clf()


