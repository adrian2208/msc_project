import os.path
import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
from scipy import stats
import numpy as np
from sklearn.utils import resample
from uncertainties import unumpy
from uncertainties import ufloat
blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'
purple = '#7047AD'

def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))

def calc_aOverr0W_errors(beta):
    temp = np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))
    out = ufloat(temp,temp*(0.0034482758621*beta - 0.0166551724138))
    return out

r_0 = 0.5

#RIGHT CLICK FILENAME AND RUN WITHOUT DEBUGGING (VS)

#Take files from this directory
topE_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\"
directory_list = [topE_dir+ "beta6_000000\\16X16X16X16\\GF\\",
                  topE_dir+ "beta6_130000\\20X20X20X20\\GF\\",
                  topE_dir+ "beta6_260000\\24X24X24X24\\GF\\",
                  topE_dir+ "beta6_460000\\32X32X32X32\\GF\\"]

fileStart_list = [0,1,0,0]
fileEnd_list = [1000,1000,148,103]#CHANGE
beta = [6.0,6.13,6.26,6.46]
V = [16**4,20**4,24**4,32**4]
extrapolate_at_t_comp = [3.2,4.75,7.3,12.6]


a = [calc_a(entry) for entry in beta]
aOverR0_list = np.array([calc_aOverr0W_errors(entry) for entry in beta])
t0_list = []
std_list = []



x_fit = unumpy.uarray(np.c_[np.array([item.nominal_value for item in aOverR0_list**2]),np.ones_like(aOverR0_list)],np.c_[np.array([item.std_dev for item in (aOverR0_list**2)]),np.zeros_like(aOverR0_list)])
y_fit = unumpy.uarray(mean_list,std_list)
inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

fit_a, fit_b = inv_mat.dot(x_fit.T.dot(y_fit))


print(fit_a.nominal_value)
print(fit_a.std_dev)
print('fit_a={}, fit_b={}'.format(fit_a, fit_b))
a_axis = np.linspace(0,0.035,100)

plt.plot(a_axis,fit_a.nominal_value*a_axis+fit_b.nominal_value,color=red)


plt.xlim(-0.0005,0.035)
#plt.ylim(175,200)
plt.errorbar([item.nominal_value for item in aOverR0_list**2],mean_list,std_list,xerr = [item.std_dev for item in aOverR0_list**2],markersize = 2.0,
                fmt='o',ecolor = purple,color = purple,capsize=2,elinewidth=1,
            markeredgewidth=1)
plt.errorbar(0.0,fit_b.nominal_value,fit_b.std_dev,markersize = 2.0,
                fmt='o',ecolor = red,color = red,capsize=2,elinewidth=1,
            markeredgewidth=1)

plt.xlabel(r'$a^2/r_0^2$')
plt.ylabel(r'$\sqrt{8t_0}/r_0$')
plt.legend()
plt.savefig('Smearing_radius_continuum_extrapolation.pdf')
