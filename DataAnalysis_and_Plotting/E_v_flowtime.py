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
from resizeFig import set_size
plt.figure(figsize=set_size(420))
blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'
cyan = '#308FAC'
purple = '#7047AD'

def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))
def calc_aOverr0W_errors(beta):
    temp = np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))
    out = ufloat(temp,temp*(0.0034482758621*beta - 0.0166551724138))
    return out

def calcE(Earray):
    df = Earray
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 10000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    return mean, std_error

root = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\"
directory_list = [
        root + "beta6_000000\\16X16X16X16\\GF\\",
        root + "beta6_130000\\20X20X20X20\\GF\\", 
        root + "beta6_260000\\24X24X24X24\\GF\\",
        root + "beta6_460000\\32X32X32X32\\GF\\"]#change
fileStart_list = [0,0,0,0] #change
fileEnd_list  = [1703,1323,1203,644]#change
beta  = [6.0,6.13,6.26,6.46] #change
V_list  = [16**4,20**4,24**4,32**4] #change
color_list = [blue,orange,green,purple]
marker_list = ['^','o','s','p']
a = [calc_a(entry) for entry in beta]
aOverR0_list = np.array([calc_aOverr0W_errors(entry) for entry in beta])


mean_list = []
std_list = []
for j in range(len(directory_list)):
    directory = directory_list[j]
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    a_val = a[j]
    aOverR0 = aOverR0_list[j]
    V = V_list[j]
    plot_color = color_list[j]
    Plot_marker = marker_list[j]
    #Turn the files into a single dataframe
    frames = []
    for i in range(fileStart,fileEnd+1,1):
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory+file,index_col='t',names = ['t','E'+str(i)], sep=",", header=None)
        frames.append(df1)
    df2 = pd.concat(frames,axis=1)

    r_0 = 0.5
    tOverR_0= []
    tsquaredE = []
    errorList = []
    for t_comp in df2.index.tolist():
        Earray = df2.loc[t_comp]
        Eval, error = calcE(Earray)
        tOverR_0.append(t_comp*aOverR0**2)
        tsquaredE.append(Eval)
        errorList.append(error)

    plt.errorbar([item.nominal_value for item in tOverR_0],tsquaredE,errorList,xerr=[item.std_dev for item in tOverR_0],markersize = 2.0,
                 fmt='o',ecolor = plot_color,color = plot_color,marker = Plot_marker, capsize=2,elinewidth=1,
                markeredgewidth=1,label = r'$a =$'+ "{:.3f}".format(a_val))


plt.xlim(0)
plt.xlabel(r'$t/r_0^2$')
plt.ylabel(r'$\langle E \rangle$')
plt.legend()
plt.savefig('E_v_flowtime.pdf', bbox_inches="tight")