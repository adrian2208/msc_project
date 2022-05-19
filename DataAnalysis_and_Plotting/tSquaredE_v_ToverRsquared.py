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

def runningCoupling_PT(t,Lambda):
    PI = 3.14159265358979
    OneOverfourPi = 1/(4*PI)
    beta_0 = OneOverfourPi*11
    beta_1 = OneOverfourPi**2*102
    beta_2 = OneOverfourPi**3*(2857/2)
    beta_3 = OneOverfourPi**4*(149753/6+3564*1.202056903)
    r_0 = 0.5
    q = 1/(np.sqrt(8*t))
    t_Q = np.log(q**2/(Lambda**2))
    term_1 = beta_1/(beta_0**2)*np.log(t_Q)/t_Q
    term_2 =(beta_1**2*((np.log(t_Q)**2)-np.log(t_Q)-1)+beta_0*beta_2)/(beta_0**4*t_Q**2)
    term_3 = (beta_1**3*((np.log(t_Q))**3-5/2*(np.log(t_Q))**2-2*np.log(t_Q)+0.5)+
    3*beta_0*beta_1*beta_2*np.log(t_Q)-0.5*beta_0**2*beta_3)/(beta_0**6*t_Q**3)
    return 1/(beta_0*t_Q)*(1-term_1+term_2-term_3)
def tSquaredE_PT(t,Lambda):
    return 3/(4*3.14159265358979)*runningCoupling_PT(t,Lambda)*(1+1.0978*runningCoupling_PT(t,Lambda))
def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))
def calc_aOverr0W_errors(beta):
    temp = np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))
    out = ufloat(temp,temp*(0.0034482758621*beta - 0.0166551724138))
    return out
#Take files from this directory
directory = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\beta6_000000\\24X24X24X24\\GF\\"

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
        tsquaredE.append(t_comp**2*Eval)
        errorList.append(t_comp**2*error)

    plt.errorbar([item.nominal_value for item in tOverR_0],tsquaredE,errorList,xerr=[item.std_dev for item in tOverR_0],markersize = 2.0,
                 fmt='o',ecolor = plot_color,color = plot_color,marker = Plot_marker, capsize=2,elinewidth=1,
                markeredgewidth=1,label = r'$a =$'+ "{:.3f}".format(a_val))





tPT = np.linspace(0.00000001,0.05,4000)
plt.plot(tPT,tSquaredE_PT(tPT,0.602),color = 'grey')
plt.fill_between(tPT,tSquaredE_PT(tPT,0.602+0.048),tSquaredE_PT(tPT,0.602-0.048),facecolor='#30988A',edgecolor='#30988A',alpha=0.3, label = "Perturbation theory")
plt.xlim(0)
plt.xlabel(r'$t/r_0^2$')
plt.ylabel(r'$t^2\langle E \rangle$')
plt.legend()


plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\t_squared_v_tOverR_0_squared.pdf', bbox_inches="tight")