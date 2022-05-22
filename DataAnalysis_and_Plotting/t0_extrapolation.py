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
from sklearn.linear_model import LinearRegression
from uncertainties import unumpy
from uncertainties import ufloat
from resizeFig import set_size
plt.figure(figsize=set_size(420))
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
#Take files from this directory
directory = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\beta6_000000\\24X24X24X24\\GF\\"

def calcE(Earray):
    df = Earray
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 10000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    return mean, std_error


E_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\"
directory_list = [E_dir+ "beta6_000000\\16X16X16X16\\GF\\",
                  E_dir+ "beta6_130000\\20X20X20X20\\GF\\",
                  E_dir+ "beta6_260000\\24X24X24X24\\GF\\",
                  E_dir+ "beta6_460000\\32X32X32X32\\GF\\"]
fileStart_list = [0,0,0,0]
fileEnd_list = [1703,1323,1203,644]#CHANGE
beta = [6.0,6.13,6.26,6.46]
V_list = [16**4,20**4,24**4,32**4]
a = [calc_a(entry) for entry in beta]
aOverR0_list = np.array([calc_aOverr0W_errors(entry) for entry in beta])
color_list = [blue,orange,purple,red]
marker_list = ['^','o','s','p']

#t_comp_lim_low_list = [2.7,4.0,6.5,11.5]
#t_comp_lim_high_list = [3.5,5.5,8.5,13]
t_comp_lim_low_list = [3.2,4.85,7.1,12.5]
t_comp_lim_high_list = [3.25,4.9,7.2,12.6]
r_0 = 0.5

t0_estimate_list = []
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
    t_comp_lim_low = t_comp_lim_low_list[j]
    t_comp_lim_high =t_comp_lim_high_list[j]
    PlotMarker = marker_list[j]
    #Turn the files into a single dataframe
    frames = []
    for i in range(fileStart,fileEnd+1,1):
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory+file,index_col='t',names = ['t','E'+str(i)], sep=",", header=None)
        frames.append(df1.loc[t_comp_lim_low:t_comp_lim_high])
    df2 = pd.concat(frames,axis=1)

    r_0 = 0.5
    tOverR_0= []
    tsquaredE = []
    errorList = []
    E_list = []
    for t_comp in df2.index.tolist():
        Earray = df2.loc[t_comp]
        Eval, error = calcE(Earray)
        #tOverR_0.append(t_comp*a_val**2/(r_0**2))
        tOverR_0.append(t_comp*aOverR0**2)
        tsquaredE.append(t_comp**2*Eval-0.3)
        errorList.append(t_comp**2*error)
        #E = ufloat(Eval,error)
        #E_list.append(0.3/E)

    plt.errorbar(tsquaredE,[item.nominal_value for item in tOverR_0],xerr=errorList,yerr=[item.std_dev for item in tOverR_0],markersize = 3.0,
                 fmt='o',ecolor = plot_color,color = plot_color,marker = PlotMarker,capsize=2,elinewidth=1,
                markeredgewidth=1)


    #x_fit = unumpy.uarray(np.c_[np.array(tOverR_0),np.ones_like(tOverR_0)],np.c_[np.zeros_like(tOverR_0),np.zeros_like(tOverR_0)])
    #y_fit = unumpy.uarray(np.array(tsquaredE),np.zeros_like(errorList))
    x_fit = unumpy.uarray(np.c_[np.array(tsquaredE),np.ones_like(tOverR_0)],np.c_[np.array(errorList),np.zeros_like(tOverR_0)])
    y_fit = np.array(tOverR_0)

    inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

    fit_a, fit_b= inv_mat.dot(x_fit.T.dot(y_fit))
    print("mean:")
    print(fit_b.nominal_value)
    print("std:")
    print(fit_b.std_dev)
    t0_estimate_list.append(fit_b)
    a_axis = np.linspace(-0.01,0.01,100)
    plt.plot(a_axis,fit_a.nominal_value*a_axis+fit_b.nominal_value,color = plot_color)
    
    label_text = "{}".format(fit_b)
    print(label_text)
    for char in "+/": label_text=label_text.replace(char,"")
    label_text=label_text.replace("-"," ")
    temp = label_text.split(" ")
    for i in range(len(temp[1])): 
        char = temp[1][i]
        if(char== '0' or char == '.'):
            continue
        temp[1] = temp[1][i:]
        break
    label_text = temp[0] +"(" +temp[1]+")"
    print(label_text)
    plt.errorbar(0.0,fit_b.nominal_value,yerr=fit_b.std_dev,markersize = 4.0,
                 fmt='o',ecolor = plot_color,color = plot_color,marker = PlotMarker,capsize=3,elinewidth=1,
                markeredgewidth=1.5,label = r'$t_0/r_0^2 ={}$, $a ={:.3f}$'.format(label_text,a_val))

plt.xlim(-0.004,0.004)
plt.ylim(0.103,0.118)
plt.ylabel(r'$t/r_0^2$')
plt.xlabel(r'$t^2\langle E \rangle-0.3$')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.legend()

#plt.savefig('t0_extrapolate_with_errors.pdf', bbox_inches="tight")
plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\FINALPLOTS\\t0_extrapolate_with_errors.pdf', bbox_inches="tight")
plt.clf()

t0_estimate_list = np.array(t0_estimate_list)
temp = np.array(aOverR0_list)**2/(t0_estimate_list)
means = np.array([item.nominal_value for item in temp])
std_devs = np.array([item.std_dev for item in temp])
print(t0_estimate_list)
smear_rads = unumpy.sqrt(8*t0_estimate_list)
x_fit = unumpy.uarray(np.c_[means,np.ones_like(a)],np.c_[std_devs,np.zeros_like(a)])
#y_fit = unumpy.uarray(np.array(t0_estimate_list),np.zeros_like(tOverR_0))
y_fit = smear_rads

inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))
fit_a, fit_b= inv_mat.dot(x_fit.T.dot(y_fit))
#print(fit_a)

a_axis = np.linspace(0.0,0.35,100)
plt.plot(a_axis,fit_a.nominal_value*a_axis+fit_b.nominal_value,color=plot_color)
print(std_devs)
plt.errorbar(means,[item.nominal_value for item in smear_rads],xerr=std_devs,yerr = [item.std_dev for item in smear_rads],markersize = 2.0,
                fmt='o',ecolor = purple,color = purple,capsize=2,elinewidth=1,
            markeredgewidth=1,marker=marker_list[1])

plt.errorbar(0.0,fit_b.nominal_value,yerr=fit_b.std_dev,markersize = 2.0,
                fmt='o',ecolor = plot_color,color = plot_color,capsize=2,elinewidth=1,
            markeredgewidth=1,marker=marker_list[0])
print(fit_b)
print(fit_b.nominal_value)
print(fit_b.std_dev)
plt.ylim(0.915,0.975)
plt.xlim(-0.012,0.35)
plt.ylabel(r'$\sqrt{8t_0}/r_0$')
plt.xlabel(r'$a^2/t_0$')
#plt.savefig('smearRad_cont_limit.pdf', bbox_inches="tight")
#plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\FINALPLOTS\\smearRad_cont_limit.pdf', bbox_inches="tight")
#C:\\Users\\adria\\Documents\\msc_project\\doc\\FINALPLOT\\