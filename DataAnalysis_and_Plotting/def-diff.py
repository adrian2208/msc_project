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
from resizeFig import set_size
from uncertainties import unumpy
from uncertainties import ufloat
import random
plt.figure(figsize=set_size(420))
blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'
purple = '#7047AD'
def TopSusc(topChargeArray,latticeConstant,latticeSize):#DO NOT COPY THIS FUNCTION - ITS OUTPUT IS ONLY MEANT TO WORK FOR SNS.REGPLOT
    df = topChargeArray**2
    result= stats.bootstrap([df],np.mean,n_resamples = 10000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    out = ufloat(mean,std_error)
    #mean = (mean/(latticeConstant**4*latticeSize))**(1/4)*197.3
    #return mean, propagated_error
    return out
def RatioCovariance(topChargeArray1,topChargeArray2,n_resamples,sampleSize):
    means = []
    vars = []
    for i in range(n_resamples):
        sample1 = np.array(random.sample(list(topChargeArray1),sampleSize))
        sample2 = np.array(random.sample(list(topChargeArray2),sampleSize))
        #vars.append(np.var(sample1/sample2))
        means.append(np.mean(sample1)/np.mean(sample2))
        #print(sample1)
        #print(sample2)
    return np.mean(means), np.var(means)
def calc_aOverr0W_errors(beta):
    temp = np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))
    out = ufloat(temp,temp*(0.0034482758621*beta - 0.0166551724138))
    return out
def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))

root = "C:\\Users\\adria\\Documents\\def-diff\\Topological_Charge\\"
dirs = [root + "beta6_000000\\16X16X16X16\\GF\\", 
        root + "beta6_130000\\20X20X20X20\\GF\\", 
        root + "beta6_260000\\24X24X24X24\\GF\\",
        root + "beta6_460000\\32X32X32X32\\GF\\"]

rootUnimproved = "C:\\Users\\adria\\Documents\\def-diff\\UnImprovedTopological_Charge\\"
dirsUnimproved = [rootUnimproved + "beta6_000000\\16X16X16X16\\GF\\", 
        rootUnimproved + "beta6_130000\\20X20X20X20\\GF\\", 
        rootUnimproved + "beta6_260000\\24X24X24X24\\GF\\",
        rootUnimproved + "beta6_460000\\32X32X32X32\\GF\\"]
        
        
#root = "C:\\Users\\adria\\Documents\\def-diff\\ImprovedEnergyDensity\\"
#dirs = [root + "beta6_000000\\16X16X16X16\\GF\\", 
#        root + "beta6_130000\\20X20X20X20\\GF\\", 
#        root + "beta6_260000\\24X24X24X24\\GF\\",
#        root + "beta6_460000\\32X32X32X32\\GF\\"]

#rootUnimproved = "C:\\Users\\adria\\Documents\\def-diff\\EnergyDensity\\"
#dirsUnimproved = [rootUnimproved + "beta6_000000\\16X16X16X16\\GF\\", 
#        rootUnimproved + "beta6_130000\\20X20X20X20\\GF\\", 
#        rootUnimproved + "beta6_260000\\24X24X24X24\\GF\\",
#        rootUnimproved + "beta6_460000\\32X32X32X32\\GF\\"]

fileStart_list = [1001,1001,1001,0]
fileEnd_list = [1907,1907,1907,214]
beta_list  = [6.0,6.13,6.26,6.46]
V = [16**4,20**4,24**4,32**4]
aoverR0 = np.array([calc_aOverr0W_errors(entry) for entry in beta_list])
t_comp_list = [3.5,5.5,8.5,13.0]
color_list = [blue,orange,green,red]

#RIGHT CLICK FILENAME AND RUN WITHOUT DEBUGGING (VS)
std_list = []
meanList = []
for j in range(len(dirs)):
    #Take files from this directory
    directory = dirs[j]
    directory_unimproved = dirsUnimproved[j]
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    plot_color = color_list[j]
    t_comp = t_comp_list[j]
    #Turn the files into a single dataframe
    frames = []
    for i in range(fileStart,fileEnd+1,1):
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
        frames.append(df1)
    df2 = pd.concat(frames,axis=1).loc[t_comp].to_numpy()

    frames_unimproved = []
    for i in range(fileStart,fileEnd+1,1):
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory_unimproved+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
        frames_unimproved.append(df1)
    df2_unimproved = pd.concat(frames_unimproved,axis=1).iloc[-1].to_numpy()
    print(df2_unimproved)
    print(df2)
    #diff = np.mean(abs(df2-df2_unimproved))
    #meanList.append(TopSusc(df2_unimproved,aoverR0[j]*0.5,V[j]))
    #meanList.append(TopSusc(df2_unimproved,aoverR0[j]*0.5,V[j])/TopSusc(df2,aoverR0[j]*0.5,V[j]))
    mean,std = RatioCovariance(df2**2,df2_unimproved**2,15000,(fileEnd_list[j]-fileStart_list[j])-5)
    meanList.append(mean)
    std_list.append(np.sqrt(std))
    #difflist.append(diff/np.mean(df2))
print(meanList)



x_fit = np.array(np.c_[aoverR0**4,aoverR0**2,np.ones_like(aoverR0)])
y_fit = unumpy.uarray(meanList,std_list)
inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

fit_a, fit_b, fit_c = inv_mat.dot(x_fit.T.dot(y_fit))
print("fit_c:")
print(fit_c)
print("fit_b:")
print(fit_b)
print("fit_a:")
print(fit_a)

#a_axis = np.linspace(0,0.035,100)
a_axis = np.linspace(0,0.035,100)
plt.plot(a_axis,fit_a.nominal_value*a_axis**2+fit_b.nominal_value*a_axis +fit_c.nominal_value,color = plot_color)

plt.errorbar([item.nominal_value for item in aoverR0**2],meanList,std_list,xerr=[item.std_dev for item in aoverR0**2],markersize = 2.0,
                fmt='o',ecolor = purple,color = purple,capsize=2,elinewidth=1,
            markeredgewidth=1)


plt.errorbar(0.0,fit_c.nominal_value,fit_c.std_dev,markersize = 4.0,
                fmt='o',ecolor = red,color = red,capsize=4,elinewidth=2,
            markeredgewidth=1)
#plt.scatter(np.array(aoverR0),np.array(difflist))
#plt.xlim(-0.01,0.035)
plt.xlim(-0.0025,0.036)
plt.ylim(0.94)
plt.xlabel(r'$a^2/r_0^2$')
plt.ylabel(r'$\frac{\chi_t^{\frac{1}{4} } (Q_I)}{\chi_t^{\frac{1}{4} } (Q_C)}$',rotation = 0,labelpad=14)
plt.savefig('E_def_diff.pdf', bbox_inches="tight")


