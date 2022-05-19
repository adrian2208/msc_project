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

def runningCoupling_PT(q,Lambda):
    PI = 3.14159265358979
    OneOverfourPi = 1/(4*PI)
    beta_0 = OneOverfourPi*11
    beta_1 = OneOverfourPi**2*102
    beta_2 = OneOverfourPi**3*(2857/2)
    beta_3 = OneOverfourPi**4*(149753/6+3564*1.202056903)
    r_0 = 0.5
    #q = 1/(np.sqrt(8*t))
    t_Q = np.log(q**2/(Lambda**2))
    term_1 = beta_1/(beta_0**2)*np.log(t_Q)/t_Q
    term_2 =(beta_1**2*((np.log(t_Q)**2)-np.log(t_Q)-1)+beta_0*beta_2)/(beta_0**4*t_Q**2)
    term_3 = (beta_1**3*((np.log(t_Q))**3-5/2*(np.log(t_Q))**2-2*np.log(t_Q)+0.5)+
    3*beta_0*beta_1*beta_2*np.log(t_Q)-0.5*beta_0**2*beta_3)/(beta_0**6*t_Q**3)
    return 1/(beta_0*t_Q)*(1-term_1+term_2-term_3)
def tSquaredE_PT(q,Lambda):
    return 3/(4*3.14159265358979)*runningCoupling_PT(q,Lambda)*(1+1.0978*runningCoupling_PT(q,Lambda))
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
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 1000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    return mean, std_error

root = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\"
#directory_list = [root + "beta6_000000\\16X16X16X16\\GF\\", 
#        root + "beta6_130000\\20X20X20X20\\GF\\", 
#        root + "beta6_260000\\24X24X24X24\\GF\\",
#        root + "beta6_460000\\32X32X32X32\\GF\\"]#change
#fileStart_list = [0,0,0,0] #change
#fileEnd_list  = [628,281,1203,701]#change
#beta  = [6.0,6.13,6.26,6.46] #change
#V_list  = [16**4,20**4,24**4,32**4] #change
#color_list = [blue,orange,green,purple]
directory_list = [
        root + "beta6_130000\\20X20X20X20\\GF\\", 
        root + "beta6_260000\\24X24X24X24\\GF\\",
        root + "beta6_460000\\32X32X32X32\\GF\\"]#change
fileStart_list = [0,0,0] #change
fileEnd_list  = [1323,1203,644]#change
beta  = [6.13,6.26,6.46] #change
V_list  = [20**4,24**4,32**4] #change
color_list = [orange,green,purple]
marker_list = ['^','o','s','p']

a = [calc_a(entry) for entry in beta]
aOverR0_list = np.array([calc_aOverr0W_errors(entry) for entry in beta])

Fit_tOverR0Squared_start =0.0035
Fit_tOverR0Squared_stop = 0.07#0.04#0.115
fit_nrPoints = 100
minSize_fit = 5

Lambda_min = 200
Lambda_max = 260
Nr_Lambdas = 400

contLimit_fitDomain = np.linspace(Fit_tOverR0Squared_start,Fit_tOverR0Squared_stop,fit_nrPoints)
#contLimit_fitDomain = np.linspace(2000,600,fit_nrPoints)
#contLimit_fitDomain = np.array([(1/item*197.327)**2/2 for item in contLimit_fitDomain])

mean_list = []
std_list = []
interpolated_array = []
for j in range(len(directory_list)):
    directory = directory_list[j]
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    a_val = a[j]
    aOverR0 = aOverR0_list[j]
    V = V_list[j]
    plot_color = color_list[j]
    plot_marker = marker_list[j]
    #Turn the files into a single dataframe
    frames = []
    for i in range(fileStart,fileEnd+1,1):
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory+file,index_col='t',names = ['t','E'+str(i)], sep=",", header=None)
        frames.append(df1)
    df2 = pd.concat(frames,axis=1)

    r_0 = 0.5
    smear_rad= []
    tsquaredE = []
    errorList = []
    for t_comp in df2.index.tolist()[1:]:
        Earray = df2.loc[t_comp]
        Eval, error = calcE(Earray)
        smear_rad.append(1/unumpy.sqrt(t_comp*aOverR0**2*8*0.5**2)*197.327)
        tsquaredE.append(t_comp**2*Eval)
        errorList.append(t_comp**2*error)
    #RESTORE\/\/\/\/\/\/\/\/\/
    plt.errorbar([item.nominal_value for item in smear_rad],tsquaredE,errorList,xerr=[item.std_dev for item in smear_rad],markersize = 2.0,
                    fmt='o',ecolor = plot_color,color = plot_color,marker = plot_marker, capsize=2,elinewidth=1,
                markeredgewidth=1,label = r'$a =$'+ "{:.3f}".format(a_val))
    
    list = df2.index.tolist()[1:]
    temp_interp_array = []
    for time in range(1,len(list)-1):
        tOverR_0 = [list[time]*aOverR0**2,list[time+1]*aOverR0**2]
        seq_tsquaredE = [tsquaredE[time],tsquaredE[time+1]]
        seq_errorList = [errorList[time],errorList[time+1]]
        x_fit = unumpy.uarray(np.c_[np.array([item.nominal_value for item in tOverR_0]),np.ones_like(tOverR_0)],np.c_[np.array([item.std_dev for item in tOverR_0]),np.zeros_like(tOverR_0)])
        y_fit = unumpy.uarray(seq_tsquaredE,seq_errorList)

        inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

        fit_a, fit_b= inv_mat.dot(x_fit.T.dot(y_fit))
        performFitFor = np.array([item for item in contLimit_fitDomain if( (item >= tOverR_0[0]) and (item <= tOverR_0[1]))])
        fitfunc = fit_a*performFitFor+fit_b
        temp_interp_array.extend(fitfunc)
    interpolated_array.append(temp_interp_array)
    #plt.plot(1/unumpy.sqrt(contLimit_fitDomain*8*0.5**2)*197.327,[item.nominal_value for item in temp_interp_array],color="red")
    #plt.errorbar(1/unumpy.sqrt(contLimit_fitDomain*8*0.5**2)*197.327,[item.nominal_value for item in temp_interp_array],[item.std_dev for item in temp_interp_array],markersize = 2.0,
    #            fmt='o',ecolor = "red",color = "red",capsize=2,elinewidth=1,
    #        markeredgewidth=1)
interpolated_array = np.array(interpolated_array)
cont_limit_array = []
for time, i in zip(contLimit_fitDomain,range(len(contLimit_fitDomain))):
    fit_temp = np.array(interpolated_array[:,i])
    x_fit = unumpy.uarray(np.c_[np.array([item.nominal_value for item in aOverR0_list**2]),np.ones_like(aOverR0_list)],np.c_[np.array([item.std_dev for item in aOverR0_list**2]),np.zeros_like(aOverR0_list)])
    y_fit = unumpy.uarray(np.array([item.nominal_value for item in fit_temp]),np.array([item.std_dev for item in fit_temp]))

    inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

    fit_a, fit_b= inv_mat.dot(x_fit.T.dot(y_fit))
    cont_limit_array.append(fit_b)

#[print(item) for item in interpolated_array]
#plt.xlim(0)
print(len(smear_rad))
print(len(smear_rad))
contLimit_fitDomain
fit_x_axis = 1/unumpy.sqrt(contLimit_fitDomain*8*0.5**2)*197.327

#RESTORE\/\/\/\/\/\/
plt.errorbar(fit_x_axis,[item.nominal_value for item in cont_limit_array],yerr=[item.std_dev for item in cont_limit_array],markersize = 2.0,
                fmt='o',ecolor = "black",color = "black",capsize=2,elinewidth=1,
            markeredgewidth=1,label = r'Continuum limit')

plt.xlabel(r'$1/\sqrt{8t}[MeV]$')
plt.ylabel(r'$t^2\langle E \rangle$')
plt.legend()
plt.ylim(0,0.25)
plt.xlim(300,2000)

plt.savefig('t_squared_v_q_contLimit.pdf', bbox_inches="tight")
#C:\\Users\\adria\\Documents\\msc_project\\doc\\
plt.clf()

cont_limit_datapoints_mean = np.array([item.nominal_value for item in cont_limit_array])
cont_limit_datapoints_std = np.array([item.std_dev for item in cont_limit_array])
OneOverChiSquared_list = []

LambdaList = np.linspace(Lambda_min,Lambda_max,Nr_Lambdas)
OneOverChiSquared_list=np.zeros(len(LambdaList))
NrFits = 0
for i in range(len(fit_x_axis)):
    for j in range(i, len(fit_x_axis)):
        if(abs(i-j)>minSize_fit):
            fit_interval = fit_x_axis[i:j+1]
            fit_mean = cont_limit_datapoints_mean[i:j+1]
            fit_std = cont_limit_datapoints_std[i:j+1]
            for Lambda,k in zip(LambdaList,range(len(LambdaList))):
                perturbation_theory_prediction = np.array([tSquaredE_PT(item,Lambda) for item in fit_interval])
                chi_squared = np.sum(((fit_mean-perturbation_theory_prediction)**2)/(fit_std**2))/len(perturbation_theory_prediction)
                OneOverChiSquared_list[k] += (1/chi_squared)
            NrFits += 1
    print("{} out of {}".format(i,len(fit_x_axis)))

norm = np.sum(OneOverChiSquared_list)
df = pd.DataFrame(data = OneOverChiSquared_list/norm, index = LambdaList)
df.to_csv("ChiSquaredFit.csv",header = False)
plt.plot(LambdaList,OneOverChiSquared_list/norm)
plt.xlabel(r'$\Lambda_{YM}[MeV]$')
#plt.ylabel(r'$1/\chi^2$')
#plt.legend()
#plt.ylim(0,0.25)
plt.xlim(Lambda_min,Lambda_max)

plt.savefig('t_squared_v_q_ChiSquared_fitDistribution.pdf', bbox_inches="tight")
plt.clf()
plt.hist(LambdaList,bins = int(len(LambdaList)/4),weights = OneOverChiSquared_list/norm,color = "#D5824B")
plt.xlabel(r'$\Lambda_{YM}[MeV]$')
#plt.ylabel(r'$1/\chi^2$')
#plt.legend()
#plt.ylim(0,0.25)
#plt.xlim(Lambda_min,Lambda_max)

plt.savefig('t_squared_v_q_ChiSquared_fitDistribution_HIST.pdf', bbox_inches="tight")
#C:\\Users\\adria\\Documents\\msc_project\\doc\\