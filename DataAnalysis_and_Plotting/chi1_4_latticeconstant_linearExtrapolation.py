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
purple = '#7047AD'

def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))
def calc_aOverr0W_errors(beta):
    temp = np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))
    out = ufloat(temp,temp*(0.0034482758621*beta - 0.0166551724138))
    return out
#DO NOT COPY!!!!!!!!!!!!!!!!!!!!!!!!!! \/ \/ \/ \/ \/ \/
def TopSusc_for_regplot(topChargeArray,latticeConstant,latticeSize):#DO NOT COPY THIS FUNCTION - ITS OUTPUT IS ONLY MEANT TO WORK FOR SNS.REGPLOT
    df = topChargeArray**2
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 10000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    #propagated_error = 197.327*0.25*mean**(-3.0/4.0)/((latticeConstant**4 * latticeSize)**(0.25))*std_error
    #mean = (mean/(latticeConstant**4*latticeSize))**(1/4)*197.3
    temp = ufloat(mean,std_error)
    temp = (temp/(latticeConstant**4*latticeSize))**(1/4)*197.327
    mean = temp.nominal_value
    propagated_error = temp.std_dev
    return mean, propagated_error
#DO NOT COPY!!!!!!!!!!!!!!!!!!!!!!!!!! /\ /\ /\ /\ /\ /\

r_0 = 0.5

#RIGHT CLICK FILENAME AND RUN WITHOUT DEBUGGING (VS)

#Take files from this directory
topC_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\"
directory_list = [topC_dir+ "beta6_000000\\16X16X16X16\\GF\\",
                  topC_dir+ "beta6_130000\\20X20X20X20\\GF\\",
                  topC_dir+ "beta6_260000\\24X24X24X24\\GF\\",
                  topC_dir+ "beta6_460000\\32X32X32X32\\GF\\Quncorrelated-4050Updates\\"]
fileStart_list = [0,0,0,0]
fileEnd_list = [1703,1323,613,214]#CHANGE
beta = [6.0,6.13,6.26,6.46]
V = [16**4,20**4,24**4,32**4]
extrapolate_at_t_comp = [3.2,4.9,7.2,12.2]
t_extrList = [[3.1,3.15,3.2,3.25,3.3],[4.7,4.75,4.8,4.85,4.9,4.95,5],[6.8,6.9,7,7.1,7.2,7.3,7.4],[11.9,12,12.1,12.2,12.3,12.4,12.5]]
t0OverR_squared = np.array([ufloat(0.1115,0.0009),ufloat(0.1112,0.0010),ufloat(0.1110,0.0011),ufloat(0.1117,0.0013)])
#directory_list = ["C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta6_000000\\16X16X16X16\\GF\\"]
#fileStart_list = [0]
#fileEnd_list = [500]#CHANGE
#beta = [6.0]
#V = [16**4]
#extrapolate_at_t_comp = [3.5]

a = [calc_a(entry) for entry in beta]
aOverR0_list = np.array([calc_aOverr0W_errors(entry) for entry in beta])
mean_list = []
std_list = []
for j in range(len(directory_list)):
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    directory = directory_list[j]
    t_comp = extrapolate_at_t_comp[j]
    frames = []
    extra = (j==1 or j==2)
    for i in range(fileStart,fileEnd+1,1+extra*2):
        try:
            file = 'torus_extdof4_'+str(i)+'.csv'
            df1 = pd.read_csv(directory+file,index_col='t',names = ['t','Q'], sep=",", header=None)
            frames.append(df1)
        except:
            print("Warning: file with name " + file + " is absent. Skipping it...")
            continue


    #df2 = pd.concat(frames,axis=1)
    #extrapolList = []
    #for t in t_extrList[j]:
    #    temp = df2.loc[t]
    #    mean, std = TopSusc_for_regplot(temp,aOverR0_list[j]*0.5,V[j])
    #    extrapolList.append(ufloat(mean,std))
    #tOverR_squared = aOverR0_list[j]**2*np.array(t_extrList[j])
    #x_fit = unumpy.uarray(np.c_[np.array([item.nominal_value for item in tOverR_squared]),np.ones_like(tOverR_squared)],np.c_[np.array([item.std_dev for item in tOverR_squared]),np.zeros_like(tOverR_squared)])
    #y_fit = unumpy.uarray([item.nominal_value for item in extrapolList],[item.std_dev for item in extrapolList])
    #inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

    #fit_a, fit_b = inv_mat.dot(x_fit.T.dot(y_fit))
    #chi1over4AtT0_extrapolated = fit_a*t0OverR_squared[j]+fit_b
    #mean_list.append(chi1over4AtT0_extrapolated.nominal_value)
    #std_list.append(chi1over4AtT0_extrapolated.std_dev)


    df2 = pd.concat(frames,axis=1).loc[t_comp]

    mean, std = TopSusc_for_regplot(df2,aOverR0_list[j]*0.5,V[j])
    mean_list.append(mean)
    std_list.append(std)


#coeffs,cov = np.polyfit(np.array(a)**2,mean_list,1,w=1/np.array(std_list)**2,cov=True)
#print(cov)
#x = np.linspace(0,0.02,100)
#plt.plot(x,coeffs[0]*x+coeffs[1])
#plt.savefig('chi_fourthRoot_linear_extrapolation.pdf')
#print(x0)
#print(x1)
x_fit = unumpy.uarray(np.c_[np.array([item.nominal_value for item in aOverR0_list**2]),np.ones_like(aOverR0_list)],np.c_[np.array([item.std_dev for item in aOverR0_list**2]),np.zeros_like(aOverR0_list)])
y_fit = unumpy.uarray(mean_list,std_list)
inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

fit_a, fit_b = inv_mat.dot(x_fit.T.dot(y_fit))

print("Linear Fit:")
print(fit_b)
#print(fit_a.nominal_value)
#print(fit_a.std_dev)
#print('fit_a={}, fit_b={}'.format(fit_a, fit_b))
a_axis = np.linspace(0,0.035,100)

plt.plot(a_axis,fit_a.nominal_value*a_axis+fit_b.nominal_value,color=red)
#plt.plot(a_axis,(fit_a.nominal_value-fit_a.std_dev)*a_axis+(fit_b.nominal_value+fit_b.std_dev))
#plt.plot(a_axis,(fit_a.nominal_value+fit_a.std_dev)*a_axis+(fit_b.nominal_value-fit_b.std_dev))

plt.xlim(-0.001,0.035)
#plt.ylim(150,230)
plt.errorbar([item.nominal_value for item in aOverR0_list**2],mean_list,std_list,xerr=[item.std_dev for item in aOverR0_list**2],markersize = 2.0,
                fmt='o',ecolor = purple,color = purple,capsize=2,elinewidth=1,
            markeredgewidth=1)
plt.errorbar(0.0,fit_b.nominal_value,fit_b.std_dev,markersize = 2.0,
                fmt='o',ecolor = red,color = red,capsize=2,elinewidth=1,
            markeredgewidth=1)

x_fit = unumpy.uarray(np.c_[np.ones_like(aOverR0_list)],np.c_[np.zeros_like(aOverR0_list)])
y_fit = unumpy.uarray(mean_list,std_list)
inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

fit_b = inv_mat.dot(x_fit.T.dot(y_fit))
print("constant fit:")
print(fit_b[0])
print("data means:")
print(mean_list)
print("data std:")
print(std_list)
print("datapoint x-coor:")
[print(item) for item in aOverR0_list**2]
a_axis = np.linspace(-0.0005,a[0]**2/(0.5**2),100)
plt.plot(a_axis,np.full_like(a_axis,fit_b[0].nominal_value),color="black")
plt.errorbar(-0.0005,fit_b[0].nominal_value,fit_b[0].std_dev,markersize = 2.0,
                fmt='o',ecolor = "black",color = "black",capsize=2,elinewidth=1,
            markeredgewidth=1)
plt.xlabel(r'$a^2/r_0^2$')
plt.ylabel(r'$\chi_t^{1/4}[MeV]$')
#plt.legend()
plt.savefig('chi_fourthRoot_continuum_extrapolation.pdf', bbox_inches="tight")
#C:\\Users\\adria\\Documents\\msc_project\\doc\\FINALPLOTS\\

#--------------\/\/\/\/\/\/   FIT USING ONLY THREE FINEST LATTICE SPACINGS     \/\/\/\/\/--------------------------------------------------------------------------------------------------------------------------------------

#x_fit = unumpy.uarray(np.c_[np.array([item.nominal_value for item in aOverR0_list[1:]**2]),np.ones_like(aOverR0_list[1:])],np.c_[np.array([item.std_dev for item in aOverR0_list[1:]**2]),np.zeros_like(aOverR0_list[1:])])
#y_fit = unumpy.uarray(mean_list[1:],std_list[1:])
#inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

#fit_a, fit_b = inv_mat.dot(x_fit.T.dot(y_fit))

#print("Linear Fit:")
#print(fit_b)
##print(fit_a.nominal_value)
##print(fit_a.std_dev)
##print('fit_a={}, fit_b={}'.format(fit_a, fit_b))
#a_axis = np.linspace(0,0.035,100)

#plt.plot(a_axis,fit_a.nominal_value*a_axis+fit_b.nominal_value,color=red)
##plt.plot(a_axis,(fit_a.nominal_value-fit_a.std_dev)*a_axis+(fit_b.nominal_value+fit_b.std_dev))
##plt.plot(a_axis,(fit_a.nominal_value+fit_a.std_dev)*a_axis+(fit_b.nominal_value-fit_b.std_dev))

#plt.xlim(-0.001,0.035)
##plt.ylim(150,230)
#plt.errorbar([item.nominal_value for item in aOverR0_list**2],mean_list,std_list,xerr=[item.std_dev for item in aOverR0_list**2],markersize = 2.0,
#                fmt='o',ecolor = purple,color = purple,capsize=2,elinewidth=1,
#            markeredgewidth=1)
#plt.errorbar(0.0,fit_b.nominal_value,fit_b.std_dev,markersize = 2.0,
#                fmt='o',ecolor = red,color = red,capsize=2,elinewidth=1,
#            markeredgewidth=1)

#x_fit = unumpy.uarray(np.c_[np.ones_like(aOverR0_list[1:])],np.c_[np.zeros_like(aOverR0_list[1:])])
#y_fit = unumpy.uarray(mean_list[1:],std_list[1:])
#inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

#fit_b = inv_mat.dot(x_fit.T.dot(y_fit))
#print("constant fit:")
#print(fit_b[0])
#print("data means:")
#print(mean_list)
#print("data std:")
#print(std_list)
#print("datapoint x-coor:")
#[print(item) for item in aOverR0_list**2]
#a_axis = np.linspace(-0.0005,a[0]**2/(0.5**2),100)
#plt.plot(a_axis,np.full_like(a_axis,fit_b[0].nominal_value),color="black")
#plt.errorbar(-0.0005,fit_b[0].nominal_value,fit_b[0].std_dev,markersize = 2.0,
#                fmt='o',ecolor = "black",color = "black",capsize=2,elinewidth=1,
#            markeredgewidth=1)
#plt.xlabel(r'$a^2/r_0^2$')
#plt.ylabel(r'$\chi_t^{1/4}[MeV]$')
##plt.legend()
#plt.savefig('chi_fourthRoot_continuum_extrapolation.pdf', bbox_inches="tight")
##C:\\Users\\adria\\Documents\\msc_project\\doc\\FINALPLOTS\\