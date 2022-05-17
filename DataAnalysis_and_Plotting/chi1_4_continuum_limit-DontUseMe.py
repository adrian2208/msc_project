import os.path
import module1 as mod
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

def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))

#DO NOT COPY!!!!!!!!!!!!!!!!!!!!!!!!!! \/ \/ \/ \/ \/ \/
def TopSusc_for_regplot(topChargeArray,latticeConstant,latticeSize,t0):#DO NOT COPY THIS FUNCTION - ITS OUTPUT IS ONLY MEANT TO WORK FOR SNS.REGPLOT
    df = topChargeArray**2
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 10000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    out = t0**2*ufloat(mean,std_error)/ufloat(latticeConstant**4 * latticeSize,0)
    #mean = (mean/(latticeConstant**4*latticeSize))**(1/4)*197.3
    #return mean, propagated_error
    return out
#DO NOT COPY!!!!!!!!!!!!!!!!!!!!!!!!!! /\ /\ /\ /\ /\ /\

r_0 = 0.5

#RIGHT CLICK FILENAME AND RUN WITHOUT DEBUGGING (VS)

topC_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\"
directory_list = [topC_dir+ "beta6_000000\\16X16X16X16\\GF\\",
                  topC_dir+ "beta6_130000\\20X20X20X20\\GF\\",
                  topC_dir+ "beta6_260000\\24X24X24X24\\GF\\",
                  topC_dir+ "beta6_460000\\32X32X32X32\\GF\\For_susceptibility\\tau_int_leq0.5-2400updates\\"]

fileStart_list = [0,1,0,0]
fileEnd_list = [1000,1000,1226,103]#CHANGE
beta = [6.0,6.13,6.26,6.46]
V = [16**4,20**4,24**4,32**4]
extrapolate_at_t_comp = [3.2,4.75,7.3,12.6]
#t0_list = [0.027725,0.027575,0.0278,0.027875]
t0_list = np.array([0.1108,0.1108,0.1109,0.1115])
t0_list *=0.5**2

a = [calc_a(entry) for entry in beta]
aOverR0_list = np.array([calc_aOverr0W_errors(entry) for entry in beta])
mean_list = []
std_list = []
for j in range(len(directory_list)):
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    directory = directory_list[j]
    t_comp = extrapolate_at_t_comp[j]
    t0 = t0_list[j]
    frames = []
    for i in range(fileStart,fileEnd+1,1):
        try:
            file = 'torus_extdof4_'+str(i)+'.csv'
            df1 = pd.read_csv(directory+file,index_col='t',names = ['t','Q'], sep=",", header=None)
            frames.append(df1)
        except:
            print("Warning: file with name " + file + " is absent. Skipping it...")
            continue

    df2 = pd.concat(frames,axis=1).loc[t_comp]

    temp = TopSusc_for_regplot(df2,aOverR0_list[j]*0.5**2,V[j],t0)
    mean = temp.nominal_value
    std = temp.std_dev
    mean_list.append(mean)
    std_list.append(std)


#coeffs,cov = np.polyfit(np.array(a)**2,mean_list,1,w=1/np.array(std_list)**2,cov=True)
#print(cov)
#x = np.linspace(0,0.02,100)
#plt.plot(x,coeffs[0]*x+coeffs[1])
#plt.savefig('chi_fourthRoot_linear_extrapolation.pdf')
#print(x0)
#print(x1)
x_fit = unumpy.uarray(np.c_[np.array(a)**2/np.array(t0_list),np.ones_like(a)],np.c_[np.zeros_like(a),np.zeros_like(a)])
y_fit = unumpy.uarray(mean_list,std_list)
inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

print(x_fit)
print(y_fit)
print(inv_mat)

fit_a, fit_b = inv_mat.dot(x_fit.T.dot(y_fit))


print(fit_a.nominal_value)
print(fit_a.std_dev)
print('fit_a={}, fit_b={}'.format(fit_a, fit_b))
a_axis = np.linspace(0,0.01,100)

#plt.plot(a_axis,fit_a.nominal_value*a_axis+fit_b.nominal_value)
#plt.plot(a_axis,(fit_a.nominal_value-fit_a.std_dev)*a_axis+(fit_b.nominal_value+fit_b.std_dev))
#plt.plot(a_axis,(fit_a.nominal_value+fit_a.std_dev)*a_axis+(fit_b.nominal_value-fit_b.std_dev))

#plt.xlim(-0.0005,0.01)
#plt.ylim(170,200)
plt.errorbar(np.array(a)**2/np.array(t0_list),mean_list,std_list)
#plt.errorbar(0.0,fit_b.nominal_value,fit_b.std_dev)
plt.xlabel(r'$t/r_0^2$')
plt.ylabel(r'$t^2\langle E \rangle$')
plt.legend()
plt.savefig('chi_fourthRoot_linear_extrapolation.pdf')