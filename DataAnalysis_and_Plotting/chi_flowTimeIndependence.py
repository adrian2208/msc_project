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
import random
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
    out = ufloat(mean,std_error)
    #mean = (mean/(latticeConstant**4*latticeSize))**(1/4)*197.3
    #return mean, propagated_error
    return out
#DO NOT COPY!!!!!!!!!!!!!!!!!!!!!!!!!! /\ /\ /\ /\ /\ /\
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
extrapolate_at_t_comp_list = [#3/4 : [2.4,3.6,5.4,9.5],
                              #1/2 : [1.6,2.4,3.6,6.3],
                              [0.8,1.2,1.8,3.0],#1/4
                              [0.55,0.8,1.2,2.0],#1/6
                              [0.25,0.4,0.6,1.0]]#1/12
t0_list = [3.2,4.85,7.1,12.2]
color_list = [blue,green,orange,purple,red]
marker_list = ['^','o','s','p']
label_list = [r'$t = 1/4t_0$',r'$t = 1/6t_0$',r'$t = 1/12t_0$']
a = [calc_a(entry) for entry in beta]
aOverR0 = [calc_aOverr0W_errors(entry) for entry in beta]
t0_list_phys = np.array([ufloat(0.1108,0.0009),ufloat(0.1107,0.0010),ufloat(0.1109,0.0011),ufloat(0.1115,0.0013)])
t0_list_phys *=0.5**2

at_t0_list = []
at_t0_Qarraylist = []
for j in range(len(directory_list)):
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    directory = directory_list[j]
    t_comp = t0_list[j]
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

    at_t0_list.append(TopSusc_for_regplot(df2,a[j],V[j]))
    at_t0_Qarraylist.append(df2)
    #mean_list_t0.append(mean)
    #std_list_t0.append(std)

rat_array = []
for k in range(len(extrapolate_at_t_comp_list)):
    extrapolate_at_t_comp = extrapolate_at_t_comp_list[k]
    rat_temp = []
    for j in range(len(directory_list)):
        fileStart = fileStart_list[j]
        fileEnd = fileEnd_list[j]
        directory = directory_list[j]
        t_comp = extrapolate_at_t_comp[j]
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

        #temp_num = TopSusc_for_regplot(df2,a[j],V[j])
        temp_den = at_t0_list[j]
        temp_denArray = at_t0_Qarraylist[j]


        temp_mean,temp_std_dev = RatioCovariance(df2**2,temp_denArray.to_numpy()**2,10000,fileEnd_list[j]-20)
        temp_std_dev = np.sqrt(temp_std_dev)

        #temp_mean = temp_num.nominal_value/temp_den.nominal_value
        #temp_std_dev = abs(temp_mean)*np.sqrt((temp_num.std_dev/temp_num.nominal_value)**2+(temp_den.std_dev/temp_den.nominal_value)**2)#- 2*0.0*(temp_num.std_dev*temp_den.std_dev/(temp_den.nominal_value*temp_num.nominal_value))
        rat_temp.append(ufloat(temp_mean,temp_std_dev))
        #print(directory)
        #print(temp_num)
        #print(temp_den)
        #print("---------")
    rat_array.append(rat_temp)

print(rat_array)
label_idx = 0
for elems,plot_color,plotMarker in zip(rat_array,color_list,marker_list):
    mean_list = [p.nominal_value for p in elems]
    std_list = [p.std_dev for p in elems]
    #plt.errorbar(np.array(a)**2/(0.5**2),mean_list,std_list,markersize = 4.0,
    #             fmt='o',ecolor = plot_color,color = plot_color,capsize=4,elinewidth=2,
    #            markeredgewidth=1)
    temp = (np.array(aOverR0)*0.5)**2/np.array(t0_list_phys)
    plt.errorbar([item.nominal_value for item in temp],mean_list,std_list,xerr=[item.std_dev for item in temp],markersize = 4.0,
                 fmt='o',ecolor = plot_color,color = plot_color,marker = plotMarker,capsize=4,elinewidth=2,
                markeredgewidth=1)


    #x_fit = unumpy.uarray(np.c_[np.array(a)**8/(0.5**8),np.array(a)**6/(0.5**6),np.array(a)**4/(0.5**4),np.array(a)**2/(0.5**2),np.ones_like(a)],np.c_[np.zeros_like(a),np.zeros_like(a),np.zeros_like(a),np.zeros_like(a),np.zeros_like(a)])
    #x_fit = unumpy.uarray(np.c_[np.array(a)**8/(t0_list_phys**4),np.array(a)**6/(t0_list_phys**3),np.array(a)**4/(t0_list_phys**2),np.array(a)**2/(t0_list_phys),np.ones_like(a)],np.c_[np.zeros_like(a),np.zeros_like(a),np.zeros_like(a),np.zeros_like(a),np.zeros_like(a)])
    #y_fit = unumpy.uarray(mean_list,std_list)
    #inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

    #fit_k, fit_p,fit_a, fit_b, fit_c = inv_mat.dot(x_fit.T.dot(y_fit))
    ##a_axis = np.linspace(0,0.035,100)
    #a_axis = np.linspace(0,0.35,100)
    #plt.plot(a_axis,fit_k.nominal_value*a_axis**4+fit_p.nominal_value*a_axis**3+fit_a.nominal_value*a_axis**2+fit_b.nominal_value*a_axis +fit_c.nominal_value,color = plot_color)
    
    
    #Restore \/
    #x_fit = unumpy.uarray(np.c_[np.array(a)**4/(0.5**4),np.array(a)**2/(0.5**2),np.ones_like(a)],np.c_[np.zeros_like(a),np.zeros_like(a),np.zeros_like(a)])
    x_fit = np.array(np.c_[np.array(a)**4/(t0_list_phys**2),np.array(a)**2/(t0_list_phys),np.ones_like(a)])
    y_fit = unumpy.uarray(mean_list,std_list)
    inv_mat = unumpy.ulinalg.pinv(x_fit.T.dot(x_fit))

    fit_a, fit_b, fit_c = inv_mat.dot(x_fit.T.dot(y_fit))
    print("fit_c:")
    print(fit_c)

    #a_axis = np.linspace(0,0.035,100)
    a_axis = np.linspace(0,0.35,100)
    plt.plot(a_axis,fit_a.nominal_value*a_axis**2+fit_b.nominal_value*a_axis +fit_c.nominal_value,color = plot_color)




    plt.errorbar(0.0,fit_c.nominal_value,fit_c.std_dev,markersize = 4.0,
                 fmt='o',ecolor = plot_color,color = plot_color,marker = plotMarker,capsize=4,elinewidth=2,
                markeredgewidth=1,label = label_list[label_idx])
    label_idx+=1
#print(fit_a.nominal_value)
#print(fit_a.std_dev)
#print('fit_a={}, fit_b={}'.format(fit_a, fit_b))


#plt.plot(a_axis,fit_a.nominal_value*a_axis+fit_b.nominal_value)
#plt.plot(a_axis,(fit_a.nominal_value-fit_a.std_dev)*a_axis+(fit_b.nominal_value+fit_b.std_dev))
#plt.plot(a_axis,(fit_a.nominal_value+fit_a.std_dev)*a_axis+(fit_b.nominal_value-fit_b.std_dev))

#plt.xlim(-0.001,0.04)
plt.ylim(0.4,1.5)
#plt.errorbar(np.array(a)**2/(1.0**2),mean_list,std_list)
#plt.errorbar(0.0,fit_b.nominal_value,fit_b.std_dev)
plt.ylabel(r'$\chi_t(t)/\chi_t(t_0)$')
#plt.xlabel(r'$a^2/r_0^2$')
plt.xlabel(r'$a^2/t_0$')
plt.legend()
#plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\chi_flowTime_independence.pdf')
#plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\chi_flowTime_independence_a^4_extrapolation.pdf')
plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\chi_flowTime_independence.pdf', bbox_inches="tight")