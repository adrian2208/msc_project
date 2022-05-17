import module1 as mod
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import numpy as np
from resizeFig import set_size
plt.figure(figsize=set_size(420))

blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'
purple = '#7047AD'

#data = [[data]]
def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))
r_0 = 0.5
topC_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\"
directory_list = [[topC_dir+ "beta6_000000\\16X16X16X16\\GF\\"],
                  [topC_dir+ "beta6_130000\\20X20X20X20\\GF\\"],
                  [topC_dir+ "beta6_260000\\24X24X24X24\\GF\\"],
                  [topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp0\\",
                   topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp1\\",
                   topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp2\\",
                   topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp3\\",
                   topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp4\\"]]

topE_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\"
directory_list_E = [[topE_dir+ "beta6_000000\\16X16X16X16\\GF\\"],
                  [topE_dir+ "beta6_130000\\20X20X20X20\\GF\\"],
                  [topE_dir+ "beta6_260000\\24X24X24X24\\GF\\"],
                  [topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp0\\",
                   topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp1\\",
                   topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp2\\",
                   topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp3\\",
                   topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp4\\"]]

fileStart_list = [0,0,0,0]#50
fileEnd_list = [1703,1323,613,42]#CHANGE

fileStart_list_E = [0,0,0,0]#50
fileEnd_list_E = [1703,1323,1203,128]#CHANGE

beta = [6.0,6.13,6.26,6.46]
autoCorr_at_t_comp = [3.5,5.5,8.5,13]

a = [calc_a(entry) for entry in beta]
#Turn the files into a single dataframe
t0_overR_0Squared = []
tau_int_Q = []
tau_int_Qsquared = []
error_Q = []
error_Qsquared = []
for j in range(len(fileStart_list)):
    directory = directory_list[j]
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    t_comp = autoCorr_at_t_comp[j]
    #plot_color = color_list[j]
    autocorrArray = []
    for k in range(len(directory)):
        frames = []
        dir = directory[k]
        for i in range(fileStart,fileEnd+1,1):#4
            file = 'torus_extdof4_'+str(i)+'.csv'
            df1 = pd.read_csv(dir+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
            frames.append(df1)
        df2 = pd.concat(frames,axis=1)
        autocorrArray.append(df2.loc[t_comp].to_numpy())

    t0_overR_0Squared.append(a[j]/(r_0))

    QArray = [autocorrArray]
    mean, var, tau_Estimate, error_temp = mod.tauint(QArray,0)
    tau_int_Q.append(tau_Estimate)
    error_Q.append(error_temp)

    QArray = [np.array(autocorrArray)**2]
    mean, var, tau_Estimate, error_temp = mod.tauint(QArray,0)
    tau_int_Qsquared.append(tau_Estimate)
    error_Qsquared.append(error_temp)

tau_int_E = []
error_E = []
for j in range(len(fileStart_list)):
    directory = directory_list_E[j]
    fileStart = fileStart_list_E[j]
    fileEnd = fileEnd_list_E[j]
    t_comp = autoCorr_at_t_comp[j]
    #plot_color = color_list[j]
    autocorrArray = []
    for k in range(len(directory)):
        frames = []
        dir = directory[k]
        for i in range(fileStart,fileEnd+1,2):#4
            file = 'torus_extdof4_'+str(i)+'.csv'
            df1 = pd.read_csv(dir+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
            frames.append(df1)
        df2 = pd.concat(frames,axis=1)
        autocorrArray.append(df2.loc[t_comp].to_numpy())

    EArray = [autocorrArray]
    mean, var, tau_Estimate, error_temp = mod.tauint(EArray,0)
    tau_int_E.append(tau_Estimate)
    error_E.append(error_temp)


plt.errorbar(t0_overR_0Squared,tau_int_Q,error_Q,markersize = 2.0,
                fmt='o',ecolor = blue,color = blue,capsize=2,elinewidth=1,
            markeredgewidth=1,label = r'$Q$')

plt.errorbar(t0_overR_0Squared,tau_int_Qsquared,error_Qsquared,markersize = 2.0,
                fmt='o',ecolor = red,color = red,capsize=2,elinewidth=1,
            markeredgewidth=1,label = r'$Q^2$')

plt.errorbar(t0_overR_0Squared,tau_int_E,error_E,markersize = 2.0,
                fmt='o',ecolor = purple,color = purple,capsize=2,elinewidth=1,
            markeredgewidth=1,label = r'$E$')
#plt.xlim(0)
#plt.ylim(100,350)
plt.xlabel(r'$a/r_0$')
plt.ylabel(r'$\tau_{int}$')
plt.legend()

plt.savefig('Statistical_independece.pdf', bbox_inches="tight")