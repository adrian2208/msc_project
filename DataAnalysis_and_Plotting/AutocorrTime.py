import pandas as pd
import puwr as mod
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
#topC_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\"
#directory_list = [[topC_dir+ "beta6_000000\\16X16X16X16\\GF\\"],
#                  [topC_dir+ "beta6_130000\\20X20X20X20\\GF\\"],
#                  [topC_dir+ "beta6_260000\\24X24X24X24\\GF\\"],
#                  [topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp0\\",
#                   topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp1\\",
#                   topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp2\\",
#                   topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp3\\",
#                   topC_dir+ "beta6_460000\\32X32X32X32\\GF\\4050Updates_forAutocorr\\temp4\\"]]

#topE_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\"
#directory_list = [[topE_dir+ "beta6_000000\\16X16X16X16\\GF\\"],
#                  [topE_dir+ "beta6_130000\\20X20X20X20\\GF\\"],
#                  [topE_dir+ "beta6_260000\\24X24X24X24\\GF\\"],
#                  [topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp0\\",
#                   topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp1\\",
#                   topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp2\\",
#                   topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp3\\",
#                   topE_dir+ "beta6_460000\\32X32X32X32\\GF\\1350Updates_forAutocorr\\temp4\\"]]

topC_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\"
directory_list = [[topC_dir+ "beta6_000000\\12X12X12X12\\GF\\"],
                  [topC_dir+ "beta6_130000\\14X14X14X14\\GF\\"],
                  [topC_dir+ "beta6_260000\\16X16X16X16\\GF\\autocorr0\\",
                   topC_dir+ "beta6_260000\\16X16X16X16\\GF\\autocorr1\\"],
                  [topC_dir+ "beta6_460000\\22X22X22X22\\GF\\autocorr0\\",
                   topC_dir+ "beta6_460000\\22X22X22X22\\GF\\autocorr1\\"]]
#topC_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\"
#directory_list = [[topC_dir+ "beta6_000000\\12X12X12X12\\GF\\"],
#                  [topC_dir+ "beta6_130000\\14X14X14X14\\GF\\"],
#                  [topC_dir+ "beta6_260000\\16X16X16X16\\GF\\autocorr0\\",
#                   topC_dir+ "beta6_260000\\16X16X16X16\\GF\\autocorr1\\"],
#                  [topC_dir+ "beta6_460000\\22X22X22X22\\GF\\autocorr0\\",
#                   topC_dir+ "beta6_460000\\22X22X22X22\\GF\\autocorr1\\"]]

fileStart_list = [0,0,0,0]#50
fileEnd_list = [2000,2000,2000,5000]#CHANGE
colorList = [green,orange,blue,purple]
marker_list = ['^','o','s','p']
#fileStart_list_E = [0,0,0,0]#50
#fileEnd_list_E = [1703,1323,1203,128]#CHANGE

beta = [6.0,6.13,6.26,6.46]
autoCorr_at_t_comp = [3.5,5.5,8.5,13]

a = [calc_a(entry) for entry in beta]
#Turn the files into a single dataframe

for j in range(len(fileStart_list)):
    t0_overR_0Squared = []
    tau_int_Q = []
    tau_int_Qsquared = []
    error_Q = []
    error_Qsquared = []
    directory = directory_list[j]
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    t_comp = autoCorr_at_t_comp[j]
    Plot_marker = marker_list[j]
    #plot_color = color_list[j]
    autocorrArray = []
    for k in range(len(directory)):
        frames = []
        dir = directory[k]
        for i in range(fileStart,fileEnd+1,1):
            file = 'torus_extdof4_'+str(i)+'.csv'
            df1 = pd.read_csv(dir+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
            frames.append(df1)
        df2 = pd.concat(frames,axis=1)
        autocorrArray.append(df2)

    t0_overR_0Squared.append(a[j]/(r_0))
    t_comp_list = []

    for t_comp in df2.index.tolist():
        QArray = [[item.loc[t_comp].to_numpy() for item in autocorrArray]]
        mean, var, tau_Estimate, error_temp = mod.tauint(QArray,0)
        
        tau_int_Q.append(tau_Estimate)
        error_Q.append(error_temp)
        t_comp_list.append(t_comp)
    #print(len(t_comp_list))
    #print(len(tau_int_Q))
    plt.errorbar(t_comp_list,tau_int_Q,error_Q,markersize = 2.0,
                 fmt='o',capsize=2,elinewidth=1,ecolor = colorList[j],color = colorList[j],marker = Plot_marker,
                 markeredgewidth=1,label = r'$\beta = {:.2f}$'.format(beta[j]))#

    #QArray = [np.array(autocorrArray)**2]
    #mean, var, tau_Estimate, error_temp = mod.tauint(QArray,0)
    #tau_int_Qsquared.append(tau_Estimate)
    #error_Qsquared.append(error_temp)

    [print(item) for item in zip(tau_int_Q,error_Q)]
    print("---------------------------------------------")

#plt.xlim(0)
#plt.ylim(100,350)
plt.xlabel(r'$t/a^2$')
plt.ylabel(r'$\tau_{int}$')
plt.legend()

plt.savefig('IntegratedAutoCorrTime_Q.pdf', bbox_inches="tight")