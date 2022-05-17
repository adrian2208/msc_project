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
from resizeFig import set_size
plt.figure(figsize=set_size(420))
blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'
def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))

def TopSusc(topChargeArray,latticeConstant,latticeSize):
    df = topChargeArray**2
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 10000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    propagated_error = 197.327*0.25*mean**(-3.0/4.0)/((latticeConstant**4 * latticeSize)**(0.25))*std_error
    mean = (mean/(latticeConstant**4*latticeSize))**(1/4)*197.3
    print("Chi^1/4 = "+str(mean)+"("+str(propagated_error)+") Mev")
    return mean, propagated_error

def CompFlowTime_to_SmearingRadius(t_comp,latticeConstant):
    r_0 = 0.5
    return np.sqrt(8*t_comp*latticeConstant**2)/r_0

#Take files from this directory
topC_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\"
directory_list = [topC_dir+ "beta6_000000\\16X16X16X16\\GF\\",
                  topC_dir+ "beta6_130000\\20X20X20X20\\GF\\",
                  topC_dir+ "beta6_260000\\24X24X24X24\\GF\\",
                  topC_dir+ "beta6_460000\\32X32X32X32\\GF\\Quncorrelated-4050Updates\\"]
fileStart_list = [0,0,0,0]
fileEnd_list = [1703,1323,613,214]#CHANGE
beta_list  = [6.0,6.13,6.26,6.46]
V_list  = [16**4,20**4,24**4,32**4]
color_list = [blue,orange,green,red]
marker_list = ['^','o','s','p']

#RIGHT CLICK FILENAME AND RUN WITHOUT DEBUGGING (VS)
for j in range(len(directory_list)):
    #Take files from this directory
    directory = directory_list[j]
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    a_val = calc_a(beta_list[j])
    V = V_list[j]
    plot_color = color_list[j]
    Plot_marker = marker_list[j]

    #Turn the files into a single dataframe
    frames = []
    for i in range(fileStart,fileEnd+1,1):#97->282
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
        frames.append(df1)
    df2 = pd.concat(frames,axis=1)




    Smear_rad = []
    topological_susceptibility = []
    error = []
    for t_comp in df2.index.tolist():
        QArray = df2.loc[t_comp]
        chi, e = TopSusc(QArray,a_val,V)
        Smear_rad.append(CompFlowTime_to_SmearingRadius(t_comp,a_val))
        topological_susceptibility.append(chi)
        error.append(e)

    plt.errorbar(Smear_rad,topological_susceptibility,error,markersize = 2.0,
                 fmt='o',ecolor = plot_color,color = plot_color,marker = Plot_marker ,capsize=2,elinewidth=1,
                markeredgewidth=1,label = r'$a =$'+ "{:.3f}".format(round(a_val,3)),alpha =0.7)


plt.xlim(0)
plt.ylim(100,350)
plt.xlabel(r'$\frac{\sqrt{8t}}{r_0}$')
plt.ylabel(r'$\chi_t^{1/4}[MeV]$')
plt.legend()
plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\chi_fourthRoot_v_smearRad.pdf', bbox_inches="tight")




