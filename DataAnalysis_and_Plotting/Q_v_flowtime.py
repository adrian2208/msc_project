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
#plt.figure(figsize=set_size(420))
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

root = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\"
dirs = [root + "beta6_000000\\16X16X16X16\\GF\\", 
        root + "beta6_130000\\20X20X20X20\\GF\\", 
        root + "beta6_260000\\24X24X24X24\\GF\\",
        root + "beta6_460000\\32X32X32X32\\GF\\For_susceptibility\\tau_int_leq0.5-2400updates\\"]
fileStart_list = [0,0,0,0]
fileEnd_list = [20,20,20,20]
beta_list  = [6.0,6.13,6.26,6.46]
V_list  = [16**4,20**4,24**4,32**4]
t0OveraList = [3.2,4.85,7.1,12.2]
color_list = [blue,orange,green,red]

fig, axs = plt.subplots(nrows=4, sharex=True, subplot_kw=dict(frameon=False),figsize=(5.535,8))
plt.subplots_adjust(hspace=.1)

#RIGHT CLICK FILENAME AND RUN WITHOUT DEBUGGING (VS)
for j in range(len(dirs)):
    #Take files from this directory
    directory = dirs[j]
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    a_val = calc_a(beta_list[j])
    V = V_list[j]
    plot_color = color_list[j]

    #Turn the files into a single dataframe
    frames = []
    for i in range(fileStart,fileEnd+1,1):#97->282
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
        frames.append(df1)
    df2 = pd.concat(frames,axis=1)

    axs[j].plot(df2)
    if (j==0):
        axs[j].axvline(x=t0OveraList[j],color = red,linestyle = 'dashed',label=r"$t_0/a^2$")
        axs[j].legend()
    else:
        axs[j].axvline(x=t0OveraList[j],color = red,linestyle = 'dashed')
    ymax = np.ceil(max(abs(df2.iloc[-1])))+1
    #ymin = np.floor(min(df2.iloc[-1]))-1
    axs[j].grid(which='both',color='grey', linestyle='-', linewidth=0.5,axis='y')
    axs[j].set_yticks(np.arange(-ymax,ymax+1,1))
    axs[j].set_ylim(-ymax-4,ymax+4)
    axs[j].tick_params(axis='both', which='major', labelsize=7)
    

plt.xlabel(r'$t/a^2$')
plt.ylabel(r'$Q$') 
plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\FINALPLOTS\\Q_v_flowtime.pdf', bbox_inches="tight")
plt.clf()

#plt.xlim(0)
#plt.ylim(100,350)




