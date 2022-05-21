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
plt.figure(figsize=set_size(300))
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
directory_list = [topC_dir+ "beta6_000000\\12X12X12X12\\GF\\",
                  topC_dir+ "beta6_130000\\14X14X14X14\\GF\\",
                  topC_dir+ "beta6_260000\\16X16X16X16\\GF\\autocorr0\\",
                  topC_dir+ "beta6_460000\\22X22X22X22\\GF\\autocorr0\\"]

fileStart_list = [0,0,0,0]#50
fileEnd_list = [2000,2000,2000,3000]#CHANGE

beta = [6.0,6.13,6.26,6.46]
V = [16**4,20**4,24**4,32**4]
extrapolate_at_t_comp = [3.5,5.5,4,4]
#t0_list = [0.027725,0.027575,0.0278,0.027875]
t0_list = np.array([0.1108,0.1108,0.1109,0.1115])
t0_list *=0.5**2
blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'
purple = '#7047AD'
color_list = [blue,green,orange,purple,red]
a = [calc_a(entry) for entry in beta]
label_list = [r"$a=$"+ "{:.3f}".format(a[0]),r"$a=$"+ "{:.3f}".format(a[1]),r"$a=$"+ "{:.3f}".format(a[2]),r"$a=$"+ "{:.3f}".format(a[3])]


fig, axs = plt.subplots(nrows=4, sharex=True, subplot_kw=dict(frameon=False))
plt.subplots_adjust(hspace=.3)

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
            df1 = pd.read_csv(directory+file,index_col='t',names = ['t',str(i)], sep=",", header=None)
            frames.append(df1)
        except:
            print("Warning: file with name " + file + " is absent. Skipping it...")
            continue


    df2 = pd.concat(frames,axis=1).loc[t_comp]
    axs[j].plot(df2,label = label_list[j],color=color_list[j])
    axs[j].set_yticks(np.arange(-3,3+1,1))
    axs[j].tick_params(axis='both', which='major', labelsize=7)
    axs[j].grid(which='both',color='grey', linestyle='-', linewidth=0.5,axis='y')
    axs[j].legend()


#plt.figure(figsize=(86,13))

plt.xticks(np.arange(0,3200,200),rotation = 90,fontsize = 10)


plt.xlabel("MC time")
plt.ylabel(r'$Q$')
#plt.legend()
plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\FINALPLOTS\\Q_v_mcTime.pdf', bbox_inches="tight")
