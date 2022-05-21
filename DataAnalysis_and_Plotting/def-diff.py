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
plt.figure(figsize=set_size(420))
blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'

def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))

#root = "/mnt/gs18/scratch/users/f0103741/Observables/Topological_Charge/"
#root6_46 = "/mnt/gs18/scratch/users/f0103741/6_46_ensemble/1/Observables/Topological_Charge/"
#dirs = [root + "beta6_000000/16X16X16X16/GF/", 
#        root + "beta6_130000/20X20X20X20/GF/", 
#        root + "beta6_260000/24X24X24X24/GF/",
#        root6_46 + "beta6_460000/32X32X32X32/GF/"]

#rootUnimproved = "/mnt/gs18/scratch/users/f0103741/Observables/UnImprovedTopological_Charge/"
#rootUnimproved6_46 = "/mnt/gs18/scratch/users/f0103741/6_46_ensemble/1/Observables/UnImprovedTopological_Charge/"
#dirsUnimproved = [rootUnimproved + "beta6_000000/16X16X16X16/GF/", 
#        rootUnimproved + "beta6_130000/20X20X20X20/GF/", 
#        rootUnimproved + "beta6_260000/24X24X24X24/GF/",
#        rootUnimproved6_46 + "beta6_460000/32X32X32X32/GF/"]
        
        
root = "C:\\Users\\adria\\Documents\\def-diff\\ImprovedEnergyDensity\\"
dirs = [root + "beta6_000000\\16X16X16X16\\GF\\", 
        root + "beta6_130000\\20X20X20X20\\GF\\", 
        root + "beta6_260000\\24X24X24X24\\GF\\",
        root + "beta6_460000\\32X32X32X32\\GF\\"]

rootUnimproved = "C:\\Users\\adria\\Documents\\def-diff\\EnergyDensity\\"
dirsUnimproved = [rootUnimproved + "beta6_000000\\16X16X16X16\\GF\\", 
        rootUnimproved + "beta6_130000\\20X20X20X20\\GF\\", 
        rootUnimproved + "beta6_260000\\24X24X24X24\\GF\\",
        rootUnimproved + "beta6_460000\\32X32X32X32\\GF\\"]

fileStart_list = [1001,1001,1001,0]
fileEnd_list = [1907,1907,1907,128]
beta_list  = [6.0,6.13,6.26,6.46]
aoverR0 = [calc_a(entry)/0.5 for entry in beta_list]
t_comp_list = [3.5,5.5,8.5,13.0]
color_list = [blue,orange,green,red]

#RIGHT CLICK FILENAME AND RUN WITHOUT DEBUGGING (VS)
difflist = []
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
        if (len(df1) == 172):
            df1 = df1.iloc[0:86]
            df1.to_csv(directory+file,header=None)
    df2 = pd.concat(frames,axis=1).loc[t_comp].to_numpy()

    frames_unimproved = []
    for i in range(fileStart,fileEnd+1,1):
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory_unimproved+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
        frames_unimproved.append(df1)
    df2_unimproved = pd.concat(frames_unimproved,axis=1).loc[t_comp].to_numpy()

    diff = np.mean(abs(df2-df2_unimproved))
    meanList.append(np.mean(df2))
    difflist.append(diff/np.mean(df2))
print(meanList)

plt.scatter(np.array(aoverR0),np.array(difflist))
plt.xlabel(r'$a/r_0$')
plt.ylabel(r'$\frac{\left\langle|E_{c}-E_{i}|\right\rangle}{\left\langle E_{i} \right\rangle}$',rotation = 0,fontsize = 10,loc='center')
plt.savefig('E_def_diff.pdf', bbox_inches="tight")


