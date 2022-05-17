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
#directory_list = [topC_dir+"beta6_460000\\22X22X22X22\\GF\\autocorr0\\"]
directory_list = [topC_dir+"beta6_260000\\16X16X16X16\\GF\\autocorr0\\"]

fileStart_list = [0]
fileEnd_list = [2000]#CHANGE
beta = [6.0,6.13,6.26,6.46]
V = [16**4,20**4,24**4,32**4]
extrapolate_at_t_comp = [4]
#t0_list = [0.027725,0.027575,0.0278,0.027875]
t0_list = np.array([0.1108,0.1108,0.1109,0.1115])
t0_list *=0.5**2

a = [calc_a(entry) for entry in beta]

j=0
fileStart = fileStart_list[j]
fileEnd = fileEnd_list[j]
directory = directory_list[j]
t_comp = extrapolate_at_t_comp[j]
t0 = t0_list[j]
frames = []
for i in range(fileStart,fileEnd+1,1):
    try:
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory+file,index_col='t',names = ['t',str(i*400)], sep=",", header=None)
        frames.append(df1)
    except:
        print("Warning: file with name " + file + " is absent. Skipping it...")
        continue

df2 = pd.concat(frames,axis=1).loc[t_comp]
plt.figure(figsize=(86,13))
plt.plot(df2)
plt.xticks(np.arange(0,fileEnd_list[j],1),rotation = 90,fontsize = 10)
plt.xlabel("Every increment represents 400 sweeps")
plt.ylabel(r'$Q$')
#plt.legend()
plt.savefig('Q_v_mcTime.pdf', bbox_inches="tight")
