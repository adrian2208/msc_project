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

def TopSusc(topChargeArray,latticeConstant,latticeSize):
    df = topChargeArray**2
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 1000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    propagated_error = 197.3*0.25*mean**(-3.0/4.0)/((latticeConstant**4 * latticeSize)**(0.25))*std_error
    mean = (mean/(latticeConstant**4*latticeSize))**(1/4)*197.3
    print("Chi^1/4 = "+str(mean)+"("+str(propagated_error)+") Mev")
    return mean, propagated_error

def CompFlowTime_to_SmearingRadius(t_comp,latticeConstant):
    r_0 = 0.5
    return np.sqrt(8*t_comp*latticeConstant**2)/r_0


#Take files from this directory
directory = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta6_000000\\24X24X24X24\\GF\\"

#Turn the files into a single dataframe
frames = []
for i in range(97,282,1):#97->282
    file = 'torus_extdof4_'+str(i)+'.csv'
    df1 = pd.read_csv(directory+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
    frames.append(df1)
df2 = pd.concat(frames,axis=1)


a = 0.0934
V = 24**4

Smear_rad = []
topological_susceptibility = []
error = []
for t_comp in df2.index.tolist():
    QArray = df2.loc[t_comp]
    chi, e = TopSusc(QArray,a,V)
    Smear_rad.append(CompFlowTime_to_SmearingRadius(t_comp,a))
    topological_susceptibility.append(chi)
    error.append(e)

plt.errorbar(Smear_rad,topological_susceptibility,error,markersize = 2.0,
             fmt='o',ecolor = '#30988A',color = '#30988A',capsize=2,elinewidth=1,
            markeredgewidth=1,label = r'$a =$'+ "{:.3f}".format(a))
plt.xlim(0)
plt.ylim(100,350)
plt.xlabel(r'$\frac{\sqrt{8t}}{r_0}$')
plt.ylabel(r'$\chi_t^{1/4}[MeV]$')

plt.savefig('chi_fourthRoot_v_smearRad.pdf', bbox_inches="tight")




