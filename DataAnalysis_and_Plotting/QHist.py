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


#INPUT FLOW TIME T_COMP AND CORRESPONDING LATTICE CONSTANT
t_comp = 1.5
latticeConstant = 0.0934

x_min = -20 # <- Adjust to spread
x_max = 20  # <-  ^             again..keep constant for all flow time graphs or it's harder to see the change
y_max = 17  # <- set to largest value of y at highest flow time and keep it constant for all flow time graphs

SmearRad = CompFlowTime_to_SmearingRadius(t_comp,latticeConstant)
QArray = df2.loc[t_comp]


NrBins = (x_max-x_min)*2
binsize = (x_max-x_min)/(NrBins-1)
fig = sns.histplot(QArray,bins=NrBins,binwidth=binsize,color = "#D5824B",label=r'$\sqrt{8t}/r_0 =$'+ "{:.2f}".format(SmearRad))
fig.set_xticks(range(x_min,x_max+1,1))
plt.xlim(x_min-0.5*binsize,x_max+0.5*binsize)
plt.ylim(0,yMax)
plt.xlabel("Topological Charge",fontsize = 10)
plt.ylabel("Configuration count",fontsize = 10)
plt.xticks(rotation = 90,fontsize = 8)
plt.yticks(fontsize = 10)
plt.legend()
plt.savefig('Q_histogram_SmearRad_{:.2f}_a_{:.3f}.pdf'.format(SmearRad,latticeConstant))




