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
from sklearn.linear_model import LinearRegression
blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'
purple = '#7047AD'

#Take files from this directory
directory = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\beta6_000000\\16X16X16X16\\GF\\"

def calcE(Earray):
    df = Earray
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 1000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    return mean, std_error
#Turn the files into a single dataframe
directory = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\beta6_000000\\16X16X16X16\\GF\\"
fileStart = 0
fileEnd = 200
a = 0.0934

r_0 = 0.5
t_comp_lim_low = 3.0
t_comp_lim_high = 3.4

frames = []
for i in range(fileStart,fileEnd+1,1):
    file = 'torus_extdof4_'+str(i)+'.csv'
    df1 = pd.read_csv(directory+file,index_col='t',names = ['t','E'], sep=",", header=None)
    frames.append(df1.loc[t_comp_lim_low:t_comp_lim_high])
df2 = pd.concat(frames,axis=0)

frames_model = []
for i in range(fileStart,fileEnd+1,1):
    file_model = 'torus_extdof4_'+str(i)+'.csv'
    df1_model = pd.read_csv(directory+file,index_col='t',names = ['t','E'+str(i)], sep=",", header=None)
    frames_model.append(df1.loc[t_comp_lim_low:t_comp_lim_high])
df2_model = pd.concat(frames,axis=1)

tOverR_0= []
tsquaredE = []
errorList = []
for t_comp in df2.index.tolist()[df2.index.tolist().index(t_comp_lim_low):df2.index.tolist().index(t_comp_lim_high)]:
    Earray = df2_model.loc[t_comp]
    Eval, error = calcE(Earray)
    tOverR_0.append(t_comp*a**2/(r_0**2))
    tsquaredE.append(t_comp**2*Eval)
    errorList.append(t_comp**2*error)

model = LinearRegression()
model.fit(np.array(tOverR_0).reshape(-1,1),np.array(tsquaredE),np.array(errorList))
#plt.errorbar(tOverR_0,tsquaredE,errorList,markersize = 2.0,
#             fmt='o',ecolor = purple,color = purple,capsize=2,elinewidth=1,
#            markeredgewidth=1,label = r'$a =$'+ "{:.3f}".format(a))

#sns.regplot(tOverR_0,tsquaredE,errorList)
sns.regplot(df2.index*a**2/r_0**2,df2.index**2*df2['E'],x_estimator=np.mean,color = purple,label = r"$a=$"+ "{:.3f}".format(a))
#print(df2)


plt.xlabel(r'$t/r_0^2$')
plt.ylabel(r'$t^2\langle E \rangle$')
#plt.legend()
a = model.intercept_
b = model.coef_[0]
t_0_overR_0 = (0.3-a)/b
#plt.plot(np.array(tOverR_0),a+ b*np.array(tOverR_0),color = blue,label = "Linear fit")
#t_prediction = np.linspace(tOverR_0[0],tOverR_0[-1],1000).reshape(-1,1)
#plt.plot(t_prediction,model.predict(t_prediction),color = blue,label = "Linear fit")
plt.axhline(y=0.3,color = green,alpha=0.6,linestyle = 'dashed')
plt.axvline(x=t_0_overR_0,color = green,alpha=0.6,linestyle = 'dashed',label=r"$t_0/r_0^2=$"+ "{:.4f}".format(t_0_overR_0))
plt.legend()
plt.savefig('t_0_linearExtrapolation.pdf')

