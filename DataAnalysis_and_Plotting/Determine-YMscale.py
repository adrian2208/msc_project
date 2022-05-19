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
plt.figure(figsize=set_size(420))
blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'



df1 = pd.read_csv("ChiSquaredFit.csv",names = ['Lambda','Weight'], sep=",", header=None)
LambdaList = df1['Lambda']
OneOverChiSquared_list = df1['Weight']

Lambda_min = min(LambdaList)
Lambda_max = max(LambdaList)
Nr_Lambdas = len(LambdaList)

plt.plot(LambdaList,OneOverChiSquared_list)
plt.xlabel(r'$\Lambda_{YM}[MeV]$')
#plt.ylabel(r'$1/\chi^2$')
#plt.legend()
#plt.ylim(0,0.25)
plt.xlim(Lambda_min,Lambda_max)

plt.savefig('t_squared_v_q_ChiSquared_fitDistribution.pdf', bbox_inches="tight")
plt.clf()
plt.hist(LambdaList,bins = int(len(LambdaList)/4),weights = OneOverChiSquared_list,color = blue)
plt.xlabel(r'$\Lambda_{YM}[MeV]$')
#plt.ylabel(r'$1/\chi^2$')
#plt.legend()
#plt.ylim(0,0.25)
#plt.xlim(Lambda_min,Lambda_max)

#df1.sort_values('Weight', inplace=True)
cumsum = df1.Weight.cumsum()
cutoff = df1.Weight.sum() / 2.0
#confidence_cutoff = df1.Weight.sum()*0.68
#confidence_cutoff2 = df1.Weight.sum()*0.32
confidence_cutoff = df1.Weight.sum()*0.84
confidence_cutoff2 = df1.Weight.sum()*0.16
median = df1.Lambda[cumsum >= cutoff].iloc[0]
upperCI = df1.Lambda[cumsum >= confidence_cutoff2].iloc[0]
lowerCI = df1.Lambda[cumsum >= confidence_cutoff].iloc[0]
n = len(cumsum)
z = 0.468
#upperCI = np.ceil(n*0.5 - z*np.sqrt(n*0.5*(1-0.5)))
#lowerCI = np.ceil(n*0.5 + z*np.sqrt(n*0.5*(1-0.5)))
print(median)
print(upperCI)
print(lowerCI)


plt.axvline(x=median,color = red,linestyle = 'dashed',label=r"$median$")
plt.axvline(x=upperCI,color = orange,linestyle = 'dashed',label=r"$68\% CI$")
plt.axvline(x=lowerCI,color = orange,linestyle = 'dashed')
plt.legend()
plt.savefig('t_squared_v_q_ChiSquared_fitDistribution_HIST.pdf', bbox_inches="tight")
#C:\\Users\\adria\\Documents\\msc_project\\doc\\

#print(np.sum(df1["Weight"]))


#print(df1.Lambda.iloc[int(upperCI)])
#print(df1.Lambda.iloc[int(lowerCI)])