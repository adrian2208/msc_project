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
from sklearn.linear_model import LinearRegression
from resizeFig import set_size
plt.figure(figsize=set_size(420))
blue = '#447DAC'
green = '#30988A'
orange = '#D5824B'
red = '#DF4545'
purple = '#7047AD'

def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))



def calcE(Earray):
    df = Earray
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 10000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    return mean, std_error


root = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\"
dirs = [root + "beta6_000000\\16X16X16X16\\GF\\", 
        root + "beta6_130000\\20X20X20X20\\GF\\", 
        root + "beta6_260000\\24X24X24X24\\GF\\",
        root + "beta6_460000\\32X32X32X32\\GF\\"]#change
fileStart_list = [0,1,0,0] #change
fileEnd_list  = [1000,1000,309,701]#change
beta_list  = [6.0,6.13,6.26,6.46] #change
V_list  = [16**4,20**4,24**4,32**4] #change
t_comp_lim_low_list = [3.0,4.5,6.5,11.5]
t_comp_lim_high_list = [3.4,5.05,8,13]
color_list = [blue,orange,green,purple]
marker_list = ['^','o','s','p']



#RIGHT CLICK FILENAME AND RUN WITHOUT DEBUGGING (VS)
for j in range(len(dirs)):
    directory = dirs[j]
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    a = calc_a(beta_list[j])
    plot_color = color_list[j]

    r_0 = 0.5
    t_comp_lim_low = t_comp_lim_low_list[j]
    t_comp_lim_high =t_comp_lim_high_list[j]
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
    sns.regplot(df2.index*a**2/r_0**2,df2.index**2*df2['E'],x_estimator=np.mean,color = plot_color,label = r"$a=$"+ "{:.3f}".format(a),n_boot=10000)
    #print(df2)



    #plt.legend()
    a = model.intercept_
    b = model.coef_[0]
    t_0_overR_0 = (0.3-a)/b
    #plt.plot(np.array(tOverR_0),a+ b*np.array(tOverR_0),color = blue,label = "Linear fit")
    #t_prediction = np.linspace(tOverR_0[0],tOverR_0[-1],1000).reshape(-1,1)
    #plt.plot(t_prediction,model.predict(t_prediction),color = blue,label = "Linear fit")
    
    plt.axvline(x=t_0_overR_0,color = plot_color,alpha=0.6,linestyle = 'dashed',label=r"$t_0/r_0^2=$"+ "{:.4f}".format(t_0_overR_0))
plt.axhline(y=0.3,color = 'grey',alpha=0.6,linestyle = 'dashed')
plt.xlabel(r'$t/r_0^2$')
plt.ylabel(r'$t^2\langle E \rangle$')
plt.ylim(0.2,0.4)
#plt.xlim(0.11,0.113)
plt.legend()
plt.savefig('C:\\Users\\adria\\Documents\\msc_project\\doc\\TsquaredE_fit.pdf', bbox_inches="tight")

