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
purple = '#7047AD'
def calc_a(beta):
    return 0.5 * np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 
    0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0))


def CompFlowTime_to_SmearingRadius(t_comp,latticeConstant):
    r_0 = 0.5
    return np.sqrt(8*t_comp*latticeConstant**2)/r_0

topC_dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\"
directory_list = [topC_dir+ "beta6_000000\\16X16X16X16\\GF\\",
                  topC_dir+ "beta6_130000\\20X20X20X20\\GF\\",
                  topC_dir+ "beta6_260000\\24X24X24X24\\GF\\",
                  topC_dir+ "beta6_460000\\32X32X32X32\\GF\\Quncorrelated-4050Updates\\"]

fileStart_list = [0,0,0,0]#50
fileEnd_list = [1703,1323,613,214]#CHANGE
#colorList = [green,orange,blue,purple]

beta = [6.0,6.13,6.26,6.46]
t0_frac = "t0_over12"
#[3.2,4.85,7.1,12.2]#t0
#[0.8,1.2,1.8,3.0],#1/4
#[0.55,0.8,1.2,2.0],#1/6
#[0.25,0.4,0.6,1.0]#1/12
plot_at_t_comp = [0.25,0.4,0.6,1.0]

a = [calc_a(entry) for entry in beta]
colorList = [green,orange,blue,purple]
#Turn the files into a single dataframe

fig, axs = plt.subplots(nrows=4, sharex=True, subplot_kw=dict(frameon=False))
plt.subplots_adjust(hspace=.5)

#Take files from this directory
step = 0.125
xmin = -6
xmax = 6
binarray = np.arange(xmin-step/2,xmax+step+step/2,step )
for j in range(len(fileStart_list)):
    directory = directory_list[j]
    fileStart = fileStart_list[j]
    fileEnd = fileEnd_list[j]
    t_comp = plot_at_t_comp[j]
    latticeConstant = a[j]
    #plot_color = color_list[j]

    frames = []
    for i in range(fileStart,fileEnd+1,1):#4
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
        frames.append(df1)
    df2 = pd.concat(frames,axis=1)



    #INPUT FLOW TIME T_COMP AND CORRESPONDING LATTICE CONSTANT

    x_min = -6 # <- Adjust to spread
    x_max = 6  # <-  ^             again..keep constant for all flow time graphs or it's harder to see the change
    y_max = 17  # <- set to largest value of y at highest flow time and keep it constant for all flow time graphs

    SmearRad = CompFlowTime_to_SmearingRadius(t_comp,latticeConstant)
    QArray = df2.loc[t_comp]


    NrBins = (x_max-x_min)*2
    binsize = (x_max-x_min)/(NrBins-1)
    #plt.subplot(len(fileStart_list),1,j+1)
    #axs[j].histplot(QArray,bins=NrBins,binwidth=binsize,color = "#D5824B")
    axs[j].hist(QArray,bins=binarray,density = True,label = r'$\beta = {}$'.format(beta[j]),color = colorList[j])
    axs[j].legend()
    #

#print(np.arange(x_min,x_max,0.25))
plt.xlim(x_min-0.5*binsize,x_max+0.5*binsize)
#plt.ylim(0,yMax)
plt.xlabel("Topological Charge",fontsize = 10)
#plt.ylabel("Configuration count",fontsize = 10)
plt.xticks(range(xmin,xmax+1,1),rotation = 90,fontsize = 8)
#plt.yticks(fontsize = 10)

plt.savefig('Q_histogram_at_' + t0_frac +'.pdf')




