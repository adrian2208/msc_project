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
#directory = '.\\data\\Observables\\Topological_Charge\\beta6_450000\\8X8X8X8\\'
directory = '.\\data\\Observables\\Topological_Charge\\beta6_000000\\24X24X24X24\\'

def getPlotDirectory(dataDirectory):
    plotDir = '.\\Plots\\'+dataDirectory[dataDirectory.find('data'):]
    Path(plotDir).mkdir(parents=True, exist_ok=True)
    return plotDir

def ObservableDirectory(observable,beta,lattice_shape):
    betaStr = str(beta).split('.')[0]+'_'+str(beta).split('.')[1]
    while len(betaStr)< 8:
        betaStr += '0'
    directory = 'C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\'+observable + '\\beta' +betaStr+'\\'+ lattice_shape +'\\'
    if not os.path.isdir(directory):
        print("No such directory exists:")
        print(directory)
    return directory


def DirectoryToDataframe(directory):
    frames = []
    name_list = []
    i = 0
    for file in os.listdir(directory):
        df1 = pd.read_csv(directory+file,names = [str(i)], sep="\n", header=None)
        frames.append(df1)
        name_list.append(str(i))
        i+=1
    result = pd.concat(frames, axis=1)
    return result





#for entry in name_list:
#   fig = sns.lineplot(x=range(0,3453), y=entry, data=result)
#fig.set_yticks(range(-10,10,1))
#fig.set(ylim=(-10,10))
#plt.savefig('24_beta6.pdf')

def HistogramFromDataFrameRow(dataframe,row,x_min,x_max,Nr_bins):
    data = dataframe.iloc[row]
    mean = data.mean()
    std = data.std()
    print("mean = " + str(mean))
    print("sigma^2 = " + str(std**2))
    x = np.linspace(x_min,x_max,200)
    p = stats.norm.pdf(x,mean,std)
    
    val_width = x_max-x_min
    bin_width = val_width/Nr_bins
    fig = sns.histplot(data,binwidth=bin_width,bins = Nr_bins,binrange = (x_min,x_max),color = "grey")
    fig.set_xticks(range(x_min,x_max,1))
    fig = plt.plot(x,len(data)*p)
    textstr = '\n'.join((
    r'$\mu=%.2f$' % (mean, ),
    r'$\sigma=%.2f$' % (std**2, )))
    props = dict(boxstyle='round', facecolor='0.9', alpha=1.0)
    plt.text(7.0,11.5,textstr, fontsize=12,
        verticalalignment='top',bbox=props)
    plt.xlabel("Topological Charge",)
    plt.ylabel("Configuration count")
    plt.xticks(rotation = 90)
    plt.savefig('Hist_Q_atZeroFlowtime.pdf', bbox_inches="tight")
    #plt.savefig('test.pdf')


#HistogramFromDataFrameRow(DirectoryToDataframe("C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta6_000000\\24X24X24X24\\GF\\"),0,-13,13,103)

#AUTOCORRELATION TIME
#data = DirectoryToDataframe("C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta6_000000\\24X24X24X24\\GF\\").iloc[0].to_numpy()
#data = np.random.normal(0.0,4.58,104)
#data = [[data]]
#print(mod.tauint(data,0))
#print(mod.gamma(data,0))


def TopSusc(topChargeArray):
    df = topChargeArray**2
    result= stats.bootstrap([df.to_numpy()],np.mean,n_resamples = 1000,confidence_level=0.0)
    mean = result.confidence_interval[0]
    std_error = result.standard_error
    propagated_error = 197.3*0.25*mean**(-3.0/4.0)/((latticeConstant**4 * latticeSize)**(0.25))*std_error
    mean = (mean/(latticeConstant**4*latticeSize))**(1/4)*197.3
    print("Chi^1/4 = "+str(mean)+"("+str(propagated_error)+") Mev")
#TopSusc("C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta6_000000\\24X24X24X24\\GF\\",24**4,0.09314)


#df1 = pd.read_csv("C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\EnergyDensity\\beta6_000000\\8X8X8X8\\torus_extdof4ensemble_0LHMC.csv",names = [0], sep="\n", header=None)
#fig = sns.lineplot(data=df1)
#plt.savefig('EDvLHMCstep.pdf')
#plt.close()

dataPath ="C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta6_000000\\24X24X24X24\\GF\\" #ObservableDirectory('Topological_Charge',6.0,'24X24X24X24')
#outDir = getPlotDirectory(dataPath)
#for i in range(16):
#    fileName = "torus_extdof4Heatbath_4_ORperHB_Flowed"+ str(i) + ".csv"
#    df1 = pd.read_csv(dataPath+fileName,names = [0], sep="\n", header=None)
#    fig = sns.lineplot(data=df1,legend=None)
#    fig.set(yticks=list(range(-12,12,1)))
#plt.savefig('FLOWv.pdf')



#df = DirectoryToDataframe(dataPath)
#fig = sns.lineplot(data=df,legend=None)
#fig.set(yticks=list(range(-12,12,1)))
#plt.savefig('FLOWv.pdf')


def plotFlowedTopologicalChargeOverMcStep(dir):
    df = DirectoryToDataframe(dir)
    print(df.iloc[800])
    fig = sns.lineplot(data=df.iloc[800])
    plt.xticks(rotation = 90)
    plt.savefig('flowedQvHeatbathStep.pdf')

#plotFlowedTopologicalChargeOverMcStep("C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta6_000000\\24X24X24X24\\GF\\")

#df1 = pd.read_csv("C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta6_000000\\8X8X8X8\\GF\\torus_extdof4Heatbath_4_ORperHB83flowed.csv",names = ["a"], sep="\n", header=None)
#fig = sns.lineplot(data=df1)
#plt.xticks(rotation = 90)
#plt.savefig('flowedQvHeatbathStep_no_therm.pdf')

def changeFileNames():
    dir = "C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta6_000000\\24X24X24X24\\GF\\"
    for file in os.listdir(dir):
        tempFileName = file[:13]+"_"+file[13:]
        #tempFileName = file[:38]+"0"+file[38:]
        print("will be renamed to:")
        print(tempFileName)
    val = str(input("accept?"))
    if (val == "yes"):
        for file in os.listdir(dir):
            tempFileName = file[:13]+"_"+file[13:]
            os.rename(dir+file,dir +tempFileName )

