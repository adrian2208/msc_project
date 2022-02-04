import os.path
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')

#directory = '.\\data\\Observables\\Topological_Charge\\beta6_450000\\8X8X8X8\\'
directory = '.\\data\\Observables\\Topological_Charge\\beta6_000000\\24X24X24X24\\'

def ObservableDirectory(observable,beta,lattice_shape):
    betaStr = str(beta).split('.')[0]+'_'+str(beta).split('.')[1]
    while len(betaStr)< 8:
        betaStr += '0'
    directory = '.\\data\\Observables\\'+observable + '\\beta' +betaStr+'\\'+ lattice_shape +'\\'
    if not os.path.isdir(directory):
        print("No such directory exists:")
        print(directory)
    return directory

ObservableDirectory('Topological_Charge',6.0,'12X12X12X12')
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
    val_width = x_max-x_min
    bin_width = val_width/Nr_bins
    fig = sns.histplot(0.25*dataframe.iloc[row],binwidth=bin_width,bins = Nr_bins,binrange = (x_min,x_max))
    fig.set_xticks(range(x_min,x_max,1))
    plt.xticks(rotation = 45)
    plt.savefig('Hist_zero_flowTIme.pdf')

HistogramFromDataFrameRow(DirectoryToDataframe(ObservableDirectory('Topological_Charge',6.0,'12X12X12X12')),599,-15,15,61)