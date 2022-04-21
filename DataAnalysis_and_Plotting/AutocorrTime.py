import module1 as mod
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')

#data = [[data]]

r_0 = 0.5
latticeConstant = [0.1]
#Take files from this directory
directory = ["C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Topological_Charge\\beta5_960000\\10X10X10X10\\GF\\"]
NrFiles = [500]
#Turn the files into a single dataframe
for j in range(len(NrFiles)):
    frames = []
    for i in range(0,NrFiles[j]+1,1):
        file = 'torus_extdof4_'+str(i)+'.csv'
        df1 = pd.read_csv(directory[j]+file,index_col='t',names = ['t','Q'+str(i)], sep=",", header=None)
        frames.append(df1)
    df2 = pd.concat(frames,axis=1)

    t0_overR_0Squared = []
    tau_int = []
    error = []
    for t_comp in df2.index.tolist():
        QArray = [[df2.loc[t_comp].to_numpy()]]
        mean, var, tau_Estimate, error_temp = mod.tauint(QArray,0)
        t0_overR_0Squared.append(t_comp*latticeConstant[j]**2/(r_0**2))
        tau_int.append(tau_Estimate)
        error.append(error_temp)

    plt.errorbar(t0_overR_0Squared,tau_int,error,markersize = 2.0,
                 fmt='o',ecolor = '#30988A',color = '#30988A',capsize=2,elinewidth=1,
                markeredgewidth=1,label = r'$a =$'+ "{:.3f}".format(latticeConstant[j]))
#plt.xlim(0)
#plt.ylim(100,350)
plt.xlabel(r'$t/r_0^2$')
plt.ylabel(r'$\tau_{int}$')

plt.savefig('IntAutoCorrTime_v_Flowtime.pdf', bbox_inches="tight")