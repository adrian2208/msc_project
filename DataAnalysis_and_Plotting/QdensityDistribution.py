import os.path
import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
from mpl_toolkits import mplot3d
from scipy import stats
import numpy as np
from sklearn.utils import resample
from mayavi import mlab
import matplotlib.colors as colors

dir = 'C:\\Users\\adria\\Documents\\msc_project\\data\\Observables\\Qdensity\\beta6_000000\\24X24X24X24\\GF\\'
filename = 'torus_extdof1_170_Flowtime3.000000.bin'

with open(dir+filename,'rb') as f:
    data = np.fromfile(f, dtype=np.float64)
xmax = 24
ymax = 24
slice = data[:xmax*ymax]

sliceArray = np.reshape(slice,(xmax,ymax))

#plt.imshow(slicearray,interpolation="none")
#plt.colorbar()
#plt.savefig('QdensityDistributiont_20.pdf')

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
x = range(0,xmax)
y = range(0,ymax)
color_map = plt.get_cmap('viridis')
X, Y= np.meshgrid(x, y)

surf = ax.plot_surface(X, Y,sliceArray,cmap = color_map)
fig.colorbar(surf)
ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')
ax.set_zlabel(r'$q(\vec{x})$')
plt.show()














# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

#z = range(0,24)
#x = range(0,24)
#y = range(0,24)

#X, Y, Z = np.meshgrid(x, y, z)

### Creating figure
#fig = plt.figure()
#ax = plt.axes(projection="3d")

## Creating plot
#color_map = plt.get_cmap('PRGn')
#elev_min= np.max(sliceArray)*2
#elev_max = np.min(sliceArray)*2
#an_array = np.where(abs(sliceArray) < 0.05, 0, sliceArray)
#ax.scatter3D(X, Y, Z, c=an_array,cmap = color_map,norm=MidpointNormalize(midpoint=0))
#plt.savefig('test.pdf')
#plt.show()
