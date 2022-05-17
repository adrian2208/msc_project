import numpy as np
from scipy import stats
from mayavi import mlab


dir = 'C:\\Users\\adria\\Documents\\3D-Qdensity\\'
filename = 'torus_extdof1_3_Flowtime13.000000.bin'

with open(dir+filename,'rb') as f:
    data = np.fromfile(f, dtype=np.float64)
xmax = 32
ymax = 32
zmax = 32
slice = data[:xmax*ymax*zmax]

density = slice
#density = np.reshape(slice,(xmax,ymax,zmax))

mu, sigma = 0, 0.1 
x = list(range(xmax))
y = list(range(ymax))
z = list(range(zmax))
#print(density)
X,Y,Z = np.meshgrid(x,y,z)
x= X.ravel()
y= Y.ravel()
z= Z.ravel()

figure = mlab.figure('DensityPlot')
cutoff = max(abs(density))/2
mask = np.where(abs(density) < cutoff,0,density)
print(mask)




pts = mlab.points3d(x,y,z, density, scale_mode='none', scale_factor=0.5,colormap='seismic')

lut = pts.module_manager.scalar_lut_manager.lut.table.to_array()
zeroSpace = 20
lut[:, -1] = np.concatenate((np.linspace(255, 0, 128-int(zeroSpace/2)),np.zeros(zeroSpace),np.linspace(0,255, 128-int(zeroSpace/2)))) 
pts.module_manager.scalar_lut_manager.lut.table = lut
mlab.axes()
mlab.show()


#figure = mlab.figure('DensityPlot')

#grid = mlab.pipeline.scalar_field(x, y, z, density)
#min = density.min()
#max=density.max()
#mlab.pipeline.volume(grid, vmin=min, vmax=min + .5*(max-min))

#mlab.axes()
#mlab.show()





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
