from __future__ import division
### The goal here is to have a script that goes through writing the meshes then running them mpi via command line
import numpy as np
import os
from fenics import *

### initialize mesh attributes
numPores = 12#np.linspace(4,8,5) number of nanopores per side
radius = 0.5#e-9
bR =  np.linspace(3,8,6)#e-9
boxSize = 10#e-9 width of reservoir
length = 20#e-9  pore length
revH = 20#e-9 height of reservoir
nm = 1e-9
#aqD = np.linspace(0.2,0.6,5)
### generate HDF5 files
import nanoporousMesher as mesher
import nanoMesherNoAQ as meshNO

if 1:  # meshes are built, don't run unless need to rebuild
  name = "AQ_NO"
  meshNO.Build(numPores,radius*nm,boxSize*nm,length*nm,revH*nm,name)
  for i,num in enumerate(bR):
    name = "AQ_{}".format(int(num))
    mesher.Build(numPores,radius*nm,num*nm,boxSize*nm,length*nm,revH*nm,name)

### construct the command line commands
myPath = os.path.abspath(__file__)
path = os.path.abspath(os.path.join(myPath,'..'))
path = path+"/"
#path = "/home/AD/srbl226/aqui/"
noAQ = path+"AQ_NO.hdf5"
line = "mpirun -np 20 python noAQP_simulator.py -runner {}".format(noAQ)
print line
os.system(line)

x = np.full((1000,),-7.5e-10)
y = np.full((1000,),-7.5e-10)
zs = np.linspace(-20e-9,30e-9,1000)   ### these are the points along the centerline
points = np.stack((x,y,zs))
points=points.transpose()

mesh = Mesh()
facets = MeshFunction("size_t",mesh)
cells = MeshFunction("size_t",mesh)
hdf5=HDF5File(mesh.mpi_comm(),noAQ,'r')
hdf5.read(mesh,'mesh',False)
hdf5.read(facets,'facets')
hdf5.read(cells,'cells')
hdf5.close()
hdf5_solution=HDF5File(mesh.mpi_comm(),path+"AQ_NO_solution.hdf5",'r')
ele = FiniteElement('CG',mesh.ufl_cell(),1)

V = FunctionSpace(mesh,ele)
#V =
F = Function(V)  ### This is just the FEM fxn instance
hdf5_solution.read(F,"solution")
hdf5_solution.close()


### to get values:
u = np.zeros_like(zs)
for i, point in enumerate(points):
  u[i] = F(point)
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.use('Agg')
plt.figure()
plt.plot(zs,u)
plt.ylabel("z (m)")
plt.xlabel("concentration (mM)")
plt.savefig("figs/AQ_NO.png")
plt.close()


for i, num in enumerate(bR):
        name = path+"AQ_{}.hdf5".format(int(num))
        line = "mpirun -np 20 python AQP_realSimulator.py -runner {} {}".format(name,int(num))
        print line
        os.system(line)
        mesh = Mesh()
        facets = MeshFunction("size_t",mesh)
        cells = MeshFunction("size_t",mesh)
        hdf5=HDF5File(mesh.mpi_comm(),name,'r')
        hdf5.read(mesh,'mesh',False)
        hdf5.read(facets,'facets')
        hdf5.read(cells,'cells')
        hdf5.close()
        print "path of HDF5 file: ", path+"AQ_{}_solution.hdf5".format(int(num))
        tempPath ="AQ_{}_solution.hdf5".format(int(num))
        hdf5_solution=HDF5File(mesh.mpi_comm(),path + tempPath,'r')
        mele = FiniteElement('CG',mesh.ufl_cell(),1)

        V = FunctionSpace(mesh,mele)
        G = Function(V)  ### This is just the FEM fxn instance
        hdf5_solution.read(G,"solution")
        hdf5_solution.close()
        plist=[]
        ulist=[]
        count = 0
        for pt in points:
          try:
            ulist.append(G(pt))
            plist.append(zs[count])
            count+=1
          except RuntimeError:
            ulist.append(0)
            plist.append(zs[count])
            count+=1



        ### to get values:
        #u = np.zeros_like(np.shape(zs))
        #for i, point in enumerate(points):
        #    u[i] = ck_u(point)
        plt.figure()
        plt.plot(zs,ulist)
        plt.ylabel("concentration (mM)")
        plt.xlabel("z (m)")
        plt.savefig("figs/AQ_{}.png".format(num))
        plt.close()
### run the command line commands

### sort through the resulting .txt files

