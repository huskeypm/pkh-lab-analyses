import numpy as np
import matplotlib.pylab as plt
#import loos
#import PyTraj

#import binfuncs as bf

#cols = [""]

def BinTimeSeries(data,binMin,binMax,res=1.,plot=True):
    # define bins 
    nBins = np.int((1/res)*(binMax-binMin)+1)
    bins = np.linspace(binMin, binMax, nBins)

    #digitize into bins 
    digitized = np.digitize(data, bins)

    #here I just want to output distance belown 12A
    ratio = 12/binMax
    newbins = np.int(nBins*ratio)

    # create 2D matrix of bin versus time
    nTs = np.shape(data)[0]
    z=np.zeros((newbins,nTs))
    for i in range(newbins):
      z[i,np.where(digitized==i)]=1 # here is to locate the pot in time series
    

    # create 2D coordinates and plot
    grid_t,grid_r = np.mgrid[0:(nTs-1):(nTs*1j),binMin:12:(newbins*1j)]
    #np.shape(grid_x)
    #pcolormesh(z)
    if plot:
      plt.figure()
      plt.pcolormesh(grid_t,grid_r,z.transpose())
    return grid_t,grid_r,z


#binMin=0.; binMax=10.
#res = 0.5

#zs=[]
#for i in range(nums):
#  grid_t,grid_r,z = BinTimeSeries(allDists[i,:],binMin,binMax)#,plot=False)
#  zs.append(z)
