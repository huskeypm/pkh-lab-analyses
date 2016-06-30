import matplotlib.pylab as plt 
import numpy as np
from scipy.interpolate import griddata
from dolfin import * 

class Omega1(SubDomain):
    def inside(self, x, on_boundary):
        result = (x[0] >= self.xbounds[0] and x[0] <= self.xbounds[1])
        return result
class Omega0(SubDomain):
    def inside(self, x, on_boundary):
        return True

# Pickle reading writing 
import cPickle as pickle

def writePickle(passingLocs,passingRads,dxs=None,volFracs=None,fileName="data.pkl"):
  randomSamples = dict()
  randomSamples['passingLocs'] = passingLocs
  randomSamples['passingRads'] = passingRads
  randomSamples['dxs'] = dxs
  randomSamples['volFracs'] = volFracs    
  output = open(fileName, 'wb')
  pickle.dump(randomSamples, output)
  output.close()
    
    
def readPickle(fileName):
  pkl_file = open(fileName, 'rb')
  data1 = pickle.load(pkl_file) 

  #print data1['D']
  passingLocsi = data1['passingLocs']
  passingRadsi = data1['passingRads']
  try:   
    dxs = data1['dxs'] 
  except: 
    1    
  try:
    volFracs = data1['volFracs'] 
  except:
    1    
    
  pkl_file.close()

  return passingLocsi,passingRadsi,data1


def integrated(mesh,subdomains,x,xbounds,verbose=False):
    
    # tag 'slice' 
    subdomains0  = Omega0()
    subdomains0.mark(subdomains, 0)
    subdomains1  = Omega1()
    subdomains1.xbounds = xbounds
    marker = 1
    subdomains1.mark(subdomains, marker)
    dx = Measure('dx')[subdomains]

    # get vol of slice 
    #vol = assemble(Constant(1.)*dx(marker,domain=mesh))  
    coords = mesh.coordinates()
    bounds = np.max(coords,axis=0) - np.min(coords,axis=0)
    bounds[0] = xbounds[1]-xbounds[0]
    vol = np.prod(bounds)

    # assemble 'x' over slice 
    X = assemble(x*dx(marker,domain=mesh))
    conc = X/vol
    
    if verbose:
      print "MU Vol", assemble(Constant(1.)*dx(marker,domain=mesh))
      print "MU Conc: %f" % (conc)
    return conc

  #xbounds = [0.9,1.00]
  #integrated(mesh,x,xbounds)



def plotMyLattice(locs,s=50,locs2=[],c1='k',c2='k',s2=50,marker2='o',\
                  title="",name="test.png", figName="NONE"):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    #ax.axis('off')
    plt.tick_params(\
        #axis='x',          # changes apply to the x-axis
        #which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        left='off',         # ticks along the top edge are off
        right='off',
        labelleft='off',
        labelbottom='off') 
    
    
    ax.scatter(locs[:,0],locs[:,1],s=s,color=c1)#,facecolor="none")
    if np.shape(locs2)[0] > 0:
        ax.scatter(locs2[:,0],locs2[:,1],color=c2,s=s2, marker=marker2)#,facecolor="none")
    ax.set_title(title,fontsize=30)
    ax.set_aspect('equal')
    if figName!="NONE":
      plt.gcf().savefig(figName)


def StoreMesh(results):
    mesh = results.mesh
    dims = np.max(mesh.coordinates(),axis=0) - np.min(mesh.coordinates(),axis=0)
    #u = results.u_n.split()[0]
    #u = results.u_n.split()[0]
    up = project(u,FunctionSpace(mesh,"CG",1))
    res = 100
    (gx,gy,gz) = np.mgrid[0:dims[0]:(res*1j),
                          dims[1]/2.:dims[1]/2.:1j,
                          0:dims[2]:(res*1j)]
    img0 = griddata(mesh.coordinates(),up.vector(),(gx,gy,gz))
    return img0

def Store2DMesh(mesh,u):
#    mesh = results.mesh
    #dims = np.max(mesh.coordinates(),axis=0) - np.min(mesh.coordinates(),axis=0)
    mmin = np.min(mesh.coordinates(),axis=0)
    mmax = np.max(mesh.coordinates(),axis=0)

    #u = results.u_n.split()[0]
    #u = results.u_n
    up = project(u,FunctionSpace(mesh,"CG",1))
    res = 100
    #(gx,gy,gz) = np.mgrid[0:dims[0]:(res*1j),
    #                      dims[1]/2.:dims[1]/2.:1j,
    #                      0:dims[2]:(res*1j)]
    (gx,gy) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
                       mmin[0]:mmax[1]:(res*1j)]
    from scipy.interpolate import griddata
    img0 = griddata(mesh.coordinates(),up.vector(),(gx,gy))
    return img0


def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]



# coded to help ensure first two groups have even numbers and the last group gets shafted
def makeGroups(l,nGrps=3):
  nSpheres = np.shape(l)[0]
  grpSize =np.int(np.ceil(nSpheres/np.float(nGrps))) 
  print grpSize
  groups =  chunks(l,grpSize)
  listSums = [len(j) for i,j in enumerate(groups)]
  print "N in group ", listSums, " tot ", np.sum(listSums)
  return groups



# distributes the list of 'M' input potentials over 'N' indexedpoints (M<N) 
def DefinePotentials(idxs,Vs):
    print Vs
    nGrps = np.shape(Vs)[0]
    grps = makeGroups(idxs,nGrps=nGrps)
    nPts  = np.shape(idxs)[0]
    Vpts = np.zeros(nPts)      
    for i,Vi in enumerate(Vs):
        Vpts[ grps[i] ] = Vi

    #np.shape(grps[0]
    #Vlocs = np.reshape(Vpts,[perSide,perSide])
    #plt.pcolormesh(Vlocs)
    return Vpts 

import scipy.spatial.distance as distance

# returns true if entries are clashing 
def CheckClash(locs,rads,eps=0,idx=False):
  n=np.shape(rads)[0]
  #print n
  #print locs
  #print rads
  #print np.ones(n)
  #print rads  
    
  minDist = np.outer(np.ones(n),rads)
  minDist = minDist +  minDist.T + eps   
  #print minDist

  dist = distance.pdist(locs); 
  squareDist = distance.squareform(dist)  

  for i in range(n):
      squareDist[i,i]=1e9

  #print squareDist


  clashMatr= ((squareDist - minDist) < 0).astype(int)
  #print "clash",clashMatr        

  if idx:
    bad = np.sum(clashMatr,axis=0)
    return np.argwhere(bad>0)
  else:
    return np.sum(clashMatr)>0

# http://stackoverflow.com/questions/9081553/python-scatter-plot-size-and-style-of-the-marker/24567352#24567352
def circles(x, y, s, c='b', ax=None, vmin=None, vmax=None, **kwargs):

    """
    Make a scatter of circles plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circle in data scale (ie. in data unit)
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or
        RGBA sequence because that is indistinguishable from an array of
        values to be colormapped.  `c` can be a 2-D array in which the
        rows are RGB or RGBA, however.
    ax : Axes object, optional, default: None
        Parent axes of the plot. It uses gca() if not specified.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.  (Note if you pass a `norm` instance, your
        settings for `vmin` and `vmax` will be ignored.)

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Other parameters
    ----------------
    kwargs : `~matplotlib.collections.Collection` properties
        eg. alpha, edgecolors, facecolors, linewidths, linestyles, norm, cmap

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection
    import pylab as plt
    #import matplotlib.colors as colors

    if ax is None:
        ax = plt.gca()    

    if isinstance(c,basestring):
        color = c     # ie. use colors.colorConverter.to_rgba_array(c)
    else:
        color = None  # use cmap, norm after collection is created
    kwargs.update(color=color)

    if isinstance(x, (int, long, float)):
        patches = [Circle((x, y), s),]
    elif isinstance(s, (int, long, float)):
        patches = [Circle((x_,y_), s) for x_,y_ in zip(x,y)]
    else:
        patches = [Circle((x_,y_), s_) for x_,y_,s_ in zip(x,y,s)]
    collection = PatchCollection(patches, **kwargs)

    if color is None:
        collection.set_array(np.asarray(c))
        if vmin is not None or vmax is not None:
            collection.set_clim(vmin, vmax)

    ax.add_collection(collection)
    return collection







# Plots 2D dolfin solution with triangulation 
def Plot2DSolution(mesh,f,fileName="test.png",dpi=300,clim=False):
    
    import matplotlib.tri as tri
    n = mesh.num_vertices()
    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates().reshape((n, d))
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)

    plt.figure()
    plt.axis('equal')
    zfaces = np.asarray([f(cell.midpoint()) for cell in cells(mesh)])
    plt.tripcolor(triangulation, facecolors=zfaces, edgecolors='k')
    plt.colorbar()
    if clim!=False:
      plt.clim(clim)
    plt.tight_layout()
    plt.gcf().savefig(fileName,dpi=dpi)

# tri - Delaunay triangulation 
def PlotStats(ps,rs,length,tri=False):
    #plt.subplot(2,1,1)
    plt.scatter(ps[:,0],ps[:,1])    
    circles(ps[:,0],ps[:,1],rs)
    plt.axes().set_aspect('equal')
    plt.xlim([0,length])
    plt.ylim([0,length])
    plt.locator_params(nbins=3)

    #plt.figure()
    #mu.circles(ps[:,0],ps[:,1],rs)
    #axes().set_aspect('equal')
    #plt.triplot(points[:,0], points[:,1], tri.simplices.copy(),color='r')

    #plt.subplot(2,1,2)
    plt.figure()
    n,bins,patches=plt.hist(rs, normed=1, facecolor='green', alpha=0.75)
    plt.xlabel("R$_{Gyr}$ [$\AA$]")
    plt.ylabel("Prob")
    font = {'weight' : 'bold',
            'size'   : 22}
    import matplotlib
    matplotlib.rc('font', **font)
    plt.locator_params(nbins=3)
    plt.tight_layout()
    #plt.gcf().savefig(outputDir+"radGyrHisto.png",dpi=300)


    import cytoutil
    plt.figure()
    closest = cytoutil.closestDists(ps)
    plt.xlabel("Distance of closest approach [$\AA$]")
    plt.ylabel("Prob")
    plt.hist(closest,normed=1,bins=50)
    plt.tight_layout()
    plt.locator_params(nbins=3)
    #plt.gcf().savefig(outputDir+"pairwiseDistHisto.png",dpi=300)

    
