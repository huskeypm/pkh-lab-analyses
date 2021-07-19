"""Generate random inclusions for 2D meshes"""

#Standard library
from collections import namedtuple

#Site Packages
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from scipy.spatial.distance import pdist

#Local modules

Rect=namedtuple('Rect',['xmin','ymin','xmax','ymax'])

def reseed(randseed=0):
  """Reseed the random number generator"""
  np.random.seed(randseed)

def calc_minseps(rvals,eps=0.1):
  """Calculate array of minimum separation values
  
  The minimum separation distance of two circles is the sum of their radii,
  and the tolernace value.
  Circles closer together than this are considered to overlap,
  or be too close to properly construct the mesh.
  
  Arguments:
  
    - rvals = 1D array of circle radii
    - eps = tolerance distance.
  
  Returns:
  
    - minseps = 1D array of minimum separation distance, with one entry for each pair of circles
    
  The size of minseps, for N circles, is N*(N-1)/2"""
  npts=len(rvals)
  minseps=[]
  for i in range(npts):
    for j in range(i+1,npts):
      minseps.append(rvals[i]+rvals[j]+eps)
  minseps=np.array(minseps)
  return minseps

def is_valid(centers,rvals,box,eps=0.1):
  """Check if the set of circles in a bounding box is valid.
  
  Arguments:
  
    - centers = Nx2 array of circle centers: (x,y) for each circle
    - rvals = 1D array of circle radii, of length N
    - box = bounding rectangle, as a Rect instance
    - eps = tolerance distance
  
  Returns:
  
    - True if no circles overlap other circles or the box boundary,
        False otherwise"""
  #First, check if any circles extend beyond the box boundary
  for i,pt in enumerate(centers):
    r=rvals[i]
    if pt[0]-r-eps<box.xmin or pt[0]+r+eps>box.xmax or pt[1]-r-eps<box.ymin or pt[1]+r+eps>box.ymax:
      return False
  #Check for collisions between circles
  minseps=calc_minseps(rvals,eps)
  collisions=[pr[0]<=pr[1] for pr in zip(pdist(centers),minseps)]
  return not any(collisions)

def mk_circle(xrange,yrange,rrange):
  """Generate one random circle
  
  Arguments:
  
    - xrange = tuple of parameters for generating random x-value (uniform distribution: limits)
    - yrange = tuple of parameters for generating random y-value ( ' ' )
    - rrange = tuple of parameters for generating random radii ( ' ' )
  
  Returns:
  
    - centers = Nx2 array of circle center coordinates (x,y)
    - rvals = 1D array of circle radii, of length N"""
  xmin,xmax=xrange
  ymin,ymax=yrange
  rmin,rmax=rrange
  xvals=xmin+(xmax-xmin)*np.random.random(1)
  yvals=ymin+(ymax-ymin)*np.random.random(1)
  rvals=rmin+(rmax-rmin)*np.random.random(1)
  centers=np.vstack((xvals,yvals)).T
  return (centers,rvals)

def combine(c1,r1,c2,r2):
  """Combine two groups of center and radius values into a single group"""
  if c1.shape[0]>0:
    centers=np.vstack((c1,c2))
    rvals=np.hstack((r1,r2))
  else:
    centers=c2
    rvals=r2
  return (centers,rvals)

def mk_test(centers,rvals,xrange,yrange,rrange):
  """Add a single trial circle to the list of existing circles"""
  cadd,radd=mk_circle(xrange,yrange,rrange)
  outc,outr=combine(centers,rvals,cadd,radd)
  return outc,outr

def fixed_to_arrays(fixed):
  """Return arrays for a set of fixed (that is, nonrandom) circles"""
  xvals=np.array([p[0] for p in fixed])
  yvals=np.array([p[1] for p in fixed])
  rvals=np.array([p[2] for p in fixed])
  return xvals,yvals,rvals

def random_circles(fixed,nrandom,xrange,yrange,rrange,box,maxtries=1e4,eps=0.1):
  """Create a set of random circles that stay inside the bounding box, and don't overlap each other.
  
  Arguments:
  
    - fixed = sequence of fixed circles, each circle as a tuple (x,y,r)
    - nrandom = the number of random circles to generate
    - xrange, yrange, rrange = tuples defining random variable distributions for x,y, and r, respectively.
    - box = Rect defining the bounding box within which the circular inclusions must remain
    - maxtries = optional, maximum number of times to try adding a single point before giving up
    - eps = tolerance distance for overlap: circles must have at least this much separation
  
  Returns:
  
    - centers = Nx2 array of circle center coordinates (x,y)
    - rvals = 1D array of circle radii, of length N"""
  xvals,yvals,rvals=fixed_to_arrays(fixed)
  centers=np.vstack((xvals,yvals)).T
  assert is_valid(centers,rvals,box,eps), "Fixed circles are not valid even without random ones."
  for ptnum in range(nrandom):
    trial=0
    testc,testr=mk_test(centers,rvals,xrange,yrange,rrange)
    while (not is_valid(testc,testr,box,eps)) and trial<maxtries:
      trial += 1
      testc,testr=mk_test(centers,rvals,xrange,yrange,rrange)
    # print(trial)
    assert is_valid(testc,testr,box,eps), "Failed to find valid result."
    centers=testc
    rvals=testr
  return centers,rvals

def nonrandom_circles(fixed,ncircles,rval,box,maxtries=1e4,eps=0.1):
  """Create a set of nonrandom circles that stay inside the bounding box, and don't overlap each other.
  
  Arguments:
  
    - fixed = sequence of fixed circles, each circle as a tuple (x,y,r)
    - ncircles = the APPROXIMATE number of circles to generate; the algorithm just covers the area
    - rval = circle radius, as float
    - box = Rect defining the bounding box within which the circular inclusions must remain
    - eps = tolerance distance for overlap: circles must have at least this much separation
  
  Returns:
  
    - centers = Nx2 array of circle center coordinates (x,y)
    - rvals = 1D array of circle radii, of length N"""
  #Process fixed circles
  xvals,yvals,rvals=fixed_to_arrays(fixed)
  centers=np.vstack((xvals,yvals)).T
  assert is_valid(centers,rvals,box,eps), "Fixed circles are not valid even without adding others."
  #For zero circles, that's all
  if ncircles <= 0:
    return centers, rvals
  else:
    #Calculate circle spacing
    xspan=box.xmax-box.xmin
    yspan=box.ymax-box.ymin
    full_A=xspan*yspan
    fixed_A_list=[np.pi*r**2 for xc,yc,r in fixed]
    fixed_A=sum(fixed_A_list)
    net_A=full_A-fixed_A
    approx_step=np.sqrt(net_A/ncircles)
    # approx_step=np.sqrt(full_A/ncircles)
    Xsteps=max(1,int(np.floor(xspan/approx_step)))
    Ysteps=max(1,int(np.floor(yspan/approx_step)))
    xsize=xspan/Xsteps
    ysize=yspan/Ysteps
    #Loop over the area to be covered
    placed=[]
    ynow=box.ymin+ysize/2
    while ynow < box.ymax:
      xnow = box.xmin+xsize/2
      while xnow < box.xmax:
        #Place circle if it doesn't overlap any fixed circles
        seps=[np.sqrt((xnow-fx)**2+(ynow-fy)**2)-rval-fr-eps for fx,fy,fr in fixed]
        if len(seps)==0 or min(seps)>0:
          placed.append((xnow,ynow,rval))
        #Next space
        xnow += xsize
      ynow += ysize
    #Combine with fixed circles
    addx,addy,addr=fixed_to_arrays(placed)
    addc=np.vstack((addx,addy)).T
    newcenters,newrvals = combine(centers,rvals,addc,addr)
    return newcenters,newrvals

def display(ax,centers,rvals,box):
  boxpatch=patches.Rectangle((box.xmin,box.ymin),box.xmax-box.xmin,box.ymax-box.ymin,fill=False,ls='-')
  o=ax.add_patch(boxpatch)
  for i,cen in enumerate(centers):
    c=patches.Circle(cen,rvals[i])
    o=ax.add_patch(c)
  ax.autoscale_view()
  ax.set_aspect('equal')
  
