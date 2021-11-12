"Add spherical inclusions on a regular lattice inside a cylinder"

import numpy as np
import numpy.linalg as la

import simproc.requesthandler.logging as logging
logger=logging.getLogger(__name__)

NUDGE=1e-9

def check_collisions(trial_abc,trial_xyz,inc_eff,cyl_eff,eff_inclusions):
  """Confirm that the trial point is not too close to the cylinder wall or any static inclusions.
  This includes checking that the point is not beyond the cylinder wall.

  Arguments:

    - trial_abc = point to test in abc coordinate system, as tuple
    - trial_xyz = point to test in xyz coordinate system, as numpy array
    - inc_eff = inclusion radius plus half buffer, as float
    - cyl_eff = cylinder radius minus half buffer, as float
    - eff_inclusions = 2D numpy array, each row an inclusion [x,y,z,r]
  
  Returns:

    - valid = boolean, True if point is OK, False if not"""
  #First, confirm that the point is close enough to the cylinder center axis
  polar_rho = np.sqrt(trial_abc[0]**2 + trial_abc[1]**2)
  valid = (polar_rho + inc_eff) < cyl_eff
  #Next, check against all the static inclusions, if any
  if valid and eff_inclusions.shape[0]>0:
    distvecs = eff_inclusions[:,:3] - trial_xyz
    dists = la.norm(distvecs,axis=1)
    limits = eff_inclusions[:,3] + inc_eff
    bool_arr = dists > limits
    valid = all(bool_arr)
    # logger.debug("Completed collision check.",distvecs=distvecs,dists=dists,limits=limits,bool_arr=bool_arr,valid=valid)
  return valid

def gen_perturbation(oldpt,stdev):
  """Create a new trial position for a given point from a random perturbation

  Arguments:

    - oldpt = coordinates tuple, (x,y,z) as floats
    - stdev = standard deviation magnitude, as float

  Returns:

    - newpt = coordinates tuple, (x,y,z) as float
  """
  #Select random direction and size
  #Direction is random in spherical coordinates
  phi = np.random.random()*np.pi #Angle in interval [0, pi)
  sintheta = np.random.random()*2.0 - 1.0 #sin in interval [-1, 1)
  theta = np.pi/2.0 + np.arcsin(sintheta) #theta in interval [0, pi), not uniformly distributed
  step = np.random.normal(0,stdev) #Signed step-size
  #Convert to Cartesian coordinates
  dx = step*sintheta*np.cos(phi)
  dy = step*sintheta*np.sin(phi)
  dz = step*np.sin(theta)
  #Add deltas to original coordinates
  oldx,oldy,oldz=oldpt
  newpt = (oldx + dx, oldy + dy, oldz + dz)
  #Done
  return newpt

def check_perturbation(newpt,buffer,statics,lattice,inc_idx,cylrad,endpoints,hhat,H):
  """Check the trial position of a point against all collisions

  Arguments:

    - newpt = trial position
    - buffer = minimum surface separation distance, as float
    - statics = list of static inclusions
    - lattice = list of lattice inclusions
    - inc_idx = index of the perturbed inclusion in the lattice inclusions list
    - cylrad = cylinder radius, as float
    - endpoints = cylinder axis endpoints, as 2D array
    - hhat = numpy 1D array, coordinate of unit vector along cylinder axis
    - H = distance between endpoints, as float

  Returns:

    - valid = boolean, True if change is OK, False if not
  """
  nx,ny,nz=newpt
  myrad=lattice[inc_idx][4]
  #Check against cylinder boundary
  mvec = np.array(newpt) - endpoints[0]
  pmag = np.dot(mvec,hhat)
  qvec = mvec - pmag*hhat
  valid = (cylrad - la.norm(qvec) - myrad - buffer) >= 0
  #Confirm the point is not beyond the ends of the cylinder
  if valid:
    valid = myrad + buffer <= pmag
  if valid:
    valid = pmag <= H - myrad - buffer
  #Check against the other inclusions
  if valid:
    total_checklist = statics + [t for idx,t in enumerate(lattice) if idx != inc_idx]
    for incid,cx,cy,cz,incrad in total_checklist:
      cendist = la.norm([nx-cx,ny-cy,nz-cz])
      valid = (cendist - buffer - incrad - myrad) >= 0
      if not valid:
        break
  #Done
  return valid

def generate_lattice_inclusions_in_cylinder(idstart,incrad,cylrad,buffer,pt1,pt2,inclusions,idstep=100,perturbation_rounds=0,perturbation_stdev=0):
  """Create points the points

  Arguments:

    - idstart = ID to use fo the first inclusion, as integer
    - incrad = inclusion radius, as float
    - cylrad = cylinder radius, as float
    - buffer = minimum surface separation distance, as float
    - pt1 = coordinates of cylinder axis start point
    - pt2 = coordinates of cylinder axis end point
    - inclusions = non-lattice inclusions to avoid collision with, as list of inclusions,
      each inclusion a tuple (ID, center x, center y, center z, radius)
    - idstep = number to increment the inclusion id by, optional, defaults to 100
    - perturbation_rounds = Number of times to try perturbing the position of each inclusion, as integer
    - perturbation_stdev = Standard deviation of the perturbation magnitude, as float
  
  Returns:

    - lattice_inclusions = list of inclusions, each inclusion a tuple as for ``inclusions``."""
  #Get the cylinder axis vector, h hat, and its magnitude, H
  endpoints = np.array([pt1,pt2])
  hvec = endpoints[1] - endpoints[0]
  H = la.norm(hvec)
  hhat = hvec/H
  #Get effective dimensions after accounting for separation buffer distances
  hafbuff = buffer/2.0
  cyl_eff = cylrad - hafbuff
  inc_eff = incrad + hafbuff
  H_eff = H - buffer
  #Apply buffer separation to static inclusions by increasing their radius
  eff_inclusions=[]
  for statinc in inclusions:
    newval=list(statinc[1:4])
    newval.append(statinc[4] + hafbuff)
    eff_inclusions.append(newval)
  eff_inclusions=np.array(eff_inclusions)
  #Generate unit vectors a and b that are mutually orthogonal and orthogonal to h hat
  v1 = np.array([1.0,0.0,0.0]) if hhat[0]<0.8 else np.array[0.0,0.0,1.0]
  ahat = v1 - np.dot(v1,hhat)*hhat
  ahat /= la.norm(ahat)
  bhat = np.cross(hhat,ahat)
  #Get the number of lattice points in each dimension
  Nh = int(np.floor(H_eff/2/inc_eff))
  Na = Nb = int(np.floor((cyl_eff+inc_eff)/2/inc_eff - NUDGE))
  #Scale the unit vectors into appropriate step sizes
  ascale = 2 * cyl_eff / (2*Na-1)
  bscale = 2 * cyl_eff / (2*Nb-1)
  cscale = H_eff / Nh
  #Loop over lattice points
  incid=idstart
  lattice_inclusions=[]
  for i in range(0,2*Na-1):
    Ai = (i - Na + 1) * ascale
    for j in range(0,2*Nb-1):
      Bj = (j - Nb + 1) * bscale
      for k in range(0,Nh):
        Ck = (k + 0.5) * cscale + buffer/2
        trial_abc = (Ai,Bj,Ck)
        trial_pt = endpoints[0] + Ai*ahat + Bj*bhat + Ck*hhat
        if check_collisions(trial_abc,trial_pt,inc_eff,cyl_eff,eff_inclusions):
          this_inc=(incid,)+tuple(trial_pt)+(incrad,)
          lattice_inclusions.append(this_inc)
          incid+=idstep
  #Perturbations
  #Perform the prescibed number of rounds
  for pround in range(perturbation_rounds):
    #Try to peturb each lattice position
    for inc_idx,this_inc in enumerate(lattice_inclusions):
      incid,oldx,oldy,oldz,radval = this_inc
      #Generate peturbed position
      oldpt = (oldx,oldy,oldz)
      newpt = gen_perturbation(oldpt, perturbation_stdev)
      #Check for collision
      accepted = check_perturbation(newpt,buffer,inclusions,lattice_inclusions,inc_idx,cylrad,endpoints,hhat,H)
      #If trial position is accepted, update
      if accepted:
        newtup = (incid,)+newpt+(radval,)
        lattice_inclusions[inc_idx]=newtup
  #Done
  return lattice_inclusions