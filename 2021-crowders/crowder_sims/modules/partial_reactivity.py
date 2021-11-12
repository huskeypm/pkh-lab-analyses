"The class for a partially reactive boundary surface"

import fenics as fem
import numpy as np

class SphericalSectorBoundary(fem.SubDomain):
  def __init__(self,cen,rad,areafrac,refpt,tol):
    super(SphericalSectorBoundary,self).__init__()
    self.cen=np.array(cen)
    self.rad=rad
    assert len(cen) == len(refpt), "Inconsistent dimensions: center=%s, refpt=%s"%(str(cen),str(refpt))
    if len(cen)==3:
      angle_cos=1-2*areafrac #For 3D
    elif len(cen)==2:
      angle_cos=np.cos(np.pi*areafrac) #For 2D
    else:
      assert False, "Must have either 2 or 3 dimensions, got %d"%(len(cen))
    self.angle_cos=angle_cos
    self.axial_vec=np.array(refpt)-self.cen
    self.axial_mag=np.sqrt(sum(self.axial_vec**2))
    self.tol=tol
  def inside(self,x,on_boundary):
    #Check that point is on a boundary
    if not on_boundary:
      return False
    else:
      #Check that point is on the surface of this sphere
      this_r=np.sqrt(sum((x-self.cen)**2))
      on_sphere = abs(this_r-self.rad)<self.tol
      if not on_sphere:
        return False
      else:
        #Check that point is within the region
        rvec=x-self.cen
        rmag=np.sqrt(sum(rvec**2))
        this_cos=np.dot(rvec,self.axial_vec)/rmag/self.axial_mag
        return this_cos >= self.angle_cos