"""Geneate gmsh .geo file(s) from standard geometric layout."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#This package
from ..requesthandler import filepath
from ..requesthandler import yaml_manager
from ..requesthandler.request import Request
from ..requesthandler import schema
from ..requesthandler import locators
from ..requesthandler import logging

logger=logging.getLogger(__name__)

#Site packages
from jinja2 import Environment, FileSystemLoader

#The package directory
modpath=filepath.Path(__file__)
pkgdir=modpath.parent
#The main template file
geotemplate=pkgdir / "common.geo.jinja2"

_GeometryDefinition_props_schema_yaml="""#GeometryDefinition
dimensions:
  type: integer
  minimum: 2
  maximum: 3
tmplfile: {type: pathlike}
tmplvars: {type: array}
outvars: {type: array}
ptdict: {type: object}
geomtable: {type: object}
surfloops: {type: object}
nonplanar: {type: array}
"""

class GeometryDefinition(schema.SelfValidating):
  """Geometry definition parameters
  
  Attributes:

    - dimensions = number of spatial dimensions (i.e. 2 for 2D mesh, 3 for 3D mesh)
    - tmplfile = geometry template file
    - tmplvars = mesh parameter variables needed by the geometry template file
    - outvars = list of gmsh variables to be output to the mesh metadata file (note: template variables are not always provided directly to gmsh)
    - ptdict = dictionary of points and their corresponding mesh density parameter name
    - geomtable = mapping of surfaces to sequence points
    - surfloops = mapping of surface loops to sequence of surfaces
    - nonplanar = list of surfaces that are not planar surfaces"""
  _validation_schema=schema.SelfValidating.update_schema(_GeometryDefinition_props_schema_yaml)
  _validation_schema.required=['dimensions', 'tmplfile', 'tmplvars', 'outvars', 'ptdict', 'geomtable', 'surfloops', 'nonplanar']

schema.extra_types_dict['GeometryDefinition']=(GeometryDefinition,)

#From mapping of surfaces to points, generate:
# - mapping of loops to line and circle names
# - mapping of line and circle names to list of points
def add_entity(tup, tdict, looplist, nameprefix):
  """Create lines/circles needed for line loops, unless they already exist.

  Arguments:

    - tup = line/circle tuple of point numbers (as integers) (lines have 2 points, circles have 3)
    - tdict = dictionary storing the line/circle tuples
    - looplist = list of lines/circles in this line loop
    - nameprefix = eg 'C' for circle, 'L' for line

  No return value.

  Side effects:

    - tdict and looplist are modified in place, to include the generated (or found) lines/circles"""
  #Reversed tuple
  rtup = tuple(reversed(tup))
  found=False
  #Search through all lines/circles already created, for one with same set of points, in forward or reverse order
  for n, pts in tdict.items():
    if pts==tup:
      found=True
      looplist.append(n)
      break
    elif pts==rtup:
      found=True
      looplist.append('-'+n)
      break
  #Create new line/circle if not found
  if not found:
    nametmpl='_'.join('%d' for i in range(len(tup)))
    name=nameprefix+nametmpl%tup
    tdict[name]=tup
    looplist.append(name)
  return

_BuildGeomRequest_props_schema_yaml="""#BuildGeomRequest
geomdef:
  anyOf:
    - {type: GeometryDefinition}
    - {type: stored}
parameters: {type: object}
geofile: {type: pathlike}
searchpaths: {type: array}
"""

class BuildGeomRequest(Request):
  """Generate a gmsh .geo file from a specified geometric layout
  
  User-defined attributes:
  
    - geomdef = instance of GeometryDefinition
    - parameters = dictionary of parameter input values
    - geofile = path to output .geo file
    - searchpaths = optional list of other paths to search for supporting template files (i.e. jinja template extension files)
  """
  _self_task=True
  _outputfile_attrs=['geofile']
  _config_attrs=['geomdef','parameters']
  _validation_schema=Request.update_schema(_BuildGeomRequest_props_schema_yaml)
  _validation_schema.required=['name','geomdef','parameters','geofile']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(BuildGeomRequest, self).__init__(**kwargs)
    #Load geometry definition if Stored
    self.geomdef=self.get_stored(self.geomdef)
    #List the geometry definition template file as an input file
    self._more_inputfiles=[geotemplate,self.geomdef.tmplfile]
  # def additional_validation(self,**kwargs):
  #   """Perform additional validation of the object data, beyond just the schema check
    
  #   Check that the provided parameters are appropriate for the template."""
  #   expected=kwargs['geomdef']['tmplvars'] #Doesn't work!!!
  #   provided=kwargs['parameters'].keys()
  #   missing=[var for var in expected if not var in provided]
  #   extra=[var for var in provided if not var in expected]
  #   errlist=[]
  #   if len(missing)>0:
  #     errstr="  - Missing required parameters for the template:\n"
  #     for var in missing:
  #       errstr += "    - %s\n"%str(var)
  #     errlist.append(errstr)
  #   if len(extra)>0:
  #     errstr="  - Extra parameters provided not known by the template:\n"
  #     for var in extra:
  #       errstr += "    - %s\n"%str(var)
  #     errlist.append(errstr)
  #   return errlist
  def prepare_template_input(self):
    """Prepare the input dictionary for a template.

    No arguments.

    Returns:

      - t_input = the input dictionary for the template"""
      
    #Check that the necessary variables are defined
    missing=[var for var in self.geomdef.tmplvars if not var in self.parameters.keys()]
    extra=[var for var in self.parameters.keys() if not var in self.geomdef.tmplvars]
    assert len(missing)==0 and len(extra)==0, "Missing or extra template variables. missing=%s extra=%s"%(str(missing),str(extra))
    
    #Put geometric and mesh refinement parameters into template input
    t_input=dict(self.parameters)
    
    #Put dimensions into mesh
    t_input['dimensions']=self.geomdef.dimensions

    #Complete listing of all output variables for mesh metadata
    t_input['metadata_vars']=self.geomdef.outvars

    #Dictionary of points
    t_input['ptstrs']=dict([(str(x),y) for x,y in self.geomdef.ptdict.items()])

    #Create dictionaries of lines, circles, and line loops
    loops={}
    lines={}
    circles={}
    for surfnum, pttup in self.geomdef.geomtable.items():
      loops[surfnum]=[]
      startpt=pttup[0]
      indx=1
      while indx < len(pttup):
        if pttup[indx]=='center':
          indx += 2
          ctup=(startpt,pttup[indx-1],pttup[indx])
          add_entity(ctup,circles,loops[surfnum],'C')
        else:
          ltup=(startpt,pttup[indx])
          add_entity(ltup,lines,loops[surfnum],'L')
        #Next point
        startpt=pttup[indx]
        indx += 1

    #Generate physical lines, if needed
    if self.geomdef.dimensions == 2:
      t_input['phys_lines']={}
      for lname in lines.keys():
        pts=lname[1:].split('_')
        pts.sort()
        lnum=''.join(pts)
        t_input['phys_lines'][lnum]=lname

    #Provide mappings to template
    linemap=dict([(n,', '.join(['p%d'%p for p in pts])) for n,pts in lines.items()])
    t_input['lines']=linemap
    circmap=dict([(n,', '.join(['p%d'%p for p in pts])) for n, pts in circles.items()])
    t_input['circles']=circmap
    loopmap=dict([(n,', '.join([x for x in ents])) for n,ents in loops.items()])
    t_input['loops']=loopmap

    #Apply surface types
    # surftypes=dict([(x,'Ruled Surface' if x in self.geomdef.nonplanar else 'Plane Surface') for x in self.geomdef.geomtable.keys()]) #for older version of gmsh
    surftypes=dict([(x,'Surface' if x in self.geomdef.nonplanar else 'Plane Surface') for x in self.geomdef.geomtable.keys()])
    t_input['surftypes']=surftypes

    #Dictionary of surface loops (and volumes)
    t_input['surfloops']=dict([(n, ', '.join(['%d'%x for x in surfs])) for n,surfs in self.geomdef.surfloops.items()])

    return t_input

  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()

    #Get the input dictionary for the template
    t_input = self.prepare_template_input()

    #Construct list of all paths to search for templates
    tmplfile=filepath.Path(self.render(self.geomdef.tmplfile))
    searchpaths=[str(pkgdir),str(tmplfile.parent)]+[str(p) for p in getattr(self,'searchpaths',[])]

    #Load template
    env=Environment(loader=FileSystemLoader(searchpaths),extensions=['jinja2.ext.do'],trim_blocks=True,keep_trailing_newline=True)
    tmpl=env.get_template(tmplfile.filename)

    #Render template
    outdat = tmpl.render(t_input)

    #Output result
    with open(self.renderstr(self.geofile),'w') as fh:
      fh.write(outdat)

    #Done
    return

#Register for loading from yaml
yaml_manager.register_classes([GeometryDefinition, BuildGeomRequest])
