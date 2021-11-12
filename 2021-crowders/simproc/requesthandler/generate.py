"""Support for requests that create their own children.

NOTE: this module has generally been made obsolete by the joblist module,
which works very differently but achieves the same purpose
in a much more understandable way."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from copy import deepcopy
import itertools

#Site packages

#This package
from . import yaml_manager
from .request import Request, get_request_class
from . import nested
from . import customization

_ParametricRequestListRequest_props_schema_yaml="""#ParametricRequestListRequest
request_type:
  anyOf:
    - {type: string}
    - {type: request}
constants: {type: object}
variations:
  type: object
  additionalProperties: {type: array}
other_parents:
  type: object
  additionalProperties: {type: request}
_children: {type: array}"""


class ParametricRequestListRequest(customization.CustomizableRequest):
  """General and base class for requests that parametrically generate child requests

  The process is:
  
    1) From the parent request attributes, construct the input fields for each child request
    2) From the input fields for each child request, compute the values of the keyword arguments to pass to its constructor
  
  Step 1 involves combining constant fields, fields that vary parametrically, and fields containing other requests.
  
  Step 2 involves calling the `get_child_kwargs` method.
  If the method returns None, then no request is generated for that combination of input fields.
  This can be used to filter out invalid combinations of parameters.

  You can subclass this by overriding the `get_child_kwargs` method.
  You can also get the same effect using the customization interface
  to load a replacement `get_child_kwargs` method from a module.
  Otherwise, the input fields are assumed to be identical to the child request keyword arguments.
  
  The fields passed to `get_child_kwargs` will also include 'index',
  which is the zero-based index number of the child in the sequence.
  
  User-defined attributes:
  
    - request_type = EITHER
          a type instance to use for the child requests
        OR
          a string containing the name of a request type registered with yaml_manager
    - constants = dictionary of variables that will be the same in each dictionary used to calculate request kwargs:
        {fieldname: value, ...}
    - variations = dictionary of fields that will vary over all possible combinations:
        {fieldname: [value, ...], ...}
        For example, {a: [1,2,3], b: [1,2,3]} will generate 9 combinations of a and b: (a=1,b=1), (a=1,b=2), (a=1,b=3), (a=2,b=1), ...
    - other_parents = dictionary of fields that will come from the children of other parent requests:
       {fieldname: parent_request, ...}
       For example, {other_request: parent} will create a field called `other_request`, which for each child request
       will contain one (immediate) child of the request called `parent`.
       Also, {request1: parent1, request2: parent2}, where parent1 and parent2 each have 3 child requests,
       will generate 9 combinations of request1 and request2.
       The iteration is only over the immediate children of the parent requests.
  
  Calculated Attributes:
  
    - _children = sequence of child requests"""
  _validation_schema=customization.CustomizableRequest.update_schema(_ParametricRequestListRequest_props_schema_yaml)
  _validation_schema.required=['request_type']
  _child_seq_attrs=['_children']
  _self_task=False #This request generates doit tasks from its children, not itself
  def get_child_kwargs(self,index=None,**fields):
    """Compute a keyword arguments dictionary from the input dictionary"""
    outkwargs={}
    outkwargs.update(fields)
    return outkwargs
  def __init__(self,**kwargs):
    #Initialization from base class
    super(ParametricRequestListRequest, self).__init__(**kwargs)
    #Get a handle to the actual child request class
    childclass=get_request_class(self.request_type)
    #Fields
    const_fields=getattr(self,'constants',{})
    variations=getattr(self,'variations',{})
    variation_fieldnames=list(variations.keys())
    variation_fieldnames.sort() #The order the fields were provided in might not have been preserved, so alphabetize.
    variation_fieldvalues=[variations[k] for k in variation_fieldnames]
    other_parents=getattr(self,'other_parents',{})
    op_fieldnames=tuple(other_parents.keys())
    op_c_generators=[p.all_children() for p in other_parents.values()] #The generator of all children for each other parent
    op_c_iterator=itertools.product(*op_c_generators) #Iterator over combinations of the other parents children
    #Loop through
    self._children=[]
    index=0
    for opc_tup in op_c_iterator:
      opc_fields=dict(zip(op_fieldnames,opc_tup))
      variation_iterator=itertools.product(*variation_fieldvalues)
      for variation_values in variation_iterator:
        variation_fields=dict(zip(variation_fieldnames,variation_values))
        #Put the fields together
        fields={'index':index}
        fields.update(deepcopy(const_fields))
        fields.update(deepcopy(variation_fields))
        fields.update(opc_fields)
        #Obtain arguments for child request constructor
        child_kwargs=self.get_child_kwargs(**fields)
        #Create child and add to list
        if child_kwargs is not None:
          ch_req=childclass(**child_kwargs)
          self._children.append(ch_req)
          index += 1
    return

_GeneratedVariationsRequest_props_schema_yaml="""#GeneratedVariationsRequest
request_type:
  anyOf:
    - {type: string}
    - {type: request}
template: {type: object}
variations:
  type: array
  items:
    type: object
    properties:
      attrloc: {type: attrpath}
      values: {type: array}
    required: [attrloc, values]
    additionalProperties: False
other_parents:
  type: object
  additionalProperties: {type: request}
parents_mapping:
  type: array
  items:
    type: object
    properties:
      parent_loc: {type: attrpath}
      template_loc: {type: attrpath}
    required: [parent_loc, template_loc]
    additionalProperties: False
_children: {type: array}"""

##TODO: this borrows lots of code from the above. Can we consolidate somehow?
##TODO: it would be nice if this could allow for "sets" of properties that vary together, instead of always taking the product
##Maybe this could be done by adding a new optional attribute: a list of "groups" of parameters.
##This defaults to an empty list, meaning to take the product of all parameters.
##Every parameter in each group must have the same number of values.
##The product of parameters includes this factor only once, not once for each parameter.
##This should only apply to items in "variations", not to "other_parents".
##Unless you could figure out a way to list parameters by destination:
##Then you could combine items from the attrlocs in variations with template_locs in parents_mapping.
##But that sounds really complicated.
class GeneratedVariationsRequest(customization.CustomizableRequest):
  """Generate requests based on a template with listed variations
  
  User-defined attributes:

    - request_type = EITHER
          a type instance to use for the child requests
        OR
          a string containing the name of a request type registered with yaml_manager
    - template = template for keyword arguments dictionary to the child request type
    - variations = sequence of changes to make to the template, each item a dictionary:
        - attrloc: attribute path to the value to vary in the child request keyword arguments dictionary
        - values: list of values to place in the specified attribute path
    - other_parents = dictionary of other parent requests:
       {parent_name: parent_request, ...}
    - parents_mapping = sequence indicating how to map values from other_parents into template locations, each item a dictionary:
        - parent_loc: attribute path to location within the parent
        - template_loc: attribute path to location within the template
  
  Note that the template may contain properties overridden by the variations.
  In that case, the parameters in the template itself are not used.
  
  """
  _validation_schema=customization.CustomizableRequest.update_schema(_GeneratedVariationsRequest_props_schema_yaml)
  _validation_schema.required=['request_type','template']
  _child_seq_attrs=['_children']
  _self_task=False #This request generates doit tasks from its children, not itself
  def get_child_kwargs(self,index=None,**fields):
    """Compute a keyword arguments dictionary from the input dictionary"""
    outkwargs={}
    outkwargs.update(fields)
    outkwargs['name']+=".%04d"%index
    return outkwargs
  def __init__(self,**kwargs):
    #Initialization from base class
    super(GeneratedVariationsRequest, self).__init__(**kwargs)
    #Get a handle to the actual child request class
    childclass=get_request_class(self.request_type)
    #Get the listed variations
    variations=getattr(self,'variations',[])
    variation_attrlocs=[v['attrloc'] for v in variations]
    variation_fieldvalues=[v['values'] for v in variations]
    #Load the other requests
    other_parents=getattr(self,'other_parents',{})
    op_fieldnames=tuple(other_parents.keys())
    op_c_generators=[p.all_children() for p in other_parents.values()] #The generator of all children for each other parent
    op_c_iterator=itertools.product(*op_c_generators) #Iterator over combinations of the other parents children
    #Parents mapping
    parents_mapping=getattr(self,'parents_mapping',[])
    #Loop through
    self._children=[]
    index=0
    for opc_tup in op_c_iterator:
      opc_dict=dict(zip(op_fieldnames,opc_tup))
      variation_iterator=itertools.product(*variation_fieldvalues)
      for variation_values in variation_iterator:
        #Put the fields together
        fields=deepcopy(self.template)
        fields['_related']=deepcopy(opc_dict)
        #Information from other parents
        for pdict in parents_mapping:
          opvalue=nested.get_nested(opc_dict,pdict['parent_loc'])
          nested.set_nested(fields,pdict['template_loc'],opvalue)
        #Information from variations
        for idx, vvalue in enumerate(variation_values):
          nested.set_nested(fields,variation_attrlocs[idx],vvalue)
        #Obtain arguments for child request constructor
        child_kwargs=self.get_child_kwargs(index=index,**fields)
        #Create child and add to list
        if child_kwargs is not None:
          ch_req=childclass(**child_kwargs)
          self._children.append(ch_req)
          index += 1
    return

#Register for loading from yaml
yaml_manager.register_classes([ParametricRequestListRequest, GeneratedVariationsRequest])
