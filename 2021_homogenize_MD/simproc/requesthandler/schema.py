"""Define class that can validate itself against a json schema"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages
import jsonschema

#This package
from . import filepath
from . import yaml_manager
from . import locators
from . import nested

#Validation partial setup (some setup must wait for Request class to be defined)
ValidatorClass = jsonschema.Draft4Validator
loc_class_tup=(locators.DataFile,locators.NameDelegator,locators.Delegator)
#jsonschema 2.6
extra_types_dict={'path':filepath.Path,
                  'locator':loc_class_tup,
                  'pathlike':(str,filepath.Path)+loc_class_tup, #This isn't the same thing as "pathlike" in python.org documentation
                  'array':(list,tuple),
                  'attrpath':(str,list,tuple),
                  'stored':nested.Stored}

def validation_error_string(err):
  "Return a string explaining the given validation error"
  #Get the basic message
  s=err.message
  #Provide a path if necessary
  if len(err.path)>0:
    s="%s: %s"%('.'.join([str(itm) for itm in err.path]),s)
  return s

#Basic schema for self-validation
#Note that validation is not applied to class attributes,
#but some of these could be instance attributes in some cases.
#Note also that properties inherited from nested.WithNested must also be included in the schema.
_SelfValidating_schema_yaml="""#SelfValidating
type: object
properties:
  name:
    anyOf:
      - {type: 'null'}
      - {type: string}
  store_globally: {type: boolean}
required: []
additionalProperties: False
"""
_SelfValidating_schema=yaml_manager.readstring(_SelfValidating_schema_yaml)


class SelfValidating(nested.WithNested):
  """A class that can validate itself against a json schema
  
  Derived classes can set their attribute ``_validation_schema`` to an
  instance of nested.WithNested containing the complete json schema to be used
  in validating the dictionary form of an instance.
  
  This class itself provides a basic validation schema, which derived classes may update
  using the ``update_schema`` method.
  
  Derived classes may specify required attributes by setting the property
  ``_validation_schema.required`` to a list of required attribute names, as strings.
  
  By default, only listed properties are allowed.
  To change this behavior, set ``_validation_schema.additionalProperties`` to ``True``."""
  _validation_schema=nested.WithNested(**_SelfValidating_schema)
  def __init__(self,**kwargs):
    #Validate kwargs
    if hasattr(self,'_validation_schema'):
      self.validate_kwargs(**kwargs)
    #Load the attributes specified
    super(SelfValidating, self).__init__(**kwargs)
  @classmethod
  def update_schema(cls,yaml_str,dpath="properties",update=True):
    """Create a property schema for a class

    The intention is for subclasses to call this from their parent class.
    However, once an initial schema has been created, a class may call it on itself.

    Arguments:

      - yaml_str = string containing yaml defining updates to the (parent) class validation schema
      - dpath = optional nested attribute path to the portion of the parent schema to be updated

          Defaults to ``properties``, as that is most frequently what is changed

      - update = optional boolean, True (default) to update, False to overwrite
    
    Returns an instance of nested.WithNested containing the updated schema, ready for validation."""
    sub_schema=yaml_manager.readstring(yaml_str)
    orig_schema=cls._validation_schema
    new_schema=orig_schema.get_copy()
    if update:
      new_schema.update_nested(dpath,sub_schema)
    else:
      new_schema.set_nested(dpath,sub_schema)
    return new_schema
  def additional_validation(self,**kwargs):
    """Perform additional validation of the object data, beyond just the schema check

    The method provided in this base class is meant to be overridden by subclasses that need it.

    Arguments:

      - \*\*kwargs = arguments dictionary of attributes

    Returns:

      - list of error strings (including their indentation)"""
    return []
  def validate_kwargs(self,**kwargs):
    if hasattr(self,'_validation_schema'):
      validator=ValidatorClass(self._validation_schema,types=extra_types_dict)
      err_iter=validator.iter_errors(kwargs)
      errlist=["  - %s"%validation_error_string(err) for err in err_iter]
      errlist+=self.additional_validation(**kwargs)
      if len(errlist)>0:
        #Found errors: raise exception listing them all
        errlist.sort()
        errstr="Errors found in %s.\n"%type(self).__name__
        keylist=list(kwargs.keys())
        keylist.sort()
        errstr+='Received arguments:\n'
        errstr+='\n'.join(['  - %s: %s'%(k,kwargs[k]) for k in keylist])
        errstr+='\nErrors:\n'
        errstr+='\n'.join(errlist)
        raise Exception(errstr)
  def validate(self):
    d=self.to_dict(False)
    self.validate_kwargs(**d)
  