"""Support for building FEniCS weak forms"""

#Standard library
from collections import OrderedDict

#This package
from ..requesthandler import yaml_manager

class EquationTerm(object):
  """To keep track of additional information about each UFL form added up into the equation to solve.
  
  Attributes:
  
    - name = a string uniquely identifying the equation term
    - form = the UFL form object for the term"""

  def __init__(self,name,form,**kwargs):
    self.name=name
    self.form=form
    for k,v in kwargs.items():
      setattr(self,k,v)
  
  def asdict(self):
    d={}
    for k,v in self.__dict__.items():
      if k == 'form':
        d[k]=repr(v)
      else:
        d[k]=v
    return d

  # @classmethod
  # def to_yaml(cls, representer, node):
  #   return representer.represent_mapping("!"+cls.__name__,node.asdict())

class EquationTermDict(OrderedDict):
  """An ordered dictionary of equation terms.
  
  Attributes:
  
    - termclass = the class of which all equationterms are instances"""

  def __init__(self,termclass=EquationTerm,*args,**kwargs):
    #Initialization from base class
    super(EquationTermDict, self).__init__(*args,**kwargs)
    #Store the term class
    self.termclass=termclass

  def __getstate__(self):
    """Used for pickling, and converting to yaml"""
    #This method is only needed because I can't get the EquationTerm class to do this itself.
    #That would be a better approach, if it worked.
    return {k:v.asdict() for k,v in self.items()}

  def add(self,*args,**kwargs):
    """Create a new term and add it to this dictionary.
    
    Arguments are the same as those to create a new EquationTerm"""
    term=self.termclass(*args,**kwargs)
    self[term.name]=term

  def selectterms(self,**kwargs):
    """Return an EquationTermDict of terms with the specified properties"""
    out=EquationTermDict(self.termclass)
    for name,term in self.items():
      include=True
      for k,v in kwargs.items():
        if getattr(term,k)!=v:
          include=False
          break
      if include:
        out[name]=term
    return out

  def sumterms(self,zeroval=None,**kwargs):
    """Return the sum of terms with the specified properties, as a UFL form
    
    Arguments:
    
      - zeroval = optional, UFL form to be returned if no terms match the specified properties
      - **kwargs = keyword arguments for selectterms"""
    terms=self.selectterms(**kwargs)
    if len(terms)>0:
      return sum([t.form for t in terms.values()])
    elif zeroval is not None:
      return zeroval
    else:
      #zeroval was not provided
      raise RuntimeError("Unable to return a proper sum of zero terms; Provide keyword argument zeroval to resolve.")

#Register for dumping to yaml
# yaml_manager.register_classes([EquationTerm, EquationTermDict])
yaml_manager.register_classes([EquationTermDict])
