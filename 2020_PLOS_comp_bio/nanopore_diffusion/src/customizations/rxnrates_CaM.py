#Reaction rate functions

def initialize_module(**kwargs):
  for k,v in kwargs.items():
    globals()[k]=v

def rate_forward(self,cca,ccacam):
    return cca*(BT-ccacam)

def rate_backward(self,cca,ccacam):
    return ccacam