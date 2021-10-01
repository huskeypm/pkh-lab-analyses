"""For the collection step"""

def get_by_tag(self,newcol="other",attrpath="other_input",tagcol="tag",dfpath="df"):
  df=self.get_nested(dfpath)
  other_data=self.get_nested(attrpath)
  def do_calc(row):
    tag=row[tagcol]
    return other_data[tag]
  df[newcol]=df.apply(do_calc,axis=1)
  return

#List of functions to be bound as methods
request_methods=[get_by_tag]
