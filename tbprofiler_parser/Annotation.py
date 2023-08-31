""" 
This class represents the annotation field from TBProfiler.
"""

class Annotation:
  def __init__(self, annotation, drug):
    self.annotation = annotation
    self.drug = drug
  
  """
  This function ranks the WHO annotation based on resistance,
  with 4 being the most resistant category and 1 the least.
  """
  def rank_annotation(self, annotation):
    if annotation == "Assoc w R":
      return 4
    elif annotation == "Assoc w R - interim":
      return 3
    elif annotation == "Uncertain significance":
      return 2
    else:
      return 1
  
  
  """
  This function converts the WHO annotation and the target drug
  into returns the LIMS' report file appropriate annotation.
  """
  def annotation_to_LIMS(self, annotation, drug):
    if annotation == "Assoc w R":
      return "Mutation(s) associated with resistance to {} detected".format(drug)
    elif (annotation == "Assoc w R - interim") or (annotation == "Uncertain significance"):
      return "The detected mutation(s) have uncertain significance. Resistance to {} cannot be ruled out".format(drug)
    # "Not assoc w R" and "Not assoc w R - Interim" and anything else
    else: 
      return "No mutations associated with resistance to {} detected".format(drug)