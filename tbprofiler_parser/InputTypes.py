import os
import argparse

class InputTypes:
  """
  Class for handling the input types for the TBProfiler_Parser tool.
  """
  
  def is_json_valid(filename):
    """
    Checks if the input JSON is valid
    """
    if not os.path.exists(filename) and filename != "-":
      raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
    return filename
  
  def is_bam_valid(filename):
    """
    Checks if the input BAM is valid and if there is an associated BAI
    """
    if not os.path.exists(filename) and (filename != "-" and not filename.endswith(".bam")):
      raise argparse.ArgumentTypeError("{0} cannot be accessed or is missing the BAM extension".format(filename))
    bai_filename = filename[:-4] + ".bai"
    if not os.path.exists(bai_filename) and (bai_filename != "-"):
      raise argparse.ArgumentTypeError("Cannot find the BAI file for this BAM")
    return filename