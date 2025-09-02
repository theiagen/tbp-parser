import os
import argparse

def is_json_valid(filename):
  """Checks if the input JSON is accessible
  
  Args:
    filename (String): The name of file to check
        
  Returns:
    String: The name of the file if valid and accessible
  """
  if not os.path.exists(filename) and filename != "-":
    raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
  return filename

def is_bam_valid(filename):
  """Checks if the input BAM is valid and accessible and if there is an associated BAI
  
  Args:
    filename (String): The name of file to check
        
  Returns:
    String: The name of the file if valid and accessible
  """
  if not os.path.exists(filename) and (filename != "-" and not filename.endswith(".bam")):
    raise argparse.ArgumentTypeError("{0} cannot be accessed or is missing the BAM extension".format(filename))
  bai_filename = filename + ".bai"
  if not os.path.exists(bai_filename) and (bai_filename != "-"):
    raise argparse.ArgumentTypeError("Cannot find the BAI file for this BAM")
  return filename

def is_bed_valid(filename):
  """Checks if the coverage regions BED file is accessible
  
  Args:
    filename (String): The name of file to check
      
  Returns:
    String: The name of the file if valid and accessible
  """
  # # these two lines are needed when I run this locally
  # scripts_dir = os.path.dirname(os.path.realpath(__file__))
  # bed_file = os.path.join(scripts_dir, filename)
  if not os.path.exists(filename) and filename != "-":
    raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
  return filename

def is_config_valid(filename):
  """Checks if the configuration file is accessible

  Args:
    filename (String): The name of the file to check
      
  Returns:
    String: the name of the file if accessible
  """
  if filename != "" and not os.path.exists(filename) and filename != "-":
    raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
  return filename

    