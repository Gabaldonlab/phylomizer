#!/usr/bin/python
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
from aligner_new import *
from utils import *
import sys, getopt

def main(argv):

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  parameters = {
    "inFile": ['file', ''],
    "config": ['file', ''],
    "replace": ['bool', False],
    "outDirec": ['directory', '.']
  }
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  try: opts, args = getopt.getopt(argv, "f:d:c:r", ["file=", "dest=", "config="\
    "replace"])
  except getopt.GetoptError: sys.exit("ERROR: Check the input parameters")

  for opt, arg in opts:
    if   opt in ("-f", "--file"):     parameters["inFile"][1]   = str(arg)
    elif opt in ("-d", "--dest"):     parameters["outDirec"][1] = str(arg)
    elif opt in ("-c", "--config"):   parameters["config"][1]   = str(arg)
    elif opt in ("-r", "--replace"):  parameters["replace"][1]  = True
    else: sys.exit("ERROR: Parameter not Available")
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if not lookForFile(parameters["config"][1]):
    sys.exit("ERROR: Config file not available")
  else: readConfig(parameters["config"][1], parameters)
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  print "Pipeline Configuration\n---"
  for key in parameters:
    print key.center(20) + "\t" + str(parameters[key])
  print "---"
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  AlignerPipeline(parameters)
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

## ***** ***** ***** ***** ***** ***** ***** ***** *****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
