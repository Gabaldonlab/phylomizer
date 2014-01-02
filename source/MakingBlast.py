#!/usr/bin/python
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
import sys, os, getopt
from blaster import *
from utils import *
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
def main(argv):

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  parameters = {
    "inFile": ['file', ''],
    "config": ['file', ''],
    "replace": ['bool', False],
    "BlastDB": ['file', ''],
    "outDirec": ['directory', '.']
  }  
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  try: opts, args = getopt.getopt(argv, "f:d:b:c:r", ["file=", "dest=",\
                    "replace", "db=", "config="])
  except getopt.GetoptError: sys.exit("ERROR: Check the input parameters")

  for opt, arg in opts:
    if   opt in ("-f", "--file"):     parameters["inFile"][1]   = str(arg)
    elif opt in ("-d", "--dest"):     parameters["outDirec"][1] = str(arg)
    elif opt in ("-b", "--db"):       parameters["BlastDB"][1]  = str(arg)
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
  if not lookForFile(parameters["inFile"]):
    sys.exit("ERROR: Input file not available")

  if not lookForFile(parameters["BlastDB"]):
    sys.exit("ERROR: BlastDB file not available")

  parameters["outDirec"] = os.path.abspath(parameters["outDirec"])
  if not lookForDirectory(parameters["outDirec"]):
      sys.exit("ERROR: It is impossible to create the output directory")

  if not lookForProgram(parameters["blastpgp"]):
    sys.exit("ERROR: It is impossible to find PHYML")

  if float(parameters["coverage"]) > 1 or float(parameters["coverage"]) <= 0:
    sys.exit("ERROR: The coverage has to be defined between 0 and 1")

  if int(parameters["hits"]) < 1:
    sys.exit("ERROR: The Number of hits for the Blast Output has to be greater"\
             + "than 0")
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  blast(parameters)
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
