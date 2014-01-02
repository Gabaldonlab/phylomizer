#!/usr/bin/python
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
from tree_reconstructor import *
from aligner import *
from utils import *
import sys, getopt
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
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
  try: opts, args = getopt.getopt(argv, "f:d:c:r", ["file=", "dest=", \
    "replace", "config="])
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
  if int(parameters["numb_models"]) > len(parameters["evol_models"].split(" "))\
    or int(parameters["numb_models"]) < 1:
    sys.exit("ERROR: The number of models to be evaluated have to be defined "
             "between 1 and " + str(len(parameters["evol_models"].split(" "))))
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  parameters["inFile"] = AlignerPipeline(parameters)
  # parameters["inFile"] = SimpleAlignerPipeline(parameters)  

  PhylogeneticTrees(parameters["phyml"], parameters["nj_parameters"], "nj",
    parameters["inFile"], parameters["outDirec"], parameters["evol_models"],
    parameters["verbose"], parameters["replace"])

  rank = SortingLKs(parameters["inFile"], parameters["outDirec"], "nj",
    parameters["evol_models"], parameters["verbose"])

  rank = ' ' . join(rank[:int(parameters["numb_models"])])

  PhylogeneticTrees(parameters["phyml"], parameters["ml_parameters"], "ml",
    parameters["inFile"], parameters["outDirec"], rank, parameters["verbose"],
    parameters["replace"])

  SortingLKs(parameters["inFile"], parameters["outDirec"], "ml", rank,
    parameters["verbose"])
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

## ***** ***** ***** ***** ***** ***** ***** ***** *****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
