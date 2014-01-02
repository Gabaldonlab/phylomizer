#!/usr/bin/python
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
from tree_reconstructor import *
from utils import *
import sys, getopt
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def main(argv):

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  #~ parameters = {
    #~ "inFile":  "",
    #~ "verbose":  1,
    #~ "replace":  False,
    #~ "numb_models": 2,
    #~ "outDirec": "./",
    #~ "phyml": "phyml",
    #~ "nj_parameters": " -d aa -b 0 -f e -v e -a e -o l -c 4 ",
    #~ "ml_parameters": " -d aa -b -2 -f e -v e -a e -o tlr -c 4 ",
    #~ "evol_models": "JTT WAG MtREV VT LG Blosum62 Dayhoff"
  #~ }
  parameters = {
    "inFile": ['file', ''],
    "config": ['file', ''],
    "replace": ['bool', False],
    "outDirec": ['directory', '.']
  }
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  try:  opts, args = getopt.getopt(argv, "f:d:c:r", ["file=", "dest=", \
    "replace", "config="])
  except getopt.GetoptError: sys.exit("ERROR: Check the input parameters")

  for opt, arg in opts:
    if   opt in ("-f", "--file"):     parameters["inFile"][1]   = str(arg)
    elif opt in ("-d", "--dest"):     parameters["outDirec"][1] = str(arg)
    elif opt in ("-c", "--config"):   parameters["config"][1]   = str(arg)
    elif opt in ("-r", "--replace"):  parameters["replace"][1]  = True
    #~ elif opt in ("-p", "--prgram"):  parameters["phyml"]       = str(arg)
    #~ elif opt in ("-v", "--verbose"): parameters["verbose"]     = int(arg)
    #~ elif opt in ("-n", "--number"):  parameters["numb_models"] = int(arg)
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
  #~ if not lookForFile(parameters["inFile"]):
    #~ sys.exit("ERROR: Input file not available")

  if int(parameters["numb_models"]) > len(parameters["evol_models"].split(" "))\
    or int(parameters["numb_models"]) < 1:
    sys.exit("ERROR: The number of models to be evaluated have to be defined "
             "between 1 and " + str(len(parameters["evol_models"].split(" "))))
#~
  #~ parameters["outDirec"] = os.path.abspath(parameters["outDirec"])
  #~ if not lookForDirectory(parameters["outDirec"]):
      #~ sys.exit("ERROR: It is impossible to create the output directory")
#~
  #~ if not lookForProgram(parameters["phyml"]):
    #~ sys.exit("ERROR: It is impossible to find PHYML")
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
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
