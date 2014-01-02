import sys, os, subprocess as sp
from operator import itemgetter
from utils import *

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def PhylogeneticTrees(bin, pars, approach, inFile, oFolder, EvolModels, verb, repl):
  '''
  A phylogenetic tree is reconstructed using the parameters defined in the
  function input. Once the different files have been generated, the function
  moves them following a preestablish scheme
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  iFolder, iFile = os.path.split(inFile)
  oFile = os.path.join(oFolder, iFile.split(".")[0])
  lFile, msg = open(oFile + ".log", "a+"), ""
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  for model in EvolModels.split(" "):
    outFile = oFile + ".tree." + approach + "." + model
    if not lookForFile(outFile + ".nw") or not lookForFile(outFile + ".st") or \
       repl:
      cmd  = bin + " -i " + inFile + " "+ pars + " -m " + model
      cmd2 = "mv " + inFile  + "_phyml_tree.txt  " + outFile + ".nw"
      cmd3 = "mv " + inFile  + "_phyml_stats.txt " + outFile + ".st"

      try: process = sp.Popen(cmd, shell = True, stderr = lFile, stdout = lFile,
            stdin = sp.PIPE)
      except OSError, e:   sys.exit("ERROR: Execution failed: " + str(e))
      process.communicate("Y\n")
      if process.wait() != 0: sys.exit("ERROR: Execution failed: PhyML")  
      
      try: sp.call(cmd2, shell = True)
      except OSError: sys.exit("ERROR: Impossible to rename the output files.")

      try: sp.call(cmd3, shell = True)
      except OSError: sys.exit("ERROR: Impossible to rename the output files.")

    else:
      msg += "The files " + outFile + ".[nw, st] already exist in the folder\n"
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  lFile.close()

  if verb > 0: print msg,
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def SortingLKs(inFile, oFolder, approach, EvolModels, verbose):
  '''
  The likelihood for different evolutionary model are evaluated and sorted by
  this function. An output file is generated with a pair <model, lk> for each
  model evaluate as well a sorted list is returned.
  '''

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  ranking = []
  iFolder, iFile = os.path.split(inFile)
  oFile = os.path.join(oFolder, iFile.split(".")[0])
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  for model in EvolModels.split(" "):
    cmd = "grep 'Log-likelihood:' " + oFile + ".tree." + approach + "." + \
          model + ".st"
    try: pipe = sp.Popen(cmd, shell = True, stdout = sp.PIPE).stdout
    except OSError, e: sys.exit("ERROR: Execution failed: " + str(e))
    ranking.append((model, float(pipe.readlines()[0].split()[2])))
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if lookForFile(oFile + ".tree.rank." + approach) and verbose > 0:
    print "Replacing the file " + oFile + ".tree.rank." + approach + \
          " with updated information"

  f1le = open(oFile + ".tree.rank." + approach, "w")
  ranking = sorted(ranking, key = itemgetter(1), reverse = True)
  for item in ranking: print >> f1le, "%-10s\t%f" % (item[0], item[1])
  f1le.close()
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  return [ rank[0] for rank in ranking ]
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
