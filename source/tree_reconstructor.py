import sys, os, subprocess as sp
from operator import itemgetter
import utils

''' Module which implements the functionality for reconstructing phylogenetic
    trees. It contains wrappers to three different programs: PhyML, RAxML &
    FastTree
'''

def wrapperFastTree(cmd, statsFile, logFile):
  ''' Wrapper to call FastTree and generate phylogenetic trees either NJ or ML
      ones
  '''

  ## Run Fast Tree
  ## Depending on whether there is a file descriptor open for capturing program
  ## output or not, we will use two different calls
  try:
    if logFile:
      process = sp.Popen(cmd, shell = True, stderr = logFile, stdout = logFile)
    else:
      process = sp.Popen(cmd, shell = True)
  except OSError, e:
    sys.exit(("ERROR: Execution failed for FastTree\n###\tCMD\t%s\t###\tReport:"
      + "\t%s") % (cmd, str(e)))

  if process.wait() != 0:
    sys.exit("ERROR: Execution failed for FastTree")

  ## Get likelihood values for the current run
  logLK = None
  for line in open(statsFile, "rU"):
    f = map(strip, line.split("\t"))
    try:
      value = float(f[2])
    except:
      continue
    logLK = value if not logLK or value < logLK else logLK
  ## Return the likelihood value for the current tree
  return logLK

def wrapperPhyML(cmd, inFile, statsFile, outFile, logFile):
  ''' Wrapper to call PhyML, generate a phylogenetic tree according to the input
      parameters and rename the output files accordingly
  '''

  ## Run PhyML
  ## Depending on whether there is a file descriptor open for capturing program
  ## output or not, we will use two different calls
  try:
    if logFile:
      process = sp.Popen(cmd, shell = True, stderr = logFile, stdout = logFile)
    else:
      process = sp.Popen(cmd, shell = True)
  except OSError, e:
    sys.exit(("ERROR: Execution failed for PhyML\n###\tCMD\t%s\t###\tReport:"
      + "\t%s") % (cmd, str(e)))

  if process.wait() != 0:
    sys.exit("ERROR: Execution failed for PhyML")

  ## Rename output files
  try:
    sp.call(("mv %s_phyml_tree.txt %s") % (inFile, outFile), shell = True)
    sp.call(("mv %s_phyml_stats.txt %s") % (inFile, statsFile), shell = True)
  except OSError:
    sys.exit("ERROR: Impossible to rename PhyML output files")

  logLK = None
  ## Get likelihood values for current run
  try:
    pipe = sp.Popen(("grep 'Log-likelihood:' %s") % (statsFile), shell = True, \
      stdout = sp.PIPE)
    logLK = float(map(strip, pipe.stdout.readlines()[0].split())[2])
  except OSError, e:
    sys.exit("ERROR: Impossible to get LK values for PhyML run")
  return logLK

def wrapperRAxML(cmd, suffix, statsFile, outFile, outdirec, logFile):
  ''' Wrapper to call RAxML, generate a phylogenetic tree according to the input
      parameters, rename the output files accordingly and remove those which are
      not needed
  '''

  ## Run RAxML
  ## Depending on whether there is a file descriptor open for capturing program
  ## output or not, we will use two different calls
  try:
    process = sp.Popen(cmd, shell = True, stderr = logFile, stdout = logFile)
  except OSError, e:
    sys.exit(("ERROR: Execution failed for RAxML\n###\tCMD\t%s\t###\tReport:"
      + "\t%s") % (cmd, str(e)))

  if process.wait() != 0:
    sys.exit("ERROR: Execution failed for RAxML")

  ## Rename output files
  try:
    sp.call(("mv RAxML_bestTree%s %s") % (suffix, outFile), shell = True)
    sp.call(("mv RAxML_info%s %s") % (suffix, statsFile), shell = True)
  except OSError:
    sys.exit("ERROR: Impossible to rename RAxML output files")

  logLK = None
  ## Parse output to get the likelihood of the best tree
  for line in open(statsFile, "rU"):
    if line.lower().startswith("final") and line.lower().find("score") != -1:
      logLK = float(map(strip, line.split())[-1])

  oFile = open(statsFile, "a+")
  for oth_file in utils.listDirectory(outdirec, suffix):
    fileName = os.path.split(oth_file)[1]
    hz_line = "#" * (len(fileName) + 4)
    print >> oFile, ("%s\n%s\n%s") % (hz_line, fileName, hz_line)
    print >> oFile, ("%s") % ("".join(open(oth_file, "rU").readlines()))
    sp.call(("rm -f %s") % (oth_file), shell = True)
  oFile.close()

  return logLK

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
