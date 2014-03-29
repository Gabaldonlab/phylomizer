import os
# import re
import sys
import datetime
import subprocess as sp

# from time import sleep
# from hashlib import md5
from string import strip
# from socket import getfqdn
# from getpass import getuser
# from operator import itemgetter
from module_utils import lookForDirectory, lookForFile, splitSequence, \
  format_time

''' Module which implements the functionality for reconstructing phylogenetic
    trees. It contains wrappers to three different programs: PhyML, RAxML &
    FastTree
'''

def phylogenetic_trees(parameters):
  ''' Phylogenetic trees are reconstructed according to the input parameters.
      Once the different files have been generated, the function moves those
      files into a pre-established filename schema
  '''

  ## Get output folder/generic filename
  oFile = os.path.join(parameters["out_directory"], parameters["prefix"])

  current_directory = os.getcwd()
  ## Change current directory to the output folder. Any temporary file will be
  ## generated therefore in this folder
  os.chdir(parameters["out_directory"])

  ## Set output filename and log file
  if parameters["replace"] and parameters["step"] == 0:
    logFile = open(oFile + ".log", "w")
  else:
    logFile = open(oFile + ".log", "a+")

  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> logFile, ("###\n###\tSTEP\tPhylogenetic Tree Reconstruction\tSTART\t"
    + "%s\n###") % (date)
  logFile.flush()

  ## Get which program will be used to reconstruct phylogenetic trees. Check
  ## such program is listed among the available binaries
  if not "tree" in parameters:
    sys.exit("ERROR: Check your configuration file. There is no definition for "
      + "the Phylogenetic TREE reconstruction step")

  program = parameters["tree"][0]
  if not program in parameters:
    sys.exit(("ERROR: Selected program '%s' is not available accordding to the "
      "the configuration file") % (program))

  if not "evol_models" in parameters:
    sys.exit("ERROR: Check your configuration file. There is no definition for "
      + "the <evol_models> parameter")

  ## If the evolutionary model list is not appropiately formated, do it
  if isinstance(parameters["evol_models"], basestring):
    parameters["evol_models"] = map(strip, parameters["evol_models"].split())

  ## Check if <numb_models parameters is defined and how many models are
  ## requested to be evaluated
  if not "numb_models" in parameters or parameters["numb_models"].lower() \
    == "all":
    parameters["numb_models"] = len(parameters["evol_models"])
  parameters["numb_models"] = int(parameters["numb_models"])

  if not parameters["numb_models"] in range(1, len(parameters["evol_models"])):
    sys.exit(("ERROR: Check how many evolutionary models has been asked to re"
      + "construct '%d'") % (parameters["numb_models"]))

  ## Check which approaches should be used for the phylogenetic reconstruction
  ## and whether there are specific program's parameters for them
  if not "tree_approach" in parameters:
    parameters["tree_approach"] = ["ml"]


  for approanch in tree_approach:
    if not





#~ def PhylogeneticTrees(parameters, approach, evolutionary_models):
#~
#~
#~
  #~ results = []
  #~ ## Explore the different models and get the likelihood of each of them
  #~ for model in evolutionary_models:
    #~ outFile = ("%s.%s.tree.%s.%s.nw") % (common, program, approach, model)
    #~ statsFile = ("%s.%s.tree.%s.%s.st") % (common, program, approach, model)
#~
    #~ ## Call to the appropiate wrapper depending on the selected
    #~ if program == "phyml":
      #~ cmd = ("%s -i %s %s -m %s") % (parameters["phyml"], parameters["inFile"],\
        #~ parameters[approach], model)
#~
      #~ lk = wrapperPhyML(cmd, parameters["inFile"], statsFile, outFile,
        #~ parameters["log"], parameters["replace"])
#~
      #~ results.append((lk, model))
#~
  #~ results = [(pair[1], pair[0]) for pair in sorted(results, reverse = True)]
#~
  #~ outFile = ("%s.%s.tree.rank.%s") % (common, program, approach)
  #~ if parameters["replace"] or not utils.lookForFile(outFile):
    #~ ranking = "\n".join(["\t".join(map(str, pair)) for pair in results])
    #~ print >> open(outFile, "w"), ranking
#~
  #~ return [pair[0] for pair in results]

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

def wrapperPhyML(cmd, inFile, statsFile, outFile, logFile, replace):
  ''' Wrapper to call PhyML, generate a phylogenetic tree according to the input
      parameters and rename the output files accordingly
  '''

  ## Run PhyML depending on if the output files have been previously generated
  ## and the user is asking to overwrite them
  if not utils.lookForFile(statsFile) or not utils.lookForFile(outFile) or \
    replace:
    try:
      process = sp.Popen(cmd, shell = True, stderr = logFile, stdout = logFile)
    except OSError, e:
      sys.exit(("ERROR: Execution failed for PhyML\n###\tCMD\t%s\t###\tReport:"
        + "\t%s") % (cmd, str(e)))

    ## Press a key until any situation to continue the execution
    process.communicate("Y\n")
    if process.wait() != 0:
      sys.exit("ERROR: Execution failed for PhyML")

    ## Rename output files
    try:
      sp.call(("mv %s_phyml_tree.txt %s") % (inFile, outFile), shell = True)
      sp.call(("mv %s_phyml_stats.txt %s") % (inFile, statsFile), shell = True)
    except OSError:
      sys.exit("ERROR: Impossible to rename PhyML output files")

  logLK = None
  ## Get likelihood value for the execution
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
