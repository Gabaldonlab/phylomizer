import sys, os, subprocess as sp
from string import strip
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

def PhylogeneticTrees(parameters, approach, evolutionary_models):
  ''' A phylogenetic tree is reconstructed using the parameters defined in the
      function input. Once the different files have been generated, the function
      moves them following a preestablish scheme
  '''

  ## Set how output file names will be
  common = os.path.join(parameters["outdirec"], parameters["prefix"])
  ## Determine which program will be used to reconstruct trees
  program = "phyml" if "phyml" in parameters else "fasttree" if "fasttree" in \
    parameters else "raxml" if "raxml" in parameters else None

  ## If there is not any program for which there is already a wrapper, return
  if not program:
    return

  ## If the evolutionary model list is not appropiately formated, do it
  if type(evolutionary_models) == str:
    evolutionary_models = map(strip, evolutionary_models.split())

  results = []
  ## Explore the different models and get the likelihood of each of them
  for model in evolutionary_models:
    outFile = ("%s.%s.tree.%s.%s.nw") % (common, program, approach, model)
    statsFile = ("%s.%s.tree.%s.%s.st") % (common, program, approach, model)

    ## Call to the appropiate wrapper depending on the selected
    if program == "phyml":
      cmd = ("%s -i %s %s -m %s") % (parameters["phyml"], parameters["inFile"],\
        parameters[approach], model)

      lk = wrapperPhyML(cmd, parameters["inFile"], statsFile, outFile,
        parameters["log"], parameters["replace"])

      results.append((lk, model))

  results = [(pair[1], pair[0]) for pair in sorted(results, reverse = True)]

  outFile = ("%s.%s.tree.rank.%s") % (common, program, approach)
  if parameters["replace"] or not utils.lookForFile(outFile):
    ranking = "\n".join(["\t".join(map(str, pair)) for pair in results])
    print >> open(outFile, "w"), ranking

  return [pair[0] for pair in results]
