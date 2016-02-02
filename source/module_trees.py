"""
  phylomizer - automated phylogenetic reconstruction pipeline - it resembles the
  steps followed by a phylogenetist to build a gene family tree with error-control
  of every step

  Copyright (C) 2014 - Salvador Capella-Gutierrez, Toni Gabaldon

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
import datetime
import tempfile
import subprocess as sp

from socket import getfqdn
from random import randint
from operator import itemgetter
from string import strip, lower

from module_utils import format_time, listDirectory
from module_utils import lookForDirectory, lookForFile, splitSequence

from module_alignments import convertInputFile_Format, getFileFormat
from module_alignments import check_count_sequences, replaceRareAminoAcids

## To guarantee consistency across the pipeline, the minimum sequenes number to
## perform any analysis is defined just in one module
from module_alignments import  min_seqs_analysis

''' Module which implements the functionality for reconstructing phylogenetic
    trees. It contains wrappers to three different programs: PhyML, RAxML &
    FastTree
'''

## Exit code meanings
##  80: Not enough sequences
exit_codes = {
  "phyml":          96,    ##  96: Problems associated to PhyML execution
  "codonphyml":     97,    ##  97: Problems associated to CodonPhyML execution
  "raxml":          98,    ##  98: Problems associated to RAxML execution
  "fasttree":       99,    ##  99: Problems associated to FastTree execution

  "generic":        95,    ##  95: Program not supported
}

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
  open_mode = "w" if parameters["replace"] and parameters["step"] == 0 else "a+"
  logFile = open(oFile + ".log", open_mode)

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

  prog = parameters["tree"][0]
  if not prog in parameters:
    sys.exit(("ERROR: Selected program '%s' is not available accordding to the "
      "the configuration file") % (prog))

  ## Get binary as well as any default parameters for the selected program
  binary = parameters[prog]
  key = ("%s_params") % (prog)
  progr_params = parameters[key] if key in parameters else ""

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

  if not parameters["numb_models"] in range(1,len(parameters["evol_models"])+1):
    sys.exit(("ERROR: Check how many evolutionary models has been asked to re"
      + "construct '%d'") % (parameters["numb_models"]))

  ## Check whether "readAl" is available or not. It is useful for sequences
  ## manipulation independently of the input format.
  if not "readal" in parameters:
    sys.exit("ERROR: Check your CONFIG file. 'readAl' is not available")

  ## Create a temporary FASTA file which will be used to detect the sequence
  ## number on the input alignment and the presence of rare amino-acids
  TEMPFILE = tempfile.NamedTemporaryFile()
  convertInputFile_Format("readal", parameters["readal"], parameters["in_file"],
    TEMPFILE.name, "fasta", logFile, parameters["replace"])
  TEMPFILE.flush()

  numSeqs, selenocys, pyrrolys = check_count_sequences(TEMPFILE.name)

  ## Set the minimum number of sequences required to reconstruct an alignment
  min_seqs = int(parameters["min_seqs"] if "min_seqs" in parameters else \
    min_seqs_analysis)
  
  ## Finish when there are not enough sequences to make an alignment
  if numSeqs < min_seqs:
    print >> logFile, ("### INFO: It is necessary, at least, %d sequences to "
      + "to reconstruct an alignment (%d)") % (min_seqs, numSeqs)
    sys.exit(80)

  ## Check which approaches should be used for the phylogenetic reconstruction
  ## and whether there are specific program's parameters for them
  if not "tree_approach" in parameters:
    parameters["tree_approach"] = ["ml"]

  ## Remove potential duplicates and lowercase all approaches for the tree
  ## reconstruction
  parameters["tree_approach"] = set(map(lower, parameters["tree_approach"]))

  ## We will first loot for Neighbour Joining tree reconstruction, then for
  ## Maximum likelihood and then for any other approach defined in the config
  ## file
  tree_approaches = []
  if "nj" in parameters["tree_approach"]:
    tree_approaches.append("nj")
  if "ml" in parameters["tree_approach"]:
    tree_approaches.append("ml")
  others = parameters["tree_approach"] - set(["nj", "ml"])
  if others != set():
    tree_approaches += sorted(others)

  ## When using RAxML, it may crash when Selenocysteines or Pyrrolysines are
  ## present in the input alignment
  if prog in ["raxml"]:
    ## If Selenocysteines or Pyrrolysines are present, substitute them by "X"
    if selenocys or pyrrolys:
      out_file = ("%s.no_rare_aa") % (parameters["in_file"])

      if replaceRareAminoAcids(TEMPFILE.name, out_file, parameters["replace"],
        logFile, "U:X O:X"):
        parameters["replace"] = True
      parameters["in_file"] = out_file
    TEMPFILE.close()

  ## When using FastTree force the conversion of input alignment to FASTA format
  ## since it may crash reading standard interleave PHYLIP format files
  if prog in ["fasttree"]:

    in_file_format, aligned = getFileFormat("readal", parameters["readal"], \
      parameters["in_file"], logFile)

    if in_file_format != "fasta":
      out_file = ("%s.fa") % (parameters["in_file"])
      if (convertInputFile_Format("readal", parameters["readal"], \
        parameters["in_file"], out_file, "fasta", logFile,
        parameters["replace"])):
        parameters["replace"] = True
      parameters["in_file"] = out_file

  replace = parameters["replace"]
  selected_models = parameters["evol_models"]
  ## Reconstruct trees for each approach considering evolutionary models order
  ## according their likelihood values
  for approach in tree_approaches:

    ## Save results - we will use such data for selecting the best -if required-
    ## models fitting to the input data
    results = {}

    ## Format the choosen program's parameters according to the default ones and
    ## the specific ones for the current approach
    params = ("%s ") % (progr_params)
    params += parameters[approach] if approach in parameters else ""

    for model in selected_models:
      out_file = ("%s.tree.%s.%s.%s.nw") % (oFile, prog, approach, model)
      stats_file = ("%s.tree.%s.%s.%s.st") % (oFile, prog, approach, model)

      if prog in ["phyml"]:
        exec_params = ("%s -m %s") % (params, model)

      ## Get additional model -if any- for codons
      elif prog in ["codonphyml"]:
        exec_params = ("%s -m %s") % (params, model)

        add_model = [p.split()[1] for p in map(strip, exec_params.split("-")) \
          if p.startswith("fmodel")]

        if len(add_model) == 1:
          add_model = add_model.pop()
          model = ("%s_%s") % (model, add_model)
          out_file = ("%s.tree.%s.%s.%s.nw") % (oFile, prog, approach, model)
          stats_file = ("%s.tree.%s.%s.%s.st") % (oFile, prog, approach, model)

      elif prog in ["fasttree"]:
        ## On FastTree is selected by default JTT model for AAs - so we don't
        ## set-up that model
        exec_params = ("%s -%s") % (params, model) if model.lower() != "jtt" \
          and model.lower() != "jc" else params
        model = model.upper()

      ## In the case of RAxML, we would concatenate the model to an specific
      ## input parameter
      elif prog in ["raxml"]:
        final_model = model
        ## It is possible to add some suffixes to the evolutionary models
        ## in RAxML - There is not better/easy way to code this option
        if "raxml_model_suffix" in parameters:
          final_model += parameters["raxml_model_suffix"]
        exec_params = " ".join([("-%s%s") %(p, final_model if p.startswith("m ")
          else "") for p in map(strip, params.split("-")) if p])

      ## Build the phylogenetic tree using any of the available methods and
      ## register if any downstream file should be redone.
      if perform_tree(prog, binary, exec_params, parameters["in_file"],
        out_file, stats_file, logFile, parameters["replace"]):
          replace = True

      ## Get the likelihood for each of the reconstructed models
      log_lk = get_likelihood(prog, stats_file)

      if not log_lk:
        print >> sys.stderr, ("ERROR: Impossible to the Log likelihood values "
          + "for '%s' model using this program '%s'") % (model, prog)
        sys.exit(exit_codes[prog])

      results.setdefault(model, log_lk)

    ## Get the models sorted by their likelihood values
    records = sorted(results.iteritems(), key = itemgetter(1), reverse = True)

    ## Set the filename which stores the ranking
    rank_file = ("%s.tree.%s.rank.%s") % (oFile, prog, approach)

    update = False
    ## Check the content of the rankings file - if any.
    ## Marked the file as updatable if there is any discrepancy
    if not replace and lookForFile(rank_file):

      old_content = "\n".join(["\t".join(map(strip, line.split("\t"))) for line
        in open(rank_file, "rU")])      

      newly_generated = "\n".join([("%s\t%s") % (r[0], r[1]) for r in records])
      
      ## Decide whether ranking file should be updated after comparing current
      ## content with newly generated content
      update = old_content != newly_generated

    ## If the file containing the ranking doesn't exist, generate it.
    ## Update the file content if the replace flag is set to true or the content
    ## has changed - since the phylogenetic tree reconstruction step is the most
    ## expensive one - in terms of time/memory consumption - we are not setting
    ## replace flag to True even when this file is generated/updated. On this
    ## way, we can take adventage of any tree generated in any downstream step. 
    if not lookForFile(rank_file) or replace or update:

      out_file = open(rank_file, "w")
      print >> out_file, "\n".join([("%s\t%s") % (r[0], r[1]) for r in records])
      out_file.close()

      ## We could set the replace flag to True. However, if any tree has been
      ## generated 'de novo' during this iteration, then the flag is already set
      ## to True. 
      #~ parameters["replace"] = True

    ## Select a given number of models for the next iteration - if any
    selected_models = [pair[0] for pair in records[:parameters["numb_models"]]]

    ## Remove the Codon Frequency model from potential new iterations
    if prog in ["codonphyml"] and add_model:
      selected_models = [m.replace("_"+ add_model, "") for m in selected_models
        if m.endswith(add_model)]

  final = datetime.datetime.now()
  date = final.strftime("%H:%M:%S %m/%d/%y")
  print >> logFile, ("###\n###\tSTEP\tPhylogenetic Tree Reconstruction\tEND\t"
    + "%s") % (date)
    
  ## We return a DELTA object comparing both timestamps
  total = format_time(final - start if start else 0)
  print >> logFile, ("###\tTOTAL Time\tPhylogenetic Tree Reconstruction\t%s"
    + "\n###") % (total)
  logFile.close()

  ## Clean-up log directory from undesirable lines
  try:
    sp.call(("sed -i '/^$/d' %s.log") % (oFile), shell = True)
    sp.call(("sed -i '/^M/d' %s.log") % (oFile), shell = True)
    sp.call(("sed -i '/\r/d' %s.log") % (oFile), shell = True)
  except OSError:
    print >> sys.stderr, ("ERROR: Impossible to clean-up '%s.log' log file") \
      % (oFile)

  ## Before returning to the main program, get back to the original working
  ## directory
  os.chdir(current_directory)

  return parameters

def perform_tree(label, binary, parameters, in_file, out_file, stats_file, \
  logFile, replace):

  '''
  Function to format the command-line of different phylogenetic tree reconstruc-
  tion programs and execute such command lines.
  '''

  ## Check whether the output file already exists. If it is not set to replace
  ## it, just return to the calling function
  if lookForFile(out_file) and not replace:
    return False

  if label in ["phyml", "codonphyml"]:
    cmd = ("%s -i %s %s") % (binary, in_file, parameters)

  elif label in ["fasttree"]:
    cmd = ("%s %s -log %s -out %s %s") % (binary, parameters, stats_file, \
      out_file, in_file)

  elif label in ["raxml"]:
    random_seed = randint(1, 10000)
    suffix = ("%s_%d") % (label, random_seed)

    cmd = ("%s -n %s -p %d -s %s %s") % (binary, suffix, random_seed, in_file, \
      parameters)

  else:
    sys.exit(exit_codes["generic"])

  ## Record the time and precise command-line
  name = getfqdn()
  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")

  print >> logFile, ("###\n###\t%s - Phylogenetic Trees\t") % (label.upper()),
  print >> logFile, ("%s\n###\t[%s]\tCommand-line\t%s\n###") % (date, name, cmd)
  logFile.flush()

  try:
    ## We add a small pipeline to avoid informatin written in the same line
    proc = sp.Popen(cmd, shell = True, stderr = logFile, stdout = logFile,
      stdin = sp.PIPE)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(exit_codes[label])
  proc.stdin.write("\n\nY\n")

  if proc.wait() != 0:
    print >> sys.stderr, ("ERROR: Execution failed: %s") % (label.upper())
    sys.exit(exit_codes[label])

  final = datetime.datetime.now()
  ## We return a DELTA object comparing both timestamps
  total = format_time(final - start if start else 0)
  print >> logFile, ("###\tTime\t%s\n###") % (total)
  logFile.flush()

  ## Process program's output and rename output files according to our own
  ## scheme
  if label in ["phyml", "codonphyml"]:

    ## Since resulting tree/stats file have slightly changed between version,
    ## we have to control for that.
    tree_file = ("%s_%s_tree.txt") % (in_file, label)
    sts_file = ("%s_%s_stats.txt") % (in_file, label)
    if not lookForFile(tree_file, attempts = 2):
      tree_file = ("%s_%s_tree") % (in_file, label)   
      sts_file = ("%s_%s_stats") % (in_file, label)   

    try:
      sp.call(("mv %s %s") % (tree_file, out_file), shell = True)
      sp.call(("mv %s %s") % (sts_file, stats_file), shell = True)
    except OSError:
      print >> sys.stderr, ("ERROR: Impossible to rename '%s' output files") \
        % (label.upper())
      sys.exit(exit_codes[label])

  elif label in ["raxml"]:
    try:
      sp.call(("mv RAxML_bestTree.%s %s") % (suffix, out_file), shell = True)
      sp.call(("mv RAxML_info.%s %s") % (suffix, stats_file), shell = True)
    except OSError:
      print >> sys.stderr, ("ERROR: Impossible to rename RAxML output files")
      sys.exit(exit_codes[label])

    oFile = open(stats_file, "a+")
    for oth_file in listDirectory(os.path.split(stats_file)[0], suffix):
      fileName = os.path.split(oth_file)[1]
      hz_line = "#" * (len(fileName) + 4)
      print >> oFile, ("%s\n%s\n%s") % (hz_line, fileName, hz_line)
      print >> oFile, ("%s") % ("".join(open(oth_file, "rU").readlines()))
      sp.call(("rm -f %s") % (oth_file), shell = True)
    oFile.close()

  return True

def get_likelihood(label, stats_file):

  ## Check whether the STATS file is available or not
  if not lookForFile(stats_file):
    return None

  logLK = None
  ## PHYML/CodonPhyML
  if label in ["phyml", "codonphyml"]:
    for line in open(stats_file, "rU"):
      if not line.startswith(". Log-likelihood"):
        continue
      logLK = float(map(strip, line.split())[2])
      break

  ## FastTree
  elif label in ["fasttree"]:
    for line in open(stats_file, "rU"):
      if line.lower().find("loglk") == -1:
        continue
      f = map(strip, line.split("\t"))
      try:
        value = float(f[2])
      except:
        continue
      logLK = value if not logLK or value < logLK else logLK

  ## RAXML
  for line in open(stats_file, "rU"):
    if not line.lower().startswith("final") or line.lower().find("score") == -1:
      continue
    logLK = float(map(strip, line.split())[-1])
    break

  ## Return the likelihood value for the current tree
  return logLK
