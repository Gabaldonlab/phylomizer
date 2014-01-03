#!/usr/bin/env python

desc="""
 Phylogenetic tree reconstruction pipeline. It comprises three main steps:
   1) Homology search using tools such as BLAST.
   2) Multiple Sequence Alignment (MSA) including the usage of different
      aligners and the generation of alignments in different directions,
      the generation of a meta-alignment and the trimming of these meta-
      MSA using information from individual alignments.
   3) Phylogenetic tree reconstruction using fast model selection over NJ
      trees. It is possible to reconstruct trees using AA, NT or Codons.
"""

epilog="""Author: Salvador Capella-Gutierrez / salcagu@gmail.com
Barcelona, 02/01/2013
"""

"""
Versions:
- 1.00:
  > Initial release
"""

import sys, os, argparse, datetime
import homology, aligner, tree_reconstructor as tree, utils

def main():

  usage = "%(prog)s -i seed_sequence -c config_file -d output_directory -b " \
    + "blast_db [options]\n"

  parser = argparse.ArgumentParser(usage = usage, description = desc,
    epilog = epilog)

  ## Capture initial parameters
  parser.add_argument("-i","--in", dest = "inFile", type = str, required = True,
    help = "Input FASTA file containing the seed sequence which will be used to"
    + " perform the homology search")

  parser.add_argument("-c", "--config", type = str, required = True, help = \
    "Input <TAB> delimited file containing the pipeline configuration")

  parser.add_argument("-b", "--db", type = str, required = True, help = \
    "Input Sequences DB file containing all sequences which will be scanned "
    + "looking for homologs")

  parser.add_argument("-l", "--log", default = "", type = str, help = "Define" \
    + "a log file where all messages will be dumped. Otherwise it will be set "
    + "automatically using the output directory and the input file name")

  parser.add_argument("-o", "--outdir", default = os.getcwd(), type = str, \
    help = "Define a directory where output files will be stored")

  parser.add_argument("-r", "--replace", default = False, action = "store_true",
    help = "Overwrite any previously generated file on the output directory")

  parser.add_argument("--compress", default = False, action = "store_true",
    help = "Compress intermediate files on .gz format")

  parser.add_argument("-v", "--verbose", default = True, action = "store_false")
  parser.add_argument("--version", action = "version", version ='%(prog)s v1.0')

  args = parser.parse_args()

  parameters = {}
  ## Check input files and, if OK, assign them to an appropiate structure

  if not os.path.isfile(args.inFile):
    sys.exit(("\nERROR: Please check your input SEED Sequence file '%s'\n") % \
      (args.inFile))
  parameters.setdefault("inFile", args.inFile)

  if not os.path.isfile(args.config):
    sys.exit(("\nERROR: Please check your input Pipeline Config file '%s'\n") \
      % (args.config))
  parameters.setdefault("config", args.config)

  if not os.path.isfile(args.db):
    sys.exit(("\nERROR: Please check your input Sequences DB file '%s'\n") % \
      (args.db))
  parameters.setdefault("db", args.db)

  if not utils.lookForDirectory(args.outdir):
    sys.exit(("\nERROR: Check your output directory '%s'") % (args.outdir))
  parameters.setdefault("outdirec", args.outdir)

  parameters.setdefault("compress", args.compress)
  parameters.setdefault("replace", args.replace)
  parameters.setdefault("verbose", args.verbose)

  ## Set the log file
  if not args.log:
    ## Determine automatically the log file name using output directory and
    ## input file name
    prefixFile = os.path.split(parameters["inFile"])[1].split(".")[0]
    refFile = ("%s.log") % (os.path.join(parameters["outdirec"], prefixFile))
  else:
    refFile = args.log
  ## If log information is requested, open the output stream
  parameters["log"] = open(refFile, "w" if parameters["replace"] else "a+") \
    if parameters["verbose"] else None

  ## Read Configuration file and update current structure with the different
  ## parameters
  parameters.update(utils.readConfig(args.config))

  ## Print configuration pipeline
  if parameters["verbose"]:
    print >> parameters["log"], ("### %s\n%s") % \
      ("Pipeline General Configuration".center(90), "#" * 90)
    for key,value in sorted(parameters.iteritems()):
      print >> parameters["log"], ("### %30s\t%s") % (key.center(30), value)
    print >> parameters["log"], ("%s") % ("#" * 90)

  ## Register when each process starts
  general_start = datetime.datetime.now()
  if parameters["verbose"]:
    date = general_start.strftime("%H:%M:%S %m/%d/%y")
    print >> parameters["log"], ("### Pipeline starts\n### Homology search\n"
      + "### %s") % (date)

  ## Homology search
  step_start = general_start
  resulting_file = ""
  #~ resulting_file = homology.blast(parameters)
  step_end = datetime.datetime.now()

  if parameters["verbose"]:
    date = step_end.strftime("%H:%M:%S %m/%d/%y")
    total = utils.format_time((step_end - step_start).seconds if step_start \
      else 0)
    whole = utils.format_time((step_end - general_start).seconds \
      if general_start else 0)
    print >> parameters["log"], ("### Pipeline progression\n### Homology search"
      + "\n### %s\n### Step \t%s\n### Total\t%s") % (date, total, whole)

  ## Multiple Sequence Alignments
  parameters["inFile"] = resulting_file
  step_start = datetime.datetime.now()
  #~ resulting_file = aligner.AlignerPipeline(parameters)
  step_end = datetime.datetime.now()

  if parameters["verbose"]:
    date = step_end.strftime("%H:%M:%S %m/%d/%y")
    total = utils.format_time((step_end - step_start).seconds if step_start \
      else 0)
    whole = utils.format_time((step_end - general_start).seconds \
      if general_start else 0)
    print >> parameters["log"], ("### Pipeline progression\n### Multiple Sequen"
      "ce Alignment\n### %s\n### Step \t%s\n### Total\t%s") % (date,total,whole)

  ## Fast evolutionary model selection based on NJ-trees
  parameters["inFile"] = resulting_file
  step_start = datetime.datetime.now()
  #~ resulting_file = tree.ModelSelection(parameters)
  step_end = datetime.datetime.now()

  #~ PhylogeneticTrees(parameters["phyml"], parameters["nj_parameters"], "nj",
    #~ parameters["inFile"], parameters["outDirec"], parameters["evol_models"],
    #~ parameters["verbose"], parameters["replace"])

  if parameters["verbose"]:
    date = step_end.strftime("%H:%M:%S %m/%d/%y")
    total = utils.format_time((step_end - step_start).seconds if step_start \
      else 0)
    whole = utils.format_time((step_end - general_start).seconds \
      if general_start else 0)
    print >> parameters["log"], ("### Pipeline progression\n### Fast evolutiona"
      + "ry model selection based on Neighbour-Joining trees\n### %s\n### Step "
      + "\t%s\n### Total\t%s") % (date, total, whole)

  ## Depending on whether the model selection has been performed using NJ trees
  ## or the user wants to reconstruct all models using a Maximum-Likelihood
  ## approach
  #~ rank = SortingLKs(parameters["inFile"], parameters["outDirec"], "nj",
    #~ parameters["evol_models"], parameters["verbose"])

  #~ rank = ' ' . join(rank[:int(parameters["numb_models"])])

  ## Phylogenetic tree reconstruction based on a Maximum Likelihood framework as
  ## implemented in PhyML, RAxML or FastML
  parameters["inFile"] = resulting_file
  step_start = datetime.datetime.now()
  #~ PhylogeneticTrees(parameters["phyml"], parameters["ml_parameters"], "ml",
    #~ parameters["inFile"], parameters["outDirec"], rank, parameters["verbose"],
    #~ parameters["replace"])
  #~ resulting_file = tree.ModelSelection(parameters)
  step_end = datetime.datetime.now()

  if parameters["verbose"]:
    date = step_end.strftime("%H:%M:%S %m/%d/%y")
    total = utils.format_time((step_end - step_start).seconds if step_start \
      else 0)
    whole = utils.format_time((step_end - general_start).seconds \
      if general_start else 0)
    print >> parameters["log"], ("### Pipeline progression\n### Phylogenetic "
      + "tree reconstruction based on a Maximum-Likelihood framework\n### %s\n"
      + "### Step \t%s\n\n### Total Pipeline\t%s\n") % (date, total, whole)

  #~ SortingLKs(parameters["inFile"], parameters["outDirec"], "ml", rank,
    #~ parameters["verbose"])

  ## Close the log output stream
  if parameters["log"]:
    parameters["log"].close()

  return 0

if __name__ == "__main__":
  sys.exit(main())
