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

import sys, os, argparse
import blaster, aligner, tree_reconstructor as tree, utils

def main():

  usage = "%(prog)s -i seed_sequence -c config_file -d output_directory -b " \
    + "blast_db [options]\n"

  parser = argparse.ArgumentParser(usage = usage, description = desc,
    epilog = epilog)

  ## Capture initial parameters
  parser.add_argument("-i", "--in", dest = "inFile", type = str, required = True,
    help = "Input FASTA file containing the seed sequence which will be used to"
    + " perform the homology search")

  parser.add_argument("-c", "--config", type = str, required = True, help = \
    "Input <TAB> delimited file containing the pipeline configuration")

  parser.add_argument("-b", "--db", type = str, required = True, help = \
    "Input Sequences DB file containing all sequences which will be scanned "
    + "looking for homologs")

  parser.add_argument("-o", "--outdir", default = ".", type = str, help = \
    "Define a directory where output files will be stored")

  parser.add_argument("-r", "--replace", default = False, action = "store_true",
    help = "Overwrite any previously generated file on the output directory")

  parser.add_argument("-v", "--verbose", default = False, action = "store_true")
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

  parameters.setdefault("replace", args.replace)
  parameters.setdefault("verbose", args.verbose)

  print str(parameters)

  ## Read Configuration file and update current structure with the different
  ## parameters
  parameters.update(utils.readConfig(args.config))

  print "---"
  print str(parameters)



  #~ # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
#~
  #~ # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  #~ print "Pipeline Configuration\n---"
  #~ for key in parameters:
    #~ print key.center(20) + "\t" + str(parameters[key])
  #~ print "---"
  #~ # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
#~
  #~ # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  #~ if float(parameters["coverage"]) > 1 or float(parameters["coverage"]) <= 0:
    #~ sys.exit("ERROR: The coverage has to be defined between 0 and 1")
#~
  #~ if int(parameters["hits"]) < 1:
    #~ sys.exit("ERROR: The Number of hits for the Blast Output has to be greater"\
             #~ + "than 0")
#~
  #~ if int(parameters["numb_models"]) > len(parameters["evol_models"].split(" "))\
    #~ or int(parameters["numb_models"]) < 1:
    #~ sys.exit("ERROR: The number of models to be evaluated have to be defined "
             #~ "between 1 and " + str(len(parameters["evol_models"].split(" "))))
  #~ # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
#~
  #~ # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  #~ parameters["inFile"] = blast(parameters)
#~
  #~ parameters["inFile"] = AlignerPipeline(parameters)
#~
  #~ PhylogeneticTrees(parameters["phyml"], parameters["nj_parameters"], "nj",
    #~ parameters["inFile"], parameters["outDirec"], parameters["evol_models"],
    #~ parameters["verbose"], parameters["replace"])
#~
  #~ rank = SortingLKs(parameters["inFile"], parameters["outDirec"], "nj",
    #~ parameters["evol_models"], parameters["verbose"])
#~
  #~ rank = ' ' . join(rank[:int(parameters["numb_models"])])
#~
  #~ PhylogeneticTrees(parameters["phyml"], parameters["ml_parameters"], "ml",
    #~ parameters["inFile"], parameters["outDirec"], rank, parameters["verbose"],
    #~ parameters["replace"])
#~
  #~ SortingLKs(parameters["inFile"], parameters["outDirec"], "ml", rank,
    #~ parameters["verbose"])
  #~ # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
#~
  #~ return 0
#~
## ***** ***** ***** ***** ***** ***** ***** ***** *****
if __name__ == "__main__":
  sys.exit(main())
