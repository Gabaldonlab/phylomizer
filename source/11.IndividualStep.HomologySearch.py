#!/usr/bin/python
import os
import sys
import argparse
from module_utils import *
from module_homology import *

if __name__ == "__main__":

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", type = str, default = None,
    help = "Input file containing the query sequence/s")

  parser.add_argument("-d", "--db", dest = "dbFile", type = str, default = None,
    help = "Input file containing the target sequence database")

  parser.add_argument("-c", "--config", dest = "configFile", default = None, \
    type = str, help = "Input configuration file")

  parser.add_argument("-o", "--out", dest = "outFolder", type = str, default = \
    ".", help = "Output folder where all generated files will be dumped")

  parser.add_argument("-r", "--replace", dest = "replace", default = False, \
    action = "store_true", help = "Over-write any previously generated file")

  args = parser.parse_args()

  parameters = {}
  parameters.setdefault("replace", args.replace)

  ## Assign input parameters directly to the dictionary which will contain all
  ## current run configuration. Check all of them
  if not lookForFile(args.inFile):
    sys.exit(("ERROR: Check input QUERY SEQUENCE/s file '%s'") % (args.inFile))
  parameters.setdefault("in_file", args.inFile)

  if not lookForFile(args.dbFile):
    sys.exit(("ERROR: Check input TARGET SEQUENCES file '%s'") % (args.dbFile))
  parameters.setdefault("db_file", args.dbFile)

  if not lookForFile(args.configFile):
    sys.exit(("ERROR: Check input CONFIG file '%s'") % (args.configFile))
  parameters.setdefault("config_file", args.configFile)

  if not lookForDirectory(args.outFolder):
    sys.exit(("ERROR: Check output folder '%s'") % (args.outFolder))
  parameters.setdefault("out_directory", os.path.abspath(args.outFolder))

  ## Read the other parameters from the input config file
  parameters.update(readConfig(parameters["config_file"]))

  ## Check specific values for input parameters.
  if not "coverage" in parameters or not (0.0 < float(parameters["coverage"]) \
    <= 1.0):
    sys.exit(("ERROR: Check your 'coverage' parameter"))

  if not "hits" in parameters or int(parameters["hits"]) < 1:
    sys.exit(("ERROR: Check your 'hits' upper limit value"))

  ## Launch the whole homology process
  homology(parameters)
