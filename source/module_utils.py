"""
  phylomizer - automated phylogenetic reconstruction pipeline - it resembles the
  steps followed by a phylogenetist to build a gene family tree with error-
  control of every step

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
## To guarantee compatibility with python3.4
from __future__ import print_function

import os
import sys
import subprocess as sp

from time import sleep

## Define a substitution for strip.
strip = lambda x: x.strip()

available_tags = ["binary", "file", "directory", "mode", "parameter", "other"]

verbose_levels = {
  "none":    0, "0": 0,
  "logfile": 1, "1": 1,
  "stderr":  2, "2": 2,
}

default_verbose = "logfile"

def parseComments(string_list):
  '''
  Generator used to parse a string list containing "#" symbols for comments
  '''
  comment = False
  for e in string_list:
    if e.startswith("#"):
      comment = True
    yield None if comment else e

def readConfig(input_file):
  '''
  Read pipeline configuration file checking whether all files, binaries and
  directories exist and can be accessed/executed
  '''

  toCheck, parameters = {}, {}
  for line in open(input_file, "rU"):
    ## Parse line upon the apperance of the first '#' symbol
    f = parseComments([e for e in map(strip, line.split()) if e])
    parsed_line = [element for element in f if element]

    ## Discard empty lines or those starting by "#"
    if not parsed_line:
      continue

    ## Get additional arguments in the input line
    args = "" if len(parsed_line) < 2 else " ".join(parsed_line[2:])

    ## Depending on the TAG for the current parameter, they need to be verify
    param = parsed_line[0]
    tag = parsed_line[1].lower()

    ## Check for the binary in the current PATH if not specific path is specific
    if tag == "binary":
      args = lookForProgram(param) if not args else lookForProgram(args)
      if not args:
        sys.exit(("ERROR: Impossible to find the binary for '%s'") % (param))
      args = os.path.abspath(args)

    ## Check whether current file exist
    elif tag == "file":
      if not lookForFile(args):
        sys.exit(("ERROR: Check your input file '%s'") % (args))
      args = os.path.abspath(args)

    ## Check directories exist and are writable
    elif tag == "directory":
      if not lookForDirectory(args):
        sys.exit(("ERROR: Check your input directory '%s'") % (args))
      args = os.path.abspath(args)

    ## Since the 'mode' tag define which programs should be executed in a given
    ## step, convert the arguments line into a line
    elif tag == "mode":
      args = args.split()

    ## We allow to run specific commands - for instance to use MPI - using the
    ## "other" tag - therefore we controll whether that tag has been used or not
    elif not tag in available_tags:
      msg =  "ERROR: Invalid tag - specific commands can be executed using "
      msg += "'other' tag.\nAvailable ones [%s]" % ", ".join(available_tags)
      sys.exit(msg)

    ## Depending whether the current parameter exists, assign or add the current
    ## arguments. On this way it is possible to have multi-line parameters.
    if param in parameters:
      args = " ".join([parameters[param], args])
    parameters[param] = args

  return parameters

def lookForProgram(binary):
  ''' Return if a given binary is on the file system
  '''
  try:
    pipe = sp.Popen(("which %s") % (binary), shell = True, stdout = sp.PIPE)
  except OSError as e:
    sys.exit(("ERROR: Impossible to find '%s'\nReport: %s") % (binary, str(e)))

  ## Parse the 'binary sequence of bytes' produce by communicate and convert it
  ## to a string
  program_path = "".join(list(map(strip, pipe.communicate()[0].decode('ascii'))))

  if program_path == "" or not lookForFile(program_path):
    return None
  return program_path

def lookForFile(input_file, sleep_interval = 5, attempts = 1):
  ''' Return if a given file exists or not
  '''

  ## If input file doesn't exit, don't check it
  if not input_file:
    return False

  ## We can delay to check whether a given file exists and contains something
  while not os.path.exists(input_file) and attempts:
    sleep(sleep_interval)
    attempts -= 1

  try:
    return (os.path.isfile(input_file) and os.path.getsize(input_file) > 0)
  except:
    return False

def lookForDirectory(input_direct, create = True):
  ''' Return whether a given folder exists or it has been created.
  '''

  ## If input directory exist, check whether is writable
  if os.path.isdir(input_direct):
    return os.access(input_direct, os.W_OK)

  ## If it doesn't exist, check whether it should be created or not
  if not create:
    return False

  try:
    sp.call(("mkdir -p %s") % (input_direct), shell = True)
  except OSError as e:
    sys.exit(("ERROR: Impossible to create '%s'\nReport: %s") % (input_direct, \
      str(e)))
  return True

def format_time(delta):
  ''' Format a given amount of elapsed seconds into a human readable notation
  '''

  ## We convert delta - difference between two timestamps - into seconds.
  ## Note that when delta has as arguments days, seconds, and microseconds so
  ## any time longer than a day has to be converted before parsing it.
  total_time = delta.days * 86400 + delta.seconds

  days, remaining = divmod(total_time, 86400)
  hours, remaining = divmod(remaining, 3600)
  minutes, seconds = divmod(remaining, 60)
  return ("%dd:%02dh:%02dm:%02ds") % (days, hours, minutes, seconds)

def listDirectory(directory, fileExtList):
  '''
  Get list of file info objects for files of particular extensions
  '''

  ## Add a "." to the file extension variable if it doesn't contain it
  if fileExtList[0] != ".":
    fileExtList = "." + fileExtList

  ## Keep only files excluding directories
  fileList = [os.path.normcase(f) for f in os.listdir(directory) \
    if os.path.isfile(f)]
  return [ os.path.join(directory, f) for f in fileList \
           if os.path.splitext(f)[1] in fileExtList ]

def splitSequence(seq, length = 80):
  ''' Split a given sequence contained in one line into lines of size "length"
  '''
  return "\n".join([seq[i:i + length] for i in range(0, len(seq), length)])

def sort_blast_hits(x, y):
  ''' Return 1, 0, -1 depending on the values comparison for BLAST hits result
      generated using the -m8 format
  '''

  if float(x[10]) > float(y[10]):
    return 1
  elif float(x[10]) < float(y[10]):
    return -1
  elif float(x[11]) < float(y[11]):
    return 1
  elif float(x[11]) > float(y[11]):
    return -1
  return 0

def sort_hmmer_hits(x, y):
  ''' Return 1, 0, -1 depending on the values comparison for HMMER hits result
      generated using the --tblout format
  '''

  try:
    if float(x[4]) > float(y[4]):
      return 1
    elif float(x[4]) < float(y[4]):
      return -1
    elif float(x[7]) < float(y[7]):
      return 1
    elif float(x[7]) > float(y[7]):
      return -1
  except:
    print(("x: %s\ty: %s") % (x, y))
  return 0

def printConfig(input_parameters, dest = sys.stderr):
  ''' Function for visualizing the selected parameters for running the pipeline
      or different steps of it
  '''

  ## Show which parameters has been set-up
  output = [("| %-24s\t| %s") % (("'%s'") % (key), value) for key, value in \
    sorted(input_parameters.items())]
  maxLen = sorted([len(l) for l in output])[-1] + 18

  print(("#%s#") % ("#" * maxLen), file = dest)
  print(("#%s#") % ("Pipeline Configuration".center(maxLen)), file = dest)
  print(("#%s#") % ("#" * maxLen), file = dest)
  print(("%s") % ("\n".join(output)), file = dest)
  print(("#%s#") % ("#" * maxLen), file = dest)
