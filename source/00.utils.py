import sys, os, subprocess as sp
from operator import itemgetter
from getpass import getpass
from string import strip

def readConfig(input_file):
  ''' Read pipeline configuration file checking whether all files, binaries and
      directories exist and can be accessed/executed
  '''

  toCheck, parameters = {}, {}
  for line in open(input_file, "rU"):
    f = map(strip, line.split())

    ## Discard empty lines or with comments
    if f[0].startswith("#") or (len(f) == 1 and f[0] == ""):
      continue

    ## Comments at the end of the line are allowed
    arg = []
    for e in f[2:]:
      if e.startswith("#"):
        break
      arg.append(e)
    arg = "" if arg == [] else " ".join(arg)

    dtype = f[1].lower()
    if dtype in ["binary", "files", "directory"]:
      toCheck.setdefault(dtype, []).append((f[0], arg))
    else:
      parameters[f[0]] = arg

  ## Check whether binaries are accessible
  if "binary" in toCheck:
    for binary,path in toCheck["binary"]:
      progr_path = lookForProgram(binary) if path == "" else lookForProgram(path)
      if progr_path == None:
        sys.exit(("ERROR: Impossible to find the binary file '%s'") % (binary))
      parameters[binary] = progr_path

  ## Check for any file included in the configuration file
  if "files" in toCheck:
    for key,infile in toCheck["files"]:
      if not lookForFile(infile):
        sys.exit(("ERROR: Impossible to find the input file for '%s'") % (key))
      parameters[key] = infile

  ## Check all directories exist and are accessible
  if "directory" in toCheck:
    for key,direc in toCheck["directory"]:
      if not lookForDirectory(direc):
        sys.exit(("ERROR: Impossible to access to '%s' directory") % (key))
      parameters[key] = direc

  return parameters

def lookForProgram(binary):
  ''' Return if a given binary is on the file system
  '''
  try:
    pipe = sp.Popen(("which %s") % (binary), shell = True, stdout = sp.PIPE)
  except OSError, e:
    sys.exit(("ERROR: Impossible to find '%s'\nReport: %s") % (binary, str(e)))

  program_path = "".join(map(strip, pipe.stdout.readlines()))

  if program_path == "" or not lookForFile(program_path):
    return None
  return program_path

def lookForFile(input_file):
  ''' Return if a given file exists or not
  '''
  return (os.path.isfile(input_file) and os.path.getsize(input_file) > 0)

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
  except OSError, e:
    sys.exit(("ERROR: Impossible to create '%s'\nReport: %s") % (input_direct, \
      str(e)))
  return True

def format_time(total_time):
  ''' Format a given amount of elapsed seconds into a human readable notation
  '''

  days, remaining = divmod(total_time, 86400)
  hours, remaining = divmod(remaining, 3600)
  minutes, seconds = divmod(remaining, 60)
  return ("%dd:%02dh:%02dm:%02ds") % (days, hours, minutes, seconds)

def listDirectory(directory, fileExtList):
  '''
  Get list of file info objects for files of particular extensions
  '''
  fileList = [os.path.normcase(f) for f in os.listdir(directory)]
  return [ os.path.join(directory, f) for f in fileList \
           if os.path.splitext(f)[1] in fileExtList ]

def splitSequence(seq, length = 80):
  ''' Split a given sequence contained in one line into lines of size "length"
  '''
  return "\n".join([seq[i:i + length] for i in range(0, len(seq), length)])

def sort_hits(x, y):
  ''' Return 1, 0, -1 depending on the values comparison for BLAST hits result
      generated using the -m8 format
  '''

  if float(x[10]) > float(y[10]):
    return 1
  elif float(x[10]) < float(y[10]):
    return -1;
  elif float(x[11]) < float(y[11]):
    return 1
  elif float(x[11]) > float(y[11]):
    return -1
  return 0
