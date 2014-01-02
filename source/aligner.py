import sys, os, subprocess as sp
from operator import itemgetter
from utils import *

programs = ["muscle", "mafft", "dialigntx"]
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def AlignerPipeline(parameters):

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  iFolder, iFile = os.path.split(parameters["inFile"])
  iFolder = "." if iFolder == "" else iFolder
  oFile = os.path.join(parameters["outDirec"], iFile.split(".")[0])
  os.chdir(iFolder)
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  numSeqs, Sel = CheckSequences(parameters["inFile"])
  lrepl = False
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  if not "both_direction" in parameters:
    parameters["both_direction"] = True

  if type(parameters["both_direction"]) == "str":
    if parameters["both_direction"].lower() == "false":
      parameters["both_direction"] = False

  if not "single_program" in parameters:
    parameters["single_program"] = programs

  else:
    if not parameters["single_program"].lower() in programs:
      raise NameError(("Program not available '%s'") % parameters["single_program"])
    parameters["single_program"] = [parameters["single_program"].lower()]

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if numSeqs < 3 and parameters["verbose"] > 0:
    sys.exit("INFO: It is necessary, at least, 3 sequences")

  if Sel == True:
    if parameters["verbose"] > 0:
      print "At least, one selenocysteine residue has been detected"

    lrepl = ChangeResiduesSeq(parameters["inFile"], oFile + ".seqs.noselcys", \
      "U", parameters["in_letter"], parameters["replace"])

    if lrepl:
      parameters["replace"] = True

    if parameters["both_direction"] == True:
      lrepl = GettingReverse(parameters["readal"], oFile + ".seqs.noselcys", \
        oFile + ".seqs.noselcys.reverse", oFile + ".log", parameters["replace"])

    # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** **
    # Muscle
    if "muscle" in parameters["single_program"]:
      if parameters["both_direction"] == True:
        lrepl = AligningMuscle(parameters["muscle"], oFile + ".seqs.noselcys" +\
          ".reverse", oFile + ".alg.noselcys.reverse.msl", oFile + ".log", \
          parameters["muscle_params"], parameters["replace"])

        lrepl = ChangeResiduesSeq(oFile + ".alg.noselcys.reverse.msl", oFile + \
          ".alg.reverse.msl", parameters["in_letter"], "U", parameters["replace"])

      lrepl = AligningMuscle(parameters["muscle"], oFile + ".seqs.noselcys", \
        oFile + ".alg.noselcys.forward.msl", oFile + ".log", \
        parameters["muscle_params"], parameters["replace"])

      lrepl = ChangeResiduesSeq(oFile + ".alg.noselcys.forward.msl", oFile + \
        ".alg.forward.msl", parameters["in_letter"], "U", parameters["replace"])
    # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** **

    # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** **
    # Mafft
    if "mafft" in parameters["single_program"]:
      if parameters["both_direction"] == True:
        lrepl = AligningMafft(parameters["mafft"], oFile + ".seqs.noselcys." + \
          "reverse", oFile + ".alg.noselcys.reverse.mft", oFile + ".log", \
          parameters["mafft_params"], parameters["replace"])

        lrepl = ChangeResiduesSeq(oFile + ".alg.noselcys.reverse.mft", oFile + \
          ".alg.reverse.mft", parameters["in_letter"], "U", parameters["replace"])

      lrepl = AligningMafft(parameters["mafft"], oFile + ".seqs.noselcys", \
        oFile + ".alg.noselcys.forward.mft", oFile + ".log", \
        parameters["mafft_params"], parameters["replace"])

      lrepl = ChangeResiduesSeq(oFile + ".alg.noselcys.forward.mft", oFile + \
        ".alg.forward.mft", parameters["in_letter"], "U", parameters["replace"])
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if Sel == False:

    if parameters["both_direction"] == True:
      lrepl = GettingReverse(parameters["readal"], parameters["inFile"], oFile \
        + ".seqs.reverse", oFile + ".log", parameters["replace"])

    if lrepl:
      parameters["replace"] = True

    if "muscle" in parameters["single_program"]:
      lrepl = AligningMuscle(parameters["muscle"], parameters["inFile"], oFile \
        + ".alg.forward.msl", oFile + ".log", parameters["muscle_params"],
        parameters["replace"])

      if parameters["both_direction"] == True:
        lrepl = AligningMuscle(parameters["muscle"], oFile + ".seqs.reverse",
          oFile + ".alg.reverse.msl", oFile + ".log", parameters["muscle_params"],
          parameters["replace"])

    if "mafft" in parameters["single_program"]:
      lrepl = AligningMafft(parameters["mafft"], parameters["inFile"], oFile + \
        ".alg.forward.mft", oFile + ".log", parameters["mafft_params"],
        parameters["replace"])

      if parameters["both_direction"] == True:
        lrepl = AligningMafft(parameters["mafft"], oFile + ".seqs.reverse",
          oFile + ".alg.reverse.mft", oFile + ".log", parameters["mafft_params"],
          parameters["replace"])
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if parameters["both_direction"] == True:
    lrepl = GettingReverse(parameters["readal"], parameters["inFile"], oFile + \
      ".seqs.reverse", oFile + ".log", parameters["replace"])

  if lrepl:
    parameters["replace"] = True

  if "dialigntx" in parameters["single_program"]:
    if parameters["both_direction"] == True:
      lrepl = AligningDialgnTX(parameters["dialign-tx"], oFile + \
        ".seqs.reverse", oFile + ".alg.reverse.dtx", oFile + ".log", \
        parameters["dialigntx_params"], parameters["replace"])

    lrepl = AligningDialgnTX(parameters["dialign-tx"], parameters["inFile"], \
      oFile + ".alg.forward.dtx", oFile + ".log", parameters["dialigntx_params"],
      parameters["replace"])

  if parameters["both_direction"] == True:

    if "muscle" in parameters["single_program"]:
      lrepl = GettingReverse(parameters["readal"], oFile + '.alg.reverse.msl', \
        oFile + ".alg.reverse.forw.msl", oFile + ".log", parameters["replace"])

    if "mafft" in parameters["single_program"]:
      lrepl = GettingReverse(parameters["readal"], oFile + '.alg.reverse.mft', \
        oFile + ".alg.reverse.forw.mft", oFile + ".log", parameters["replace"])

    if "dialigntx" in parameters["single_program"]:
      lrepl = GettingReverse(parameters["readal"], oFile + '.alg.reverse.dtx', \
        oFile + ".alg.reverse.forw.dtx", oFile + ".log", parameters["replace"])
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  # If the current pipeline is used then perform all of these operations
  if len(parameters["single_program"]) > 1 or parameters["both_direction"] == True:

    pathFile = open(oFile + ".alg.paths", "w")
    direction = ['forward']
    if parameters["both_direction"] == True:
      direction.append('reverse.forw')

    for d in direction:
      for program in ['msl', 'mft', 'dtx']:
        print >> pathFile, oFile + ".alg." + d + "." + program
    pathFile.close()

    if lrepl:
      parameters["replace"] = True

    MetaAligningMCoffee(parameters["t_coffee"], parameters["inFile"], oFile + \
      ".alg.metalig", oFile + ".log",  oFile + ".alg.paths",
      parameters["mcoffee_params"], parameters["replace"])

    convertAlignment(parameters["trimal"], oFile + ".alg.metalig", oFile + \
      ".alg.metalig", "phylip", oFile + ".log", parameters["replace"])

    trimmingAlignment(parameters["trimal"], None, oFile + ".alg.clean", oFile \
      + ".log", oFile + ".alg.metalig", oFile + ".alg.paths",
      parameters["trimal_compare"], parameters["trimal_params"],
      parameters["replace"])

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  # Otherwise, perform another operations
  else:

    if "muscle" in parameters["single_program"]:
      ext = ".alg.forward.msl"
    elif "mafft" in parameters["single_program"]:
      ext = ".alg.forward.mft"
    elif "dialigntx" in parameters["single_program"]:
      ext = ".alg.forward.dtx"

    convertAlignment(parameters["trimal"], oFile + ext, oFile + \
      ".alg.metalig", "phylip", oFile + ".log", parameters["replace"])

    trimmingAlignment(parameters["trimal"], oFile + ".alg.metalig", oFile + \
      ".alg.clean", oFile + ".log", "", "", "", parameters["trimal_params"],
      parameters["replace"])

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  return oFile + ".alg.clean"
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def GettingReverse(bin, inFile, outFile, logFile, replace):
  '''
  Get the reverse alignment using readal for that
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if lookForFile(outFile) and not replace: return;
  lgFile = open(logFile, "a+ ") if logFile != "" else None
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  cmd = bin + " -in " + inFile + " -out " + outFile + " -reverse"
  try: proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:   sys.exit("ERROR: Execution failed: " + str(e))
  if proc.wait() != 0: sys.exit("ERROR: Execution failed: readAl")
  lgFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def convertAlignment(bin, inFile, outFile, format, logFile, replace):
  '''
  Convert the input alignmnet into the input format
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if lookForFile(outFile) and not replace: return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  cmd = bin + " -out " + outFile + " -in " + inFile + " -" + format
  try: proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:   sys.exit("ERROR: Execution failed: " + str(e))
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if proc.wait() != 0: sys.exit("ERROR: Execution failed: trimAl")
  lgFile.close()
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  return True

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def trimmingAlignment(bin, inFile, outFile, logFile, forceSelection,
  compareSelection, compareParameters, parameters, replace):
  '''
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if lookForFile(outFile) and not replace: return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  cmd = ""
  if compareSelection:
    cmd  = " -compareset " + compareSelection
  if forceSelection:
    cmd += " -forceselect " + forceSelection
  if compareParameters:
    cmd += " " + compareParameters
  if inFile:
    cmd += " -in " + inFile
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  cmd = bin + " -out " + outFile + " " + parameters + cmd
  try: proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:   sys.exit("ERROR: Execution failed: " + str(e))
  if proc.wait() != 0: sys.exit("ERROR: Execution failed: trimAl")
  lgFile.close()

  return True
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****


### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def AligningMuscle(bin, inFile, outFile, logFile, parameters, replace):
  '''
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if lookForFile(outFile) and not replace: return False;
  lgFile = open(logFile, "a+ ") if logFile != "" else None
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  cmd = bin + " " + parameters + " -in " + inFile + " -out " + outFile
  try: proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:    sys.exit("ERROR: Execution failed: " + str(e))
  if proc.wait() != 0: sys.exit("ERROR: Execution failed: MUSCLE")
  lgFile.close()

  return True
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def AligningMafft(bin, inFile, outFile, logFile, parameters, replace):
  '''
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if lookForFile(outFile) and not replace: return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  cmd = bin + " " + parameters + " " + inFile + " > " + outFile
  try: proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:   sys.exit("ERROR: Execution failed: " + str(e))
  if proc.wait() != 0: sys.exit("ERROR: Execution failed: MAFFT")
  lgFile.close()

  return True
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def AligningDialgnTX(bin, inFile, outFile, logFile, parameters, replace):
  '''
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if lookForFile(outFile) and not replace: return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  cmd = bin + " " + parameters + " " + inFile + " " + outFile
  try: proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:   sys.exit("ERROR: Execution failed: " + str(e))
  if proc.wait() != 0: sys.exit("ERROR: Execution failed: DialignTX")
  lgFile.close()

  return True
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def MetaAligningMCoffee(bin, inFile, outFile, logFile, pathFile, pars, replace):
  '''
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if lookForFile(outFile) and not replace: return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  ## Define some environment variables for t_coffee
  # user = os.getenv("LOGNAME")
  # try:
  #   proc = sp.Popen(("mkdir -p -m0777 /tmp/tcoffee/%s;") % (user), shell = True)
  #   os.putenv("LOCKDIR_4_TCOFFEE", ("/tmp/tcoffee/%s") % (user))
  #   os.putenv("TMP_4_TCOFFEE", ("/tmp/tcoffee/%s") % (user))
  # except OSError, e:
  #   pass

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  pFile, files = open(pathFile, "rU"), " -aln "
  for line in pFile.readlines(): files += line.strip() + " "
  pFile.close()
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  cmd = bin + " " + inFile + " -outfile " + outFile + " " + pars + files
  try: proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:   sys.exit("ERROR: Execution failed: " + str(e))
  if proc.wait() != 0:
    print ("COMMAND-LINE: %s") % (cmd)
    sys.exit("ERROR: Execution failed: T-Coffee")
  lgFile.close()

  sp.call("rm -f " + os.path.split(inFile)[1][:10] + "*.dnd", shell = True)

  return True
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def CheckSequences(inFile):
  '''
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  iFile, nSeqs, selCys = open(inFile, "rU"), 0, False
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  for line in iFile.readlines():
    if line[0] == ">": nSeqs += 1
    elif line.strip().upper().find("U") != -1: selCys = True
  iFile.close()
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  return nSeqs, selCys
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def ChangeResiduesSeq(inFile, outFile, inLetter, outLetter, replace):

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if lookForFile(outFile) and not replace: return False;
  oFile, iFile = open(outFile, "w"), open(inFile, "rU")
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  for line in iFile.readlines():
    if line[0] == ">": print >> oFile, line.strip()
    else: print >> oFile, line.replace(inLetter, outLetter).strip()
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  oFile.close()
  iFile.close()
  return True
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
