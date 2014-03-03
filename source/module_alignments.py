import sys, os, subprocess as sp
from operator import itemgetter
from string import strip
from time import sleep
from utils import *
import datetime

programs = { "muscle":"msl", "mafft":"mft", "dialign-tx":"dtx", "kalign":"kal" }

## Exit code meanings
##  80: Not enough sequences
##  81: Problems making format conversions/reversion sequences
##  82: Problems trimming input alignment
##  86: Problems aligning with MUSCLE
##  87: Problems aligning with MAFFT
##  88: Problems aligning with DiAlign-TX
##  89: Problems aligning with KAlign
##  90: Problems aligning M-Coffee

def alignment(parameters):

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
  print >> logFile, ("###\n###\tSTEP\tMultiple Sequence Alignment\tSTART\t%s"
    + "\n###") % (date)
  logFile.flush()

  ## Get some information such as number of input sequences and the presence of
  ## selenocysteine/pyrrolysine residues
  numSeqs, selenocys, pyrrolys = check_count_sequences(parameters["in_file"])




  lrepl = False







  print >> logFile, ("###\n###\tSTEP\tMultipple Sequence Alignment\tEND\t"
    + "%s") % (date)
  total = format_time((final - start).seconds if start else 0)
  print >> logFile, ("###\tTOTAL Time\tMultiple Sequence Alignment\t%s"
    + "\n###") % (total)
  logFile.close()

  ## Update the input file parameter and return the dictionary containing all
  ## parameters. Those parameters may be used in other steps
  parameters["in_file"] = outFile

  ## Before returning to the main program, get back to the original working
  ## directory
  os.chdir(current_directory)

  return parameters



def check_count_sequences(in_file):
  '''
  Given a set of sequences, return how many there are as well as if any
  selenocysteine or pyrrolysine is detected among the sequences
  '''

  numb_sequences, selenocys, pyrrolys = 0, False, False
  for record in SeqIO.parse(in_file, "fasta"):
    numb_sequences += 1

    seq = str(record.seq).upper()
    ## Detect any ocurrence of selenocysteine/pyrrolysine residues
    if seq.count("U") != 0:
      selenocys = True
    if seq.count("O") != 0:
      pyrrolys = True

  return numb_sequences, selenocys, pyrrolys













def AlignerPipeline(parameters):

  ## Get input folder and seed name

  ## Get output folder
  oFile = os.path.join(parameters["outDirec"], iFile.split(".")[0])

  ## Change current directory to the input folder
  os.chdir(iFolder)

  ## Start counting how much time will cost making an alignment
  start = datetime.datetime.now()


  ## Check whether alignments should be performed in forward or forward/reverse
  ## directions
  if not "both_direction" in parameters:
    parameters["both_direction"] = True

  elif type(parameters["both_direction"]) == "str":
    if parameters["both_direction"].lower() == "false":
      parameters["both_direction"] = False

  ## Check whether how many programs should be used to make the alignment
  ## reconstruction step. Right now we support up to 4 programs:
  ## Muscle, Mafft, DiAlign-TX & KAlign2
  if not "single_program" in parameters:
    parameters["single_program"] = []
    for param in parameters:
      if param.lower() in programs:
        parameters["single_program"].append(param.lower())

  ## In case just one program is selected, get which one among those availables
  else:
    if not parameters["single_program"].lower() in programs:
      raise NameError(("Program not available '%s'") % parameters["single_program"])
    parameters["single_program"] = [parameters["single_program"].lower()]

  ## Finish when there are not enough sequences to make an alignment
  if numSeqs < 3 and parameters["verbose"] > 0:
    date = start.strftime("%H:%M:%S %m/%d/%y")
    lgFile = open(oFile + ".log", "a+")
    print >> lgFile, ("### %s\n### INFO: It is necessary, "
      + "at least, 3 sequences (%d)") % (date, numSeqs)
    lgFile.close()
    sys.exit(80)

  ## Predefined input file/s. They may be changed depending on the presence/
  ## absence of selenocysteines
  currentInFile = parameters["inFile"]
  if parameters["both_direction"] == True:
    revCurrentInFile = ("%s.seqs.reverse") % (oFile)

    ## Get reversed input sequences
    lrepl = GettingReverse(parameters["readal"], parameters["inFile"],
      revCurrentInFile, oFile + ".log", parameters["replace"])

    ## If the file has been generated "de novo", regenerate the rest of files
    if lrepl:
      parameters["replace"] = True

  ## In those cases where a Selenocysteines have been detected, change them by
  ## wildcard characters
  if Sel == True:
    if parameters["verbose"] > 0:
      print "At least, one selenocysteine residue has been detected"
      print >> open(oFile + ".log", "a+"), ("### At least one SelenoCysteine "
        + "residue has been detected")

    lrepl = ChangeResiduesSeq(parameters["inFile"], oFile + ".seqs.noselcys", \
      "U", parameters["in_letter"], parameters["replace"])

    ## Certain programs need wildcards instead of SelenoCysteines (U) for
    ## working properly
    currentInFile = ("%s.seqs.noselcys") % (oFile)

    ## If the file has been generated "de novo", regenerate the rest of files
    if lrepl:
      parameters["replace"] = True

    if parameters["both_direction"] == True:
      lrepl = GettingReverse(parameters["readal"], oFile + ".seqs.noselcys", \
        oFile + ".seqs.noselcys.reverse", oFile + ".log", parameters["replace"])
      revCurrentInFile = ("%s.seqs.noselcys.reverse") % (oFile)

  ## Align input sequences depending on the selected programs
  ## Muscle
  if "muscle" in parameters["single_program"]:

    ## Alignment will be made in forward/reverse direction
    if parameters["both_direction"] == True:
      output =  ("%s.alg") % (oFile)
      output += ("%s.reverse.msl") % (".noselcys" if Sel == True else "")

      lrepl = AligningMuscle(parameters["muscle"], revCurrentInFile, output,
        oFile + ".log", parameters["muscle_params"], parameters["replace"])

      if Sel == True:
        realOutput = ("%s.alg.reverse.msl") % (oFile)
        lrepl = ChangeResiduesSeq(output, realOutput, parameters["in_letter"],
          "U", parameters["replace"])

      lrepl = GettingReverse(parameters["readal"], oFile + '.alg.reverse.msl', \
        oFile + ".alg.reverse.forw.msl", oFile + ".log", parameters["replace"])

    output =  ("%s.alg") % (oFile)
    output += ("%s.forward.msl") % (".noselcys" if Sel == True else "")

    lrepl = AligningMuscle(parameters["muscle"], currentInFile, output, \
      oFile + ".log", parameters["muscle_params"], parameters["replace"])

    if Sel == True:
      realOutput = ("%s.alg.forward.msl") % (oFile)
      lrepl = ChangeResiduesSeq(output, realOutput, \
        parameters["in_letter"], "U", parameters["replace"])

  ## KAlign
  if "kalign" in parameters["single_program"]:

    ## Alignment will be made in forward/reverse direction
    if parameters["both_direction"] == True:
      output =  ("%s.alg") % (oFile)
      output += ("%s.reverse.kal") % (".noselcys" if Sel == True else "")

      lrepl = AligningKAlign(parameters["kalign"], revCurrentInFile, output,
        oFile + ".log", parameters["kalign_params"], parameters["replace"])

      if Sel == True:
        realOutput = ("%s.alg.reverse.kal") % (oFile)
        lrepl = ChangeResiduesSeq(output, realOutput, parameters["in_letter"],
          "U", parameters["replace"])

      lrepl = GettingReverse(parameters["readal"], oFile + '.alg.reverse.kal', \
        oFile + ".alg.reverse.forw.kal", oFile + ".log", parameters["replace"])

    output =  ("%s.alg") % (oFile)
    output += ("%s.forward.kal") % (".noselcys" if Sel == True else "")

    lrepl = AligningKAlign(parameters["kalign"], currentInFile, output, \
      oFile + ".log", parameters["kalign_params"], parameters["replace"])

    if Sel == True:
      realOutput = ("%s.alg.forward.kal") % (oFile)
      lrepl = ChangeResiduesSeq(output, realOutput, \
        parameters["in_letter"], "U", parameters["replace"])

  ## Mafft
  if "mafft" in parameters["single_program"]:

    ## Alignment will be made in forward/reverse direction
    if parameters["both_direction"] == True:
      output =  ("%s.alg") % (oFile)
      output += ("%s.reverse.mft") % (".noselcys" if Sel == True else "")

      lrepl = AligningMafft(parameters["mafft"], revCurrentInFile, output,
        oFile + ".log", parameters["mafft_params"], parameters["replace"])

      if Sel == True:
        realOutput = ("%s.alg.reverse.mft") % (oFile)
        lrepl = ChangeResiduesSeq(output, realOutput, parameters["in_letter"],
          "U", parameters["replace"])

      lrepl = GettingReverse(parameters["readal"], oFile + '.alg.reverse.mft', \
        oFile + ".alg.reverse.forw.mft", oFile + ".log", parameters["replace"])

    output =  ("%s.alg") % (oFile)
    output += ("%s.forward.mft") % (".noselcys" if Sel == True else "")

    lrepl = AligningMafft(parameters["mafft"], currentInFile, output, \
      oFile + ".log", parameters["mafft_params"], parameters["replace"])

    if Sel == True:
      realOutput = ("%s.alg.forward.mft") % (oFile)
      lrepl = ChangeResiduesSeq(output, realOutput, parameters["in_letter"],
        "U", parameters["replace"])

  ## DiAlign-TX
  if "dialign-tx" in parameters["single_program"]:

    if parameters["both_direction"] == True:
      revInFile = ("%s.seqs.reverse") % (oFile)
      output = ("%s.alg.reverse.dtx") % (oFile)

      lrepl = AligningDialgnTX(parameters["dialign-tx"], revInFile, output,
        oFile + ".log", parameters["dialigntx_params"], parameters["replace"])

      lrepl = GettingReverse(parameters["readal"], oFile + '.alg.reverse.dtx', \
        oFile + ".alg.reverse.forw.dtx", oFile + ".log", parameters["replace"])

    lrepl = AligningDialgnTX(parameters["dialign-tx"], parameters["inFile"], \
      oFile + ".alg.forward.dtx", oFile + ".log", parameters["dialigntx_params"],
      parameters["replace"])

  if lrepl:
    parameters["replace"] = True

  ## If more than one alignment have been generated, generate an meta-alignment
  if len(parameters["single_program"]) > 1 or \
    parameters["both_direction"] == True:

    ## Decide how many directions have been use
    direction = ['forward']
    if parameters["both_direction"] == True:
      direction.append('reverse.forw')

    ## Open output file to store alignments generated previously
    pathFile = open(oFile + ".alg.paths", "w")

    ## Decide which alignments/orientations will be used to generate the
    ## meta-alignment
    for orient in direction:
      for pr in parameters["single_program"]:
        print >> pathFile, ("%s.alg.%s.%s") % (oFile, orient, programs[pr])
    pathFile.close()

    if lrepl:
      parameters["replace"] = True

    ## Generate a meta-alignment using information coming from single alignments
    ## generated with available programs
    MetaAligningMCoffee(parameters["t_coffee"], parameters["inFile"], oFile + \
      ".alg.metalig", oFile + ".log",  oFile + ".alg.paths",
      parameters["mcoffee_params"], parameters["replace"])

    ## Convert meta-alignment into phylip format
    rFile = ("%s.alg.metalig") % (oFile)
    convertAlignment(parameters["trimal"], rFile, rFile, "phylip",
      oFile + ".log", parameters["replace"])

    ## Trim the meta-alignment forcing its selection with information coming
    ## from those alignments used to generate such meta-alignment
    inFile, outFile = ("%s.alg.metalig") % (oFile), ("%s.alg.clean") % (oFile)
    trimmingAlignment(parameters["trimal"], None, outFile, oFile + ".log",
      inFile, oFile + ".alg.paths", parameters["trimal_compare"],
      parameters["trimal_params"], parameters["replace"])

  # Otherwise, trim output aligment using information coming from one alignment
  else:

    ## Select which program has been used to generate the single alignment
    if "muscle" in parameters["single_program"]:
      ext = "alg.forward.msl"
    elif "mafft" in parameters["single_program"]:
      ext = "alg.forward.mft"
    elif "kalign" in parameters["single_program"]:
      ext = "alg.forward.kal"
    elif "dialign-tx" in parameters["single_program"]:
      ext = "alg.forward.dtx"

    ## Convert to phylip format the single alignment generated by any of
    ## available programs
    inFile, outFile = ("%s.%s") % (oFile, ext), ("%s.alg.metalig") % (oFile)
    convertAlignment(parameters["trimal"], inFile, outFile, "phylip",
      oFile + ".log", parameters["replace"])

    ## Trim selected alignment using parameters indicated in the config file
    inFile, outFile = ("%s.alg.metalig") % (oFile), ("%s.alg.clean") % (oFile)
    trimmingAlignment(parameters["trimal"], inFile, outFile, oFile + ".log", "",
      "", "", parameters["trimal_params"], parameters["replace"])

  ## Report how much time has cost the alignment reconstruction part
  final = datetime.datetime.now()
  date  = final.strftime("%H:%M:%S %m/%d/%y")
  lgFile = open(oFile + ".log", "a+")
  print >> lgFile, ("###\n###\n### %s\n### Finishing ALIGNMENT step") % (date)
  total = format_time((final - start).seconds if start else 0)
  print >> lgFile, ("### Total ALIGNMENT step\t%s\n###\n###\n") % (total)
  lgFile.close()

  return oFile + ".alg.clean"

def GettingReverse(bin, inFile, outFile, logFile, replace):
  '''
  Get the reverse alignment using readal for that
  '''

  if lookForFile(outFile) and not replace:
    return
  lgFile = open(logFile, "a+ ") if logFile != "" else None

  ## Construct command-line call and print it onto log file
  cmd = ("%s -in %s -out %s -reverse") % (bin, inFile, outFile)

  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> lgFile, ("###\n### readAl - reverse seqs\n### %s\n### %s\n###") % \
    (date, cmd)
  lgFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(81)
  if proc.wait() != 0:
    print >> sys.stderr, "ERROR: Execution failed: readAl"
    sys.exit(81)

  final = datetime.datetime.now()
  total = format_time((final - start).seconds if start else 0)
  print >> lgFile, ("###\n### Total time\t%s\n###") % (total)
  lgFile.close()

def convertAlignment(bin, inFile, outFile, outFormat, logFile, replace):
  '''
  Convert the input alignmnet into the input format
  '''

  if lookForFile(outFile) and not replace:
    return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None

  ## Construct command-line call and print it onto log file
  cmd = ("%s -in %s -out %s -%s") % (bin, inFile, outFile, outFormat)

  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> lgFile, ("###\n### readAl\n### %s\n### %s###\n") % (date, cmd)
  lgFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(81)

  if proc.wait() != 0:
    print >> sys.stderr, "ERROR: Execution failed: trimAl"
    sys.exit(81)

  final = datetime.datetime.now()
  total = format_time((final - start).seconds if start else 0)
  print >> lgFile, ("###\n### Total time\t%s\n###") % (total)
  lgFile.close()

  return True

def trimmingAlignment(bin, inFile, outFile, logFile, forceSelection,
  compareSelection, compareParameters, parameters, replace):
  '''
  '''

  if lookForFile(outFile) and not replace:
    return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None

  cmd = ""
  ## Construct a customize trimAl command-line call
  if compareSelection:
    cmd = (" -compareset %s") % (compareSelection)
  if forceSelection:
    cmd = ("%s -forceselect %s") % (cmd, forceSelection)
  if compareParameters:
    cmd = ("%s %s") % (cmd, compareParameters)
  if inFile:
    cmd = ("%s -in %s") % (cmd, inFile)
  cmd = ("%s %s -out %s %s") % (bin, cmd, outFile, parameters)

  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> lgFile, ("###\n### trimAl\n### %s\n### %s\n###") % (date, cmd)
  lgFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(82)

  if proc.wait() != 0:
    print >> sys.stderr, "ERROR: Execution failed: trimAl"
    sys.exit(82)

  final = datetime.datetime.now()
  total = format_time((final - start).seconds if start else 0)
  print >> lgFile, ("###\n### Total time\t%s\n###") % (total)
  lgFile.close()

  return True

def AligningMuscle(bin, inFile, outFile, logFile, parameters, replace):
  '''
  '''

  if lookForFile(outFile) and not replace:
    return False;
  lgFile = open(logFile, "a+ ") if logFile != "" else None
  lgFile.flush()

  cmd = ("%s %s -in %s -out %s") % (bin, parameters, inFile, outFile)
  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> lgFile, ("###\n### MUSCLE\n### %s\n### %s\n###") % (date, cmd)
  lgFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(86)

  if proc.wait() != 0:
    print >> sys.stderr, "ERROR: Execution failed: Muscle"
    sys.exit(86)

  if not CheckGeneratedAlignment(inFile, outFile):
    print >> sys.stderr, "ERROR: Execution failed: Muscle - Alignment Check"
    sp.call(("rm -f %s") % (outFile), shell = True)
    sys.exit(86)

  final = datetime.datetime.now()
  total = format_time((final - start).seconds if start else 0)
  print >> lgFile, ("###\n### Total time\t%s\n###") % (total)
  lgFile.close()

  return True

def AligningKAlign(bin, inFile, outFile, logFile, parameters, replace):
  '''
  '''

  if lookForFile(outFile) and not replace:
    return False;
  lgFile = open(logFile, "a+ ") if logFile != "" else None

  cmd = ("%s %s -in %s -out %s") % (bin, parameters, inFile, outFile)
  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> lgFile, ("###\n### KAlign\n### %s\n### %s\n###") % (date, cmd)
  lgFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(89)

  if proc.wait() != 0:
    print >> sys.stderr, "ERROR: Execution failed: KAlign"
    sys.exit(89)

  ## Check whether output alignment have the same number of sequences/residues
  ## than the input unaligned sequences file
  if not CheckGeneratedAlignment(inFile, outFile):
    print >> sys.stderr, "ERROR: Execution failed: KAlign - Alignment Check"
    sp.call(("rm -f %s") % (outFile), shell = True)
    sys.exit(89)

  final = datetime.datetime.now()
  total = format_time((final - start).seconds if start else 0)
  print >> lgFile, ("###\n### Total time\t%s\n###") % (total)
  lgFile.close()

  return True

def AligningMafft(bin, inFile, outFile, logFile, parameters, replace):
  '''
  '''

  if lookForFile(outFile) and not replace:
    return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None

  cmd = ("%s %s %s > %s") % (bin, parameters, inFile, outFile)
  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> lgFile, ("###\n### MAFFT\n### %s\n### %s\n###") % (date, cmd)
  lgFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(87)

  if proc.wait() != 0:
    print >> sys.stderr, "ERROR: Execution failed: MAFFT"
    sys.exit(87)

  ## Check whether output alignment have the same number of sequences/residues
  ## than the input unaligned sequences file
  if not CheckGeneratedAlignment(inFile, outFile):
    print >> sys.stderr, "ERROR: Execution failed: MAFFT - Alignment Check"
    sp.call(("rm -f %s") % (outFile), shell = True)
    sys.exit(87)

  final = datetime.datetime.now()
  total = format_time((final - start).seconds if start else 0)
  print >> lgFile, ("###\n### Total time\t%s\n###") % (total)
  lgFile.close()

  return True

def AligningDialgnTX(bin, inFile, outFile, logFile, parameters, replace):
  '''
  '''

  if lookForFile(outFile) and not replace:
    return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None

  cmd = ("%s %s %s %s") % (bin, parameters, inFile, outFile)
  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> lgFile, ("###\n### DiAlign-TX\n### %s\n### %s\n###") % (date, cmd)
  lgFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(88)

  if proc.wait() != 0:
    print >> sys.stderr, "ERROR: Execution failed: DiAlign-TX"
    sys.exit(88)

  if not CheckGeneratedAlignment(inFile, outFile):
    print >> sys.stderr, "ERROR: Execution failed: DiAlign-TX - Alignment Check"
    sp.call(("rm -f %s") % (outFile), shell = True)
    sys.exit(88)

  final = datetime.datetime.now()
  total = format_time((final - start).seconds if start else 0)
  print >> lgFile, ("###\n### Total time\t%s\n###") % (total)
  lgFile.close()

  return True

def MetaAligningMCoffee(bin, inFile, outFile, logFile, pathFile, pars, replace):
  '''
  '''

  if lookForFile(outFile) and not replace:
    return False
  lgFile = open(logFile, "a+ ") if logFile != "" else None

  ## Define some environment variables for t_coffee
  # user = os.getenv("LOGNAME")
  # try:
  #   proc = sp.Popen(("mkdir -p -m0777 /tmp/tcoffee/%s;") % (user), shell = True)
  #   os.putenv("LOCKDIR_4_TCOFFEE", ("/tmp/tcoffee/%s") % (user))
  #   os.putenv("TMP_4_TCOFFEE", ("/tmp/tcoffee/%s") % (user))
  # except OSError, e:
  #   pass

  ## Get all alignments which will be used to construct the MetaAlignment
  files = " ".join([line.strip() for line in open(pathFile, "rU")])
  files = (" -aln %s") % (files)

  cmd = ("%s %s -outfile %s %s %s") % (bin, inFile, outFile, pars, files)
  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> lgFile, ("###\n### M-COFFEE\n### %s\n### %s\n###") % (date, cmd)
  # localVars =  ("###\n### TMP_4_TCOFFEE=/tmp/tcoffee/%s\n") % (user)
  # localVars += ("### LOCKDIR_4_TCOFFEE=/tmp/tcoffee/%s") % (user)
  # print >> lgFile, ("###\n### M-COFFEE\n### %s\n%s\n### %s\n###") % (date,
  #  localVars, cmd)
  lgFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = lgFile, stdout = lgFile)
  except OSError, e:
    print >> sys.stderr, "ERROR: Execution failed: " + str(e)
    sys.exit(90)

  if proc.wait() != 0:
    print >> sys.stderr, "ERROR: Execution failed: M-Coffee"
    sys.exit(90)

  if not CheckGeneratedAlignment(inFile, outFile):
    print >> sys.stderr, "ERROR: Execution failed: M-Coffee - Alignment Check"
    sp.call(("rm -f %s") % (outFile), shell = True)
    sys.exit(90)

  final = datetime.datetime.now()
  total = format_time((final - start).seconds if start else 0)
  print >> lgFile, ("###\n### Total time\t%s\n###") % (total)
  lgFile.close()

  ## Remove any temp file
  sp.call("rm -f " + os.path.split(inFile)[1][:10] + "*.dnd", shell = True)

  return True



def ChangeResiduesSeq(inFile, outFile, inLetter, outLetter, replace):
  '''
  '''

  if lookForFile(outFile) and not replace:
    return False;

  oFile = open(outFile, "w")
  for line in [line for line in map(strip, open(inFile, "rU")) if line]:
    if line[0] == ">":
      print >> oFile, line
    else:
      print >> oFile, line.replace(inLetter, outLetter)
  oFile.close()

  return True

def CheckGeneratedAlignment(inSeqsFile, outAligFile):
  '''
  '''

  ## We introduce a delay to ensure data is already written in the disk.
  ## With high-computing facilities, sometimes there are some problems of
  ## writing to disk the already computed results
  sleep(5)

  ## Get input sequences before aligning them
  inSeqs, currentId = {}, None
  for line in [line for line in map(strip, open(inSeqsFile, "rU")) if line]:
    if line[0] == ">":
      currentId = line[1:]
      inSeqs.setdefault(currentId, "")
    elif currentId:
      inSeqs[currentId] += line.replace("-", "")

  ## Get sequences after aligning them to check whether the same number of
  ## sequences/residues are present in the output alignment.
  outSeqs, currentId = {}, None
  for line in [line for line in map(strip, open(outAligFile, "rU")) if line]:
    if line[0] == ">":
      currentId = line[1:]
      outSeqs.setdefault(currentId, "")
    elif currentId:
      outSeqs[currentId] += line.replace("-", "")
  ## Check both files have exactly the same number of sequences
  if len(inSeqs) != len(outSeqs):
    return False

  ## ... and each sequences have exactly the same residues
  for seqId in inSeqs:
    if not seqId in outSeqs or inSeqs[seqId] != outSeqs[seqId]:
      return False

  ## If everything is OK, inform about it
  return True
