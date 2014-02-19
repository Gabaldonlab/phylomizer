#!/usr/bin/python
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
import os
import hashlib
import tempfile
import datetime
import subprocess as sp
from Bio import SeqIO
from string import strip, capitalize
from module_utils import *

def homology(parameters):

  ## Get input folder/filename
  iFolder, iFile = os.path.split(parameters["in_file"])

  ## Get output folder/generic filename
  oFile = os.path.join(parameters["out_directory"], iFile.split(".")[0])

  ## Set output filename and log file
  logFile = open(oFile + ".log", "w" if parameters["replace"] else "a+")

  start = datetime.datetime.now()
  date = start.strftime("%H:%M:%S %m/%d/%y")
  print >> logFile, ("###\n###\tSTEP\tHomology\tSTART\t%s\n###") % (date)
  logFile.close()

  ## Check which mode will be used for the homology search
  modes = set(parameters.keys()) & set(["hmmer", "legacy_blast", "blast+"])
  if len(modes) > 1:
    packs = ",".join(sorted(map(capitalize, modes)))
    sys.exit(("ERROR: Check your configuration file. There is more than one "
      + "package configured for performing the homology search '%s'") % (packs))
  elif len(modes) < 1:
    sys.exit("ERROR: Check your configuration file. There is not package select"
      + "ed for performing the homology search")

  ## If the homology search will use any program from the BLAST package, check
  ## whether the TARGET SEQUENCES file has been already formatted.
  if "legacy_blast" in parameters or "blast+" in parameters:

    ## Get database sequence type - p: protein or n:nucleotide
    dt = "p" if parameters["residue_datatype"].startswith("prot") else "n"

    ## Check if BLAST DB associated files already exist or not
    for extension in ["hr", "in", "sq"]:
      filename = ("%s.%s%s") % (parameters["db_file"], dt, extension)
      ## If the input file doesn't exist check whether input database has been
      ## split into different volumes
      if not lookForFile(filename):
        alternative = ("%s.00.%s%s") % (parameters["db_file"], dt, extension)
        if not lookForFile(alternative):
          sys.exit(("ERROR: Check your input TARGET SEQUENCES file '%s' has "
            + "been formated using 'formatdb'") % (parameters["db_file"]))

    ## If the homology search step should be perfomed using BLAST, call the
    ## appropiate function
    blast(parameters)

  elif "hmmer" in parameters:
    hmmer(parameters)

  final = datetime.datetime.now()
  date  = final.strftime("%H:%M:%S %m/%d/%y")
  logFile = open(oFile + ".log", "a+")
  print >> logFile, ("###\n###\tSTEP\tHomology\tEND\t%s") % (date)
  total = format_time((final - start).seconds if start else 0)
  print >> logFile, ("###\tTOTAL Time\tHomology\t%s\n###") % (total)
  logFile.close()
#~
#~
#~
  #~ db, mappedIDs = None, None
  #~ outFile = ("%s.blast.filter") % (oFile)
#~
  #~ ## If BLAST filtered results, read BLAST output and filter it according to
  #~ ## config options
  #~ if not lookForFile(outFile) or parameters["replace"]:
    #~ parameters["replace"] = True
#~
    #~ if not start:
      #~ start = datetime.datetime.now()
#~
    #~ output = [line for line in open(("%s.blast.out") % (oFile), "rU")]
    #~ db, mappedIDs = getMappedDB(parameters["BlastDB"])
#~
    #~ proteinID = iFile.split(".")[0]
    #~ lenQuery = len(lookForSequence(proteinID, db[mappedIDs[proteinID]]))
#~
    #~ filteredResults = filterBLAST(proteinID, output, lenQuery,
      #~ parameters["coverage"], parameters['hits'], db,  mappedIDs)
#~
    #~ generateBlastOutput(outFile, filteredResults)
#~
  #~ ## If BLAST filtered results already exist, read information.
  #~ elif lookForFile(outFile):
    #~ filteredResults = [map(strip, l.split("\t")) for l in open(outFile, "rU")]
#~
  #~ outFile = ("%s.seqs") % (oFile)
  #~ ## Generate output sequences files
  #~ if not lookForFile(outFile) or parameters["replace"]:
    #~ parameters["replace"] = True
#~
    #~ if not start:
      #~ start = datetime.datetime.now()
#~
    #~ ## If Blast database has not been loaded in a previous step, do it
    #~ if not db or not mappedIDs:
      #~ db, mappedIDs = getMappedDB(parameters["BlastDB"])
#~
    #~ generateSequences(outFile, filteredResults, db, mappedIDs)
#~
  #~ outFile = ("%s.seqs.md5") % (oFile)
  #~ ## Generate MD5 file if already does not exist or any of previous files
  #~ ## has been modified
  #~ if not lookForFile(outFile) or parameters["replace"]:
    #~ parameters["replace"] = True
#~
    #~ if not start:
      #~ start = datetime.datetime.now()
#~
    #~ key = iFile.split(".")[0]
    #~ generateMD5(outFile, parameters["outDirec"], key, filteredResults)
#~
  #~ final = datetime.datetime.now()
  #~ date  = final.strftime("%H:%M:%S %m/%d/%y")
  #~ print >> logFile, ("### %s\n### Finishing BLAST step") % (date)
  #~ total = format_time((final - start).seconds if start else 0)
  #~ print >> logFile, ("### Total BLAST step\t%s") % (total)
  #~ logFile.close()
#~
  #~ return ("%s.seqs") % (oFile)

def blast(parameters):

  ## Get input folder/filename
  iFolder, iFile = os.path.split(parameters["in_file"])

  ## Get output folder/generic filename
  oFile = os.path.join(parameters["out_directory"], iFile.split(".")[0])

  ## Set output filename and log file
  logFile = open(oFile + ".log", "a+")

  ## Get output file name and check whether has been previously generated or
  ## not. It will also affect whether the variable REPLACE is set or not
  outFile = ("%s.homology.out") % (oFile)

  ## If BLAST output does not exist, generate it and, therefore, replace any
  ## file generated downstream
  if not lookForFile(outFile) or parameters["replace"]:
    parameters["replace"] = True

  ## Generate command-line depending on which BLAST package is being used.
  if "legacy_blast" in parameters:
    binary = parameters["legacy_blast"][0]
    params = parameters[binary +"_params"]
    cmd = ("%s %s -e %s -d %s -i %s -o %s") % (parameters[binary], params, \
      str(parameters["e_value"]), parameters["db_file"], parameters["in_file"],\
      outFile)

  elif "blast+" in parameters:
    binary = parameters["blast+"][0]
    params = parameters[binary +"_params"]
    cmd = ("%s %s -evalue %s -db %s -query %s -out %s") % (parameters[binary], \
      params, str(parameters["e_value"]), parameters["db_file"], \
      parameters["in_file"], outFile)

  print >> logFile, ("###\n###\tCommand-line\t%s\n###\n") % (cmd)
  logFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = logFile)
  except OSError, e:
    sys.exit("ERROR: Execution failed: " + str(e))

  if proc.wait() != 0:
    sys.exit(("ERROR: Execution failed: '%s'") % (parameters[binary]))

  logFile.close()

def hmmer(parameters):

  ## Get input folder/filename
  iFolder, iFile = os.path.split(parameters["in_file"])

  ## Get output folder/generic filename
  oFile = os.path.join(parameters["out_directory"], iFile.split(".")[0])

  ## Set output filename and log file
  logFile = open(oFile + ".log", "a+")

  ## Get output file name and check whether has been previously generated or
  ## not. It will also affect whether the variable REPLACE is set or not
  outFile = ("%s.homology.out") % (oFile)

  ## If output file does not exist, generate it and, therefore, replace any
  ## file generated downstream
  if not lookForFile(outFile) or parameters["replace"]:
    parameters["replace"] = True

  ## If we are ask to perform a HMM search using a Multiple Sequence Alignment
  ## as input rather than a single sequence, we need first to construct a HMM
  ## to perfom the search
  if parameters["hmmer"][0] == "hmmsearch":
    if not "readal" in parameters or not "hmmbuild" in parameters:
      sys.exit(("ERROR: Check your CONFIG file to search whether 'readAl' and "
        + "'hmmbuild' are available"))

    ## Create a temporary FASTA file which will be used as input for HMMBuild
    TEMPFILE = tempfile.NamedTemporaryFile()
    cmd = ("%s -in %s -out %s -fasta") % (parameters["readal"], \
      parameters["in_file"], TEMPFILE.name)
    sp.call(cmd, shell = True)
    TEMPFILE.flush()

    ## Generate the profile
    ## Set the current residues type to amino-acids if search is performed using
    ## proteins, otherwise, allow the program to guess it
    dt = "--amino" if parameters["residue_datatype"].startswith("prot") else ""
    hmmFile = ("%s.homology.hmm") % (oFile)

    cmd = ("%s --informat afa %s %s %s") % (parameters["hmmbuild"], dt, hmmFile,
      TEMPFILE.name)

    print >> logFile, ("###\n###\tCommand-line\t%s\n###\n") % (cmd)
    logFile.flush()

    try:
      proc = sp.Popen(cmd, shell = True, stderr = logFile, stdout = logFile)
    except OSError, e:
      sys.exit("ERROR: Execution failed: " + str(e))

    if proc.wait() != 0:
      sys.exit(("ERROR: Execution failed: '%s'") % ("hmmbuild"))

    ## We update the input file for performing the HMM-based homology search
    parameters["in_file"] = hmmFile
    TEMPFILE.close()

  ## Generate command-line depending on HMMER specific program and parameters
  binary = parameters["hmmer"][0]
  params = parameters["hmmer_params"]

  cmd = ("%s %s -E %s --tblout %s %s %s") % (parameters[binary], params, \
    str(parameters["e_value"]), outFile, parameters["in_file"], \
    parameters["db_file"])

  print >> logFile, ("###\n###\tCommand-line\t%s\n###\n") % (cmd)
  logFile.flush()

  try:
    proc = sp.Popen(cmd, shell = True, stderr = logFile, stdout = logFile)
  except OSError, e:
    sys.exit("ERROR: Execution failed: " + str(e))

  if proc.wait() != 0:
    sys.exit(("ERROR: Execution failed: '%s'") % (parameters[binary]))

  logFile.close()

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def filterBLAST(seed_proteinID, blastResult, lenQuery, coverage, numberHits, \
  db, mappedIDs):
  '''
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  upperLimit, lowerLimit = (3 * lenQuery) + 1,(lenQuery / 3) - 1
  lines, seedLine, acceptedHits = [], [], []
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  for line in blastResult:
    fields = line.strip().split() if line != "" else None
    protid = fields[1].split("|")[0]

    if fields and protid == seed_proteinID and not seedLine:
      seedLine = [fields]

    if fields and protid not in acceptedHits:
      protid = fields[1].split("|")[0]
      seq = len(lookForSequence(protid, db[mappedIDs[protid]]))
      covQueryHit = ((int(fields[7]) - int(fields[6])) + 1) / float(lenQuery)

      if covQueryHit > float(coverage) and seq < upperLimit and seq > lowerLimit:
        acceptedHits.append(protid)
        lines.append(fields)
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  lines.sort(sort_hits)
  if not seed_proteinID in acceptedHits:
    if len(lines) >= int(numberHits):
      lines = seedLine + lines[:int(numberHits) - 1]
    else:
      lines = seedLine + lines
  return lines[:int(numberHits)]
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def generateBlastOutput(outFile, outLines):
  '''
  Given an output file and a certain lines, print these lines in the given file
  '''
  oFile = open(outFile, "w")
  for line in outLines: print >> oFile, "\t". join(line)
  oFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def generateSequences(outFile, outLines, db, mappedIDs):
  '''
  '''
  oFile = open(outFile, "w")
  for hit in outLines:
    proteinID = hit[1].split("|")[0]
    seq = lookForSequence(proteinID, db[mappedIDs[proteinID]])
    print >> oFile, ">" + proteinID
    for line in split_len(seq, 60): print >> oFile, line
  oFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def generateMD5(outFile, outFolder, ID, outLines):
  '''
  '''
  oFile = open(outFile, "w")
  hits = ''.join(sorted([ hit[1].split("|")[0] for hit in outLines]))
  print >> oFile, ID + "\t" + hashlib.md5(hits).hexdigest()
  oFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def getMappedDB(BlastDB):
  '''
  Returns a dictionnary with Seq objetcs contained in a multiFASTA file,
  under "longest" or "others" key depending if the sequence is the longest
  isoform of a given gene or not
  '''
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  MultiFASTA, db, record, prot, gen = open(BlastDB, "rU"), {}, {}, "", ""
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  for line in MultiFASTA.readlines():
    if line[0] == ">":
      if record != {}:
        # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
        if db.has_key(record['gen']):
          if len(db[record['gen']]['longest']['seq']) < len(record['seq']):
            db[record['gen']]['others'].append(db[record['gen']]['longest'])
            db[record['gen']]['longest'] = record
          else:
            db[record['gen']]['others'].append(record)
        else:
          db[record['gen']] = { 'longest': record, 'others': []}
      # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
      prot = line[1:].strip().split('|')[0]
      try: gen = line[1:].strip().split('|')[1] + '-' + prot[0:3]
      except: gen = prot  + '-' + prot[0:3]
      gen = gen.replace(" ", "_") if len(gen.split()) > 1 else gen
      # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
      record = {'id': prot, 'gen': gen, 'seq': "" }
      # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
    else: record['seq'] += line.strip().replace("*", "")
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  if db.has_key(record['gen']):
    if len(db[record['gen']]['longest']['seq']) < len(record['seq']):
      db[record['gen']]['others'].append(db[record['gen']]['longest'])
      db[record['gen']]['longest'] = record
    else:
      db[record['gen']]['others'].append(record)
  else:
    db[record['gen']] = { 'longest': record, 'others': []}
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  mappedIDs = {}
  for gen in db:
    mappedIDs[db[gen]['longest']['id']] = gen
    for protein in db[gen]['others']:
      mappedIDs[protein['id']] = gen
  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

  # ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
  MultiFASTA.close()
  return db, mappedIDs
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
def lookForSequence(key, entryDB):
  '''
  '''
  if entryDB['longest']['id'] == key:
    return entryDB['longest']['seq']
  for member in entryDB['others']:
    if member['id'] == key: return member['seq']
  return ''
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
