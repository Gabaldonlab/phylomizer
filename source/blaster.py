#!/usr/bin/python
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
import hashlib, os, subprocess as sp
from string import strip
from Bio import SeqIO
from utils import *
import datetime

def blast(parameters):

  ## Get input folder/filename
  iFolder, iFile = os.path.split(parameters["inFile"])
  ## Get output folder/generic filename
  oFile = os.path.join(parameters["outDirec"], iFile.split(".")[0])
  ## Set output filename and log file
  logFile = open(oFile + ".log", "a+")
  ## Set variables to monitor time execution
  start = None

  ## Check if BLAST db associated files already exist or not
  if not lookForFile(parameters["BlastDB"] + ".phr") or \
     not lookForFile(parameters["BlastDB"] + ".pin") or \
     not lookForFile(parameters["BlastDB"] + ".psq"):

    if not lookForFile(parameters["BlastDB"] + ".00.phr") or \
       not lookForFile(parameters["BlastDB"] + ".00.pin") or \
       not lookForFile(parameters["BlastDB"] + ".00.psq"):
      sys.exit("ERROR: You should format your BlastDB")

  outFile = ("%s.blast.out") % (oFile)
  ## If BLAST output does not exist, generate it
  if not lookForFile(outFile) or parameters["replace"]:
    parameters["replace"] = True

    if not lookForFile(parameters["inFile"]):
      sys.exit(("ERROR: Execution failed: Input FASTA File not found '%s'") \
        % (parameters["inFile"]))

    ## Generate command-line
    cmd = ("%s %s -e %s -d %s -i %s -o %s") % (parameters["blastpgp"], \
      parameters["blast_params"], str(parameters["e_value"]),
      parameters["BlastDB"], parameters["inFile"], outFile)

    start = datetime.datetime.now()
    date = start.strftime("%H:%M:%S %m/%d/%y")
    print >> logFile, ("### %s\n### %s") % (date, cmd)

    try:
      proc = sp.Popen(cmd, shell = True, stderr = logFile)
    except OSError, e:
      sys.exit("ERROR: Execution failed: " + str(e))
    if proc.wait() != 0:
      sys.exit(("ERROR: Execution failed: '%s'") % (parameters["blastpgp"]))

  db, mappedIDs = None, None
  outFile = ("%s.blast.filter") % (oFile)

  ## If BLAST filtered results, read BLAST output and filter it according to
  ## config options
  if not lookForFile(outFile) or parameters["replace"]:
    parameters["replace"] = True

    if not start:
      start = datetime.datetime.now()

    output = [line for line in open(("%s.blast.out") % (oFile), "rU")]
    db, mappedIDs = getMappedDB(parameters["BlastDB"])

    proteinID = iFile.split(".")[0]
    lenQuery = len(lookForSequence(proteinID, db[mappedIDs[proteinID]]))

    filteredResults = filterBLAST(proteinID, output, lenQuery,
      parameters["coverage"], parameters['hits'], db,  mappedIDs)

    generateBlastOutput(outFile, filteredResults)

  ## If BLAST filtered results already exist, read information.
  elif lookForFile(outFile):
    filteredResults = [map(strip, l.split("\t")) for l in open(outFile, "rU")]

  outFile = ("%s.seqs") % (oFile)
  ## Generate output sequences files
  if not lookForFile(outFile) or parameters["replace"]:
    parameters["replace"] = True

    if not start:
      start = datetime.datetime.now()

    ## If Blast database has not been loaded in a previous step, do it
    if not db or not mappedIDs:
      db, mappedIDs = getMappedDB(parameters["BlastDB"])

    generateSequences(outFile, filteredResults, db, mappedIDs)

  outFile = ("%s.seqs.md5") % (oFile)
  ## Generate MD5 file if already does not exist or any of previous files
  ## has been modified
  if not lookForFile(outFile) or parameters["replace"]:
    parameters["replace"] = True

    if not start:
      start = datetime.datetime.now()

    key = iFile.split(".")[0]
    generateMD5(outFile, parameters["outDirec"], key, filteredResults)

  final = datetime.datetime.now()
  date  = final.strftime("%H:%M:%S %m/%d/%y")
  print >> logFile, ("### %s\n### Finishing BLAST step") % (date)
  total = format_time((final - start).seconds if start else 0)
  print >> logFile, ("### Total BLAST step\t%s") % (total)
  logFile.close()

  return ("%s.seqs") % (oFile)

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


