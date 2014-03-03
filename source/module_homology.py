#!/usr/bin/python

import os
import sys
import tempfile
import datetime
import subprocess as sp
from Bio import SeqIO
from hashlib import md5
from string import strip, capitalize, ljust
from module_utils import lookForDirectory, lookForFile, splitSequence, \
  parseComments, format_time, sort_blast_hits, sort_hmmer_hits

def homology(parameters):

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
  print >> logFile, ("###\n###\tSTEP\tHomology\tSTART\t%s\n###") % (date)
  logFile.flush()

  ## Get which tool will be used to perform the homology search. Check such tool
  ## is listed among the available binaries
  if not "homology" in parameters:
    sys.exit("ERROR: Check your configuration file. There is not tool set for "
      + "the homology search")

  if not parameters["homology"][0] in parameters:
    sys.exit("ERROR: Check your configuration file. This tool '%s' is not among"
      + " available methods")

  ## If the homology search will use any program from the BLAST package, check
  ## whether the TARGET SEQUENCES file has been already formatted.
  if parameters["homology"][0] in ["legacy_blast", "blast+"]:

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
          db_file = parameters["db_file"]
          sys.exit(("ERROR: Check your input TARGET SEQUENCES file '%s' has "
            + "been formated using 'formatdb'/'makeblastdb'") % (db_file))

    ## If the homology search step should be perfomed using BLAST, call the
    ## appropiate function
    blast(parameters, logFile)

  elif parameters["homology"][0] in ["phmmer", "jackhmmer", "hmmer_search"]:
    hmmer(parameters, logFile)

  ## Filter homology search data. A dictionary containing selected sequences,
  ## including the sequences themselves
  selected_sequences = filter_results(parameters, logFile)

  ## Generate a MD5 file containing selected sequences for the current run.
  ## MD5s are used to recompute the same phylogenetic tree starting from other
  ## seqs - with identical similarity search results - in the set of homologs
  outFile = ("%s.seqs.md5") % (oFile)

  ## Check whether the file already exists or not.
  if not lookForFile(outFile) or parameters["replace"]:
    parameters["replace"] = True

    seqs_md5 = md5("".join(sorted(selected_sequences.keys()))).hexdigest()
    print >> open(outFile, "w"), ("%s\t%s") % (parameters["prefix"], seqs_md5)

  ## Generate a file containing the selected sequences after performing the
  ## homology search and filtering its output according to a set of parameters.
  outFile = ("%s.seqs") % (oFile)

  ## Check whether the file already exists or not.
  if not lookForFile(outFile) or parameters["replace"]:
    parameters["replace"] = True

    output_file = open(outFile, "w")
    for seqId in sorted(selected_sequences):
      print >> output_file, (">%s\n%s") % (seqId, selected_sequences[seqId][1])
    output_file.close()

  ## Print how much time was needed to perform the whole homology search step
  final = datetime.datetime.now()
  date  = final.strftime("%H:%M:%S %m/%d/%y")

  print >> logFile, ("###\n###\tSTEP\tHomology\tEND\t%s") % (date)
  total = format_time((final - start).seconds if start else 0)
  print >> logFile, ("###\tTOTAL Time\tHomology\t%s\n###") % (total)
  logFile.close()

  ## Update the input file parameter and return the dictionary containing all
  ## parameters. Those parameters may be used in other steps
  parameters["in_file"] = outFile

  ## Before returning to the main program, get back to the original working
  ## directory
  os.chdir(current_directory)

  return parameters

def blast(parameters, logFile):
  '''
  Perform the homology search using the different BLAST package programs. This
  module offers retrocompatibility to legacy blast.
  '''

  ## Get output folder/generic filename
  oFile = os.path.join(parameters["out_directory"], parameters["prefix"])

  ## Get output file name and check whether has been previously generated or
  ## not. It will also affect whether the variable REPLACE is set or not
  outFile = ("%s.homology.blast.out") % (oFile)

  ## If output file exist and it is not set to replace it, just go back to the
  ## main function. Otherwise, set the replace parameter to TRUE in other to
  ## replace any already generated file downstream
  if lookForFile(outFile) and not parameters["replace"]:
    return
  parameters["replace"] = True

  ## Generate command-line depending on which BLAST package is being used.
  if parameters["homology"][0] == "legacy_blast":
    binary = parameters["legacy_blast"][0]
    params = parameters[binary +"_params"]
    cmd = ("%s %s -e %s -d %s -i %s -o %s") % (parameters[binary], params, \
      str(parameters["e_value"]), parameters["db_file"], parameters["in_file"],\
      outFile)

  elif parameters["homology"][0] == "blast+":
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

  ## Remove any error file generated during the legacy_blast execution - We try
  ## to delete this file only if it is empty
  if not lookForFile("error.log"):
    sp.call(("rm -f error.log"), shell = True)

def hmmer(parameters, logFile):
  '''
  Perform the homology search using three different approximations implemented
  in the HMMER package.
  '''

  ## Get output folder/generic filename
  oFile = os.path.join(parameters["out_directory"], parameters["prefix"])

  ## Get output file name and check whether has been previously generated or
  ## not. It will also affect whether the variable REPLACE is set or not
  outFile = ("%s.homology.hmmer.out") % (oFile)

  ## If output file exist and it is not set to replace it, just go back to the
  ## main function. Otherwise, set the replace parameter to TRUE in other to
  ## replace any already generated file downstream
  if lookForFile(outFile) and not parameters["replace"]:
    return
  parameters["replace"] = True

  ## If we are ask to perform a HMM search using a Multiple Sequence Alignment
  ## as input rather than a single sequence, we need first to construct a HMM
  ## to perfom the search
  if parameters["homology"][0] == "hmmsearch" :
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
    hmmFile = ("%s.homology.hmmer.hmm") % (oFile)

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
  binary = parameters["homology"][0]
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

def filter_results(parameters, logFile):
  '''
  Filter Homology search results taking into account which package was used to
  perform the search. Depending on the package only e-values (HMMER) or e-value
  plus coverage -ratio of aligned region between query and target sequences vs.
  query sequence lenght- (BLAST) are used.
  '''

  ## Get output folder/generic filename
  oFile = os.path.join(parameters["out_directory"], parameters["prefix"])

  ## Get tag for the input/output file. It will depend on which method has been
  ## used to perform the homology seach
  tag = "hmmer" if parameters["homology"][0] in ["phmmer", "jackhmmer", \
    "hmmer_search"] else "blast" if parameters["homology"][0] in \
    ["legacy_blast",  "blast+"] else ""

  ## Get input file
  inFile = ("%s.homology.%s.out") % (oFile, tag)
  ## If input file doesn't exist, just go back to the main function
  if not lookForFile(inFile):
    sys.exit(("ERROR: Check previously generated file '%s'") % (inFile))

  ## Get input file name and check whether has been previously generated or
  ## not. It will also affect whether the variable REPLACE is set or not
  outFile = ("%s.homology.%s.filter") % (oFile, tag)

  ## If output file exist and it is not set to replace it, just go back to the
  ## main function. Otherwise, set the replace parameter to TRUE in other to
  ## replace any already generated file downstream
  if lookForFile(outFile) and not parameters["replace"]:
    return
  parameters["replace"] = True

  input_lines, target_sequences, query_line = [], set(), None
  for line in open(inFile, "rU"):
    ## Parse line
    parsed_line = [element for element in parseComments([e for e in map(strip, \
      line.split()) if e]) if element]

    ## Discard empty lines or those starting by "#"
    if not parsed_line:
      continue

    ## Detect the target sequence which is placed at different columns depending
    ## whether it is blast or hmmer package which generated the output
    target = parsed_line[0] if tag == "hmmer" else parsed_line[1]

    ## We also include the query sequence, it is only important for the BLAST-
    ## based search
    query = parsed_line[2] if tag == "hmmer" else parsed_line[0]

    ## Discard previously found target sequences
    if target in target_sequences:
      continue
    input_lines.append(parsed_line)
    target_sequences|= set([target, query])

    if target == query and not query_line:
      query_line = parsed_line

  sequences = read_database(parameters["db_file"], target_sequences)

  ## Depending on how the search was performed, we will filter-out data
  ## by e-values and coverage (BLAST only) or not
  e_value = float(parameters["e_value"])
  coverage = float(parameters["coverage"])
  hits = -1 if not "hits" in parameters or parameters["hits"] == "no_limit" \
    else int(parameters["hits"])

  accepted_lines, accepted_targets = [], set()
  for line in input_lines:
    ## If the current target has been already found, move to next hit
    if (line[0] if tag == "hmmer" else line[1]) in accepted_targets:
      continue
    ## Depending on the package, filter just by two e-values (sequence and best
    ## found domain) or by sequence e-value + coverage between sequences
    if tag == "hmmer":
      if float(line[4]) > e_value or float(line[7]) > e_value:
        continue
    elif tag == "blast":
      covTarget = ((int(line[7]) - int(line[6]))+1)/float(sequences[line[0]][0])
      if covTarget < coverage or float(line[-2]) > e_value:
        continue

    ## Store current line and target sequence
    accepted_lines.append(line)
    accepted_targets.add(line[0] if tag == "hmmer" else line[1])

  ## Sort by e-values (and bit-score for BLAST only) accepted lines
  accepted_lines.sort(sort_blast_hits if tag == "blast" else sort_hmmer_hits)

  if hits != -1 and len(accepted_lines) >= hits:
    if not query in accepted_targets and query_line:
      accepted_lines = query_line + accepted_lines[:hits-1]
    else:
      accepted_lines = accepted_lines[:hits]

  ## Get selected sequences. It will be used to produce MD5s key as well as to
  ## generate the sequences FASTA file
  selected_sequences = {}
  for line in accepted_lines:
    sequence_id = line[0] if tag == "hmmer" else line[1]
    selected_sequences.setdefault(sequence_id, sequences[sequence_id])

  out = ["\t".join(map(lambda x: str(x).ljust(6), l)) for l in accepted_lines]
  print >> open(outFile, "w"), "\n".join(out)

  return selected_sequences

def read_database(input_db_file, sequences):
  '''
  Read input TARGET Sequences database returning those sequences which may be
  potentially selected during the filtering step
  '''
  output = {}
  for record in SeqIO.parse(input_db_file, "fasta"):
    if not record.id in sequences:
      continue
    seq = str(record.seq) if record.seq[-1] != "*" else record.seq[:-1]
    output.setdefault(record.id, (len(seq), splitSequence(seq)))
  return output
