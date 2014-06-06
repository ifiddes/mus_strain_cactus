"""
convenience library for assisting filters.
"""
import os


class Sequence(object):
  """ Represents a sequence of DNA.
  """
  def __init__(self, name, sequence):
    self.name = name  # chromosome or scaffold name
    self._sequence = sequence  # ACGTs
  def setSequence(self, seq):
    self._sequence = seq
  def getSequence(self):
    return self._sequence
  def getLength(self):
    return len(self._sequence)


class PslRow(object):
  """ Represents a single row in a PSL file.
  http://genome.ucsc.edu/FAQ/FAQformat.html#format2
  """
  def __init__(self, line):
    data = line.split()
    assert(len(data) == 20)
    self.matches = int(data[0])
    self.misMatches = int(data[1])
    self.nCount = int(data[2])
    self.qNumInsert = int(data[3])
    self.qBaseInsert = int(data[4])
    self.tNumInsert = int(data[5])
    self.tBaseInsert = int(data[6])
    self.strand = data[7]
    self.qName = data[8]
    self.qSize = int(data[9])
    self.qStart = int(data[10])
    self.qEnd = int(data[11])
    self.tName = data[12]
    self.tSize = int(data[13])
    self.tStart = int(data[14])
    self.tEnd = int(data[15])
    self.blockCount = int(data[16])
    self.blockSizes = int(data[17])
    self.qStarts = int(data[18])
    self.tStarts = int(data[19])


def initializeArguments(parser):
  """ given an argparse ArgumentParser object, add in the default arguments.
  """
  parser.add_argument('--refGenome')
  parser.add_argument('--genome')
  parser.add_argument('--geneCheckBed')
  parser.add_argument('--geneCheckBedDetails')
  parser.add_argument('--alignment')
  parser.add_argument('--sequence')
  parser.add_argument('--chromSizes')
  parser.add_argument('--outDir')


def checkArguments(args, parser):
  """ Make sure all of the args are properly set for the default arguments.
  """
  # setting
  pairs = tuple((item, getattr(args, item)) for item in
                ['refGenome', 'genome',
                 'geneCheckBed', 'geneCheckBedDetails',
                 'alignment', 'sequence', 'chromSizes',
                 'outDir'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)
  # existence
  pairs = tuple((item, getattr(args, item)) for item in
                ['geneCheckBed', 'geneCheckBedDetails',
                 'alignment', 'sequence', 'chromSizes',
                 'outDir'])
  for name, value in pairs:
    if not os.path.exists(value):
      parser.error('--%s=%s does not exist' % (name, value))
  # regular file
  pairs = tuple((item, getattr(args, item)) for item in
                ['geneCheckBed', 'geneCheckBedDetails',
                 'alignment', 'sequence', 'chromSizes',
                 ])
  for name, value in pairs:
    if not os.path.isfile(value):
      parser.error('--%s=%s is not a file' % (name, value))
  # directory
  if not os.path.isdir(args.outDir):
    parser.error('--outDir=%s is not a directory' % args.outDir)


def getSequences(infile):
  """ Given a path to a fasta file, return a dictionary of Sequence objects
  keyed on the sequence name
  """
  seqDict = {}
  seq = None
  with open(infile, 'r') as f:
    for seq in readSequence(f):
      seqDict[seq.name] = seq
  return seqDict


def readSequence(infile):
  """ provide an iterator that reads through fasta files.
  """
  buff = None
  eat_buffer = True
  while True:
    if eat_buffer:
      # new record
      if buff is None:
        header = ''
        while not header.startswith('>'):
          header = infile.readline().strip()
      else:
        header = buff
      assert(header.startswith('>'))
      name = header.replace('>', '').strip()
      seq = ''
    line = infile.readline().strip()
    if line:
      if line.startswith('>'):
        # stop processing the record, store this line.
        buff = line
        eat_buffer = True
        yield Sequence(name, seq)
      else:
        eat_buffer = False
        seq += line
    else:
      # eof
      if buff is not None:
        buff = None
        yield Sequence(name, seq)
      else:
        if seq != '':
          yield Sequence(name, seq)
          name = ''
          seq = ''
        else:
          return

def getChromSizes(infile):
  """ read a chrom sizes file and return a dict keyed by names valued by ints.
  """
  chromDict = {}
  with open(infile, 'r') as f:
    for line in f:
      line = line.strip()
      if line == '':
        continue
      data = line.split()
      chromDict[data[0]] = int(data[1])
  return chromDict


def getAlignment(infile):
  """ read a PSL file and return a list of PSL
  """
  pass
