"""
convenience library for assisting filters.
"""
import os

class Sequence(object):
  """ Represents a sequence of DNA.
  """
  def __init__(self, name):
    self.name = name  # chromosome or scaffold name
    self._sequence = ''  # ACGTs
  def addSequence(self, seq):
    self._sequence += seq
  def setSequence(self, seq):
    self._sequence = seq
  def getSequence(self):
    return self._sequence
  def getLength(self):
    return len(self._sequence)


def InitializeArguments(parser):
  parser.add_argument('--refGenome')
  parser.add_argument('--genome')
  parser.add_argument('--geneCheckBed')
  parser.add_argument('--geneCheckBedDetails')
  parser.add_argument('--alignment')
  parser.add_argument('--sequence')
  parser.add_argument('--chromSizes')
  parser.add_argument('--outDir')


def CheckArguments(args, parser):
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


def getSequences(seqFile):
  """ Given a path to a fasta file, return a dictionary of Sequence objects
  keyed on the sequence name
  """
  seqDict = {}
  seq = None
  with open(seqFile, 'r') as f:
    for line in f:
      line = line.strip()
      if line == '':
        continue
      if line.startswith('>'):
        name = line.replace('>', '').strip()
        seq = Sequence(name)
        seqDict[name] = seq
        continue
      seq.addSequence(line)
  return seqDict
