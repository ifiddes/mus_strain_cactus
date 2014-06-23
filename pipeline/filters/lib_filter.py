"""
convenience library for assisting filters.
"""
from argparse import ArgumentTypeError
import os
import sys


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
  def setUpper(self):
    self._sequence = self._sequence.upper()
  def sliceSequence(self, start, stop):
    """ return the proper slice of the sequence.
    BED format coordinates: 0 based start, stop is exclusive
    [start, stop). E.g. Sequence.sliceSequence(0, 3) returns a string length 3.
    """
    return self._sequence[start:stop]


class PslRow(object):
  """ Represents a single row in a PSL file.
  http://genome.ucsc.edu/FAQ/FAQformat.html#format2
  """
  def __init__(self, line):
    data = line.split()
    assert(len(data) == 21)
    self.matches = int(data[0])
    self.misMatches = int(data[1])
    self.repMatches = int(data[2])
    self.nCount = int(data[3])
    self.qNumInsert = int(data[4])
    self.qBaseInsert = int(data[5])
    self.tNumInsert = int(data[6])
    self.tBaseInsert = int(data[7])
    self.strand = data[8]
    self.qName = data[9]
    self.qSize = int(data[10])
    self.qStart = int(data[11])
    self.qEnd = int(data[12])
    self.tName = data[13]
    self.tSize = int(data[14])
    self.tStart = int(data[15])
    self.tEnd = int(data[16])
    self.blockCount = int(data[17])
    # lists of ints
    self.blockSizes = map(int, [x for x in data[18].split(',') if x])
    self.qStarts = map(int, [x for x in data[19].split(',') if x])
    self.tStarts = map(int, [x for x in data[20].split(',') if x])
  def hashkey(self):
    """ return a string to use as dict key.
    """
    return '%s_%s_%d_%d' % (self.qName, self.tName, self.tStart, self.tEnd)


"""The following data types are used for iterating over gene-check-detail and
gene-check bed files.
An example of entries from such files:

[benedict@hgwdev tracks]$ more Rattus.coding.gene-check-details.bed
1       2812370 2812372 noStop/ENSMUST00000065527.4
1       2812370 2812372 noStop/ENSMUST00000095795.4

[benedict@hgwdev tracks]$ more Rattus.coding.gene-check.bed
1       2812346 3113743 ENSMUST00000178026.1    0       -       2812370 3038729
128,0,0 9       54,2,89,249,90,165,105,13,45    0,58,62,698,1209,1305,226292,301
050,301352
1       2812346 3113783 ENSMUST00000095795.4    0       -       2812370 3038729
128,0,0 9       54,2,89,249,197,52,105,13,85    0,58,62,698,1209,1418,226292,301
050,301352
"""

class ChromosomeInterval(object):
  """ Represents an interval of a chromosome. BED coordinates, strand is True,
  False or None (if no strand)
  """
  def __init__(self, chromosome, start, stop, strand):
    self.chromosome = str(chromosome)
    self.start = int(start)
    self.stop = int(stop)
    self.strand = strand

  def __eq__(self, other):
    return (self.chromosome == other.chromosome and
            self.start == other.start and
            self.stop == other.stop and
            self.strand == other.strand)

  def __cmp__(self, cI):
    return cmp((self.chromosome, self.start, self.stop, self.strand),
               (cI.chromosome, cI.start, cI.stop, cI.strand))


class TranscriptAnnotation(object):
  """ Represents an annotation of a transcript, from one of the
  classification bed files
  """
  def __init__(self, chromosomeInterval, name, label):
    self.chromosomeInterval = chromosomeInterval
    self.name = str(name)
    self.labels = list(label)  # list of strings

  def addLabel(self, label):
    """ maintain the ordering of self.labels but prevent duplicates
    """
    if label not in self.labels:
      self.labels.append(label)

  def bedString(self):
      return "\t".join([self.chromosomeInterval.chromosome,
                        str(self.chromosomeInterval.start),
                        str(self.chromosomeInterval.stop),
                        "/".join(self.labels + [self.name])])

  def __eq__(self, other):
    return (self.chromosomeInterval == other.chromosomeInterval and
            self.name == other.name and
            self.labels == other.labels)

  def __cmp__(self, annotation):
    """ Sort by chromosome interval, then name
    """
    return cmp((self.chromosomeInterval, self.name),
               (annotation.chromosomeInterval, annotation.name))

class Transcript(object):
  """ Represent a transcript and its annotations
  """
  def __init__(self, chromosomeInterval, name, exons, annotations,
               score, thickStart, thickEnd, itemRgb):
    self.chromosomeInterval = chromosomeInterval
    self.name = str(name)
    self.exons = exons # list of chromosome intervals
    self.annotations = annotations # list of TranscriptAnnotation() objects
    # Bed fields
    self.score = score
    self.thickStart = thickStart
    self.thickEnd = thickEnd
    self.itemRgb = itemRgb

  def __eq__(self, other):
    return (self.chromosomeInterval == other.chromosomeInterval and
            self.name == other.name and
            self.exons == other.exons and
            self.annotations == other.annotations and
            self.score == other.score and
            self.thickStart == other.thickStart and
            self.thickEnd == other.thickEnd and
            self.itemRgb == other.itemRgb)

  def hashkey(self):
    """ return a string to use as dict key.
    """
    return '%s_%s_%d_%d' % (self.name, self.chromosomeInterval.chromosome,
                            self.chromosomeInterval.start,
                            self.chromosomeInterval.stop)

  def bedString(self):
    """ Write a transcript object to the given file.
    """
    strandChar = "-"
    if self.chromosomeInterval.strand:
        strandChar = "+"
    return "\t".join(
      [self.chromosomeInterval.chromosome,
       str(self.chromosomeInterval.start),
       str(self.chromosomeInterval.stop),
       self.name, str(self.score), strandChar,
       str(self.thickStart), str(self.thickEnd),
       self.itemRgb, str(len(self.exons)),
       ",".join([str(exon.stop - exon.start) for exon in self.exons]),
       ",".join([str(exon.start - self.chromosomeInterval.start)
                 for exon in self.exons])])

  def __cmp__(self, transcript):
    return cmp((self.chromosomeInterval, self.name),
               (transcript.chromosomeInterval, transcript.name))


def DirType(d):
  """ given a string path to a directory, D, verify it can be used.
  """
  d = os.path.abspath(d)
  if not os.path.exists(d):
    raise ArgumentTypeError('DirType:%s does not exist' % d)
  if not os.path.isdir(d):
    raise ArgumentTypeError('DirType:%s is not a directory' % d)
  if os.access(d, os.R_OK):
    return d
  else:
    raise ArgumentTypeError('DirType:%s is not a readable dir' % d)


def FileType(f):
  """ given a string path to a file, F, verify it can be used.
  """
  f = os.path.abspath(f)
  if not os.path.exists(f):
    raise ArgumentTypeError('FileType:%s does not exist' % f)
  if not os.path.isfile(f):
    raise ArgumentTypeError('FileType:%s is not a regular file' % f)
  if os.access(f, os.R_OK):
    return f
  else:
    raise ArgumentTypeError('FileType:%s is not a readable file' % f)


def initializeArguments(parser):
  """ given an argparse ArgumentParser object, add in the default arguments.
  """
  parser.add_argument('--refGenome', type=str)
  parser.add_argument('--genome', type=str)
  parser.add_argument('--geneCheckBed', type=FileType)
  parser.add_argument('--geneCheckBedDetails', type=FileType)
  parser.add_argument('--alignment', type=FileType)
  parser.add_argument('--sequence', type=FileType)
  parser.add_argument('--chromSizes', type=FileType)
  parser.add_argument('--outDir', type=DirType)


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
  # record the issuing command
  with open(os.path.join(args.outDir, 'command.log'), 'w') as f:
    f.write('%s %s %s %s %s %s %s %s %s\n'
            % (sys.argv[0],
               args.refGenome,
               args.genome,
               args.geneCheckBed,
               args.geneCheckBedDetails,
               args.alignment,
               args.sequence,
               args.chromSizes,
               args.outDir,
               ))


def reverseComplement(seq):
  """ Given a sequence, return the reverse complement.
  """
  pairs = [('a', 't'), ('g', 'c')]
  complement = {}
  for a, b in pairs:
    complement[a], complement[b] = b, a
    complement[a.upper()], complement[b.upper()] = b.upper(), a.upper()
  complement['-'] = '-'
  seq = seq[::-1]  # reverse
  seq = ''.join(map(lambda s: complement[s], seq))
  return seq


def getSequences(infile, upper=False):
  """ Given a path to a fasta file, return a dictionary of Sequence objects
  keyed on the sequence name
  """
  seqDict = {}
  seq = None
  with open(infile, 'r') as f:
    for seq in readSequence(f):
      seqDict[seq.name] = seq
      if upper:
        seqDict[seq.name].setUpper()
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


def getAlignments(infile):
  """ read a PSL file and return a list of PslRow objects
  """
  psls = []
  with open(infile, 'r') as f:
    for psl in readPsls(f):
      psls.append(psl)
  return psls


def readPsls(infile):
  """ provide an iterator that reads through psl files.
  """
  while True:
    line = infile.readline().strip()
    if line == '':
      return
    yield PslRow(line)


def getTranscripts(bedFile, bedDetailsFile):
  """ Given a path to a standard BED file and a details BED, return a list of
  Transcript objects.
  """
  transcripts = []
  bf = open(bedFile, 'r')
  bdf = open(bedDetailsFile, 'r')
  for t in transcriptIterator(bf, bdf):
    transcripts.append(t)
  return transcripts


def tokenizeBedStream(bedStream):
  """ Iterator through bed file, returning lines as list of tokens
  """
  for line in bedStream:
    if line != '':
      tokens = line.split()
      yield tokens

def normalizeAnnotation(transcriptAnnotation):
  """ Normalizes the transcript annotation labels.
  This is meant to munge the labels of transcript annotations according to
  hacky needs of the input annotation type labels.
  """
  if (len(transcriptAnnotation.labels) > 1 and
      (transcriptAnnotation.labels[0].endswith('Splice') and
       len(transcriptAnnotation.labels[1].split('.')) == 3)):
    # try to find lists like ['unknownUtrSplice', 'CC..AC']
    newLabels = [ "_".join(transcriptAnnotation.labels[:2])]
    newLabels += transcriptAnnotation.labels[2:]
    transcriptAnnotation.labels = newLabels
  elif (len(transcriptAnnotation.labels) > 1 and
        transcriptAnnotation.labels[0] == 'orfStop' and
        transcriptAnnotation.labels[1] in ['TAA', 'TAG', 'TGA']):
    # try to find lists like ['orfStop', 'TAG']
    newLabels = [ "_".join(transcriptAnnotation.labels[:2])]
    newLabels += transcriptAnnotation.labels[2:]
    transcriptAnnotation.labels = newLabels

def transcriptIterator(transcriptsBedStream, transcriptDetailsBedStream):
  """ Iterates over the transcripts detailed in the two streams, producing
  Transcript objects. Streams are any iterator that returns bedlines or empty
  strings.
  """
  transcriptsAnnotations = {}
  for tokens in tokenizeBedStream(transcriptDetailsBedStream):
    assert len(tokens) == 4
    tA = TranscriptAnnotation(
      ChromosomeInterval(tokens[0], tokens[1], tokens[2], None),
      tokens[3].split('/')[-1], tokens[3].split('/')[:-1])
    normalizeAnnotation(tA)
    key = (tA.name, tA.chromosomeInterval.chromosome)
    if key not in transcriptsAnnotations:
      transcriptsAnnotations[key] = []
    transcriptsAnnotations[key].append(tA)

  for tokens in tokenizeBedStream(transcriptsBedStream):
    assert len(tokens) == 12
    # Transcript
    name = tokens[3]
    # Get the chromosome interval
    assert tokens[5] in ['+', '-']
    cI = ChromosomeInterval(tokens[0], tokens[1], tokens[2], tokens[5] == '+')
    # Get the exons
    def getExons(exonNumber, blockSizes, blockStarts):
      assert exonNumber == len(blockSizes)
      assert exonNumber == len(blockStarts)
      return [ChromosomeInterval(
          cI.chromosome, cI.start + int(blockStarts[i]),
          cI.start + int(blockStarts[i]) + int(blockSizes[i]), cI.strand)
              for i in range(exonNumber)]
    exons = getExons(int(tokens[9]),
                     tokens[10].split(','), tokens[11].split(','))
    # Get the name annotations
    annotations = []
    key = (name, cI.chromosome)
    if key in transcriptsAnnotations:
      annotations = transcriptsAnnotations[key]
    yield Transcript(
      cI, name, exons, annotations,
      int(tokens[4]), int(tokens[6]),
      int(tokens[7]), tokens[8])


def writeAllBeds(transcripts, args):
  """ Convenience function to take a list of TRANSCRIPTS and write out the
  standard and details beds to the expected location in args.outDir.
  """
  out_bed, out_bed_details = getBedOutFiles(args)
  writeTranscriptBedFile(transcripts, out_bed)
  writeDetailsBedFile(transcripts, out_bed_details)


def getBedOutFiles(args):
  """ Return the paths of the new bed and bedDetail files.
  """
  bed = os.path.join(args.outDir, 'out.bed')
  bedDetails = os.path.join(args.outDir, 'out_details.bed')
  return bed, bedDetails


def writeDetailsBedFile(transcripts, detailsBedFile):
  """Writes out a details bed file for a set of transcripts - that is the set of
  annotations of the transcripts. The bed file must be in chromosome order.
  """
  annotations = reduce(lambda x, y : x + y,
                       [transcript.annotations for transcript in transcripts])
  annotations.sort()
  annotationsFileHandle = open(detailsBedFile, 'w')
  for annotation in annotations:
    annotationsFileHandle.write(annotation.bedString() + "\n")
  annotationsFileHandle.close()


def writeTranscriptBedFile(transcripts, bedFile):
  """Writes out an bed file for a set of transcripts.
  """
  transcripts = transcripts[:]
  transcripts.sort()
  bedFileHandle = open(bedFile, 'w')
  for transcript in transcripts:
    bedFileHandle.write(transcript.bedString() + "\n")
  bedFileHandle.close()

