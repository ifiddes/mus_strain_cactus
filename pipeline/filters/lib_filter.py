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
    assert(len(data) == 21)
    self.matches = int(data[0])
    self.misMatches = int(data[1])
    self.repMatches = int(data[2])
    self.nCount = int(data[3])
    self.qNumInsert = int(data[4])
    self.qBaseInsert = int(data[5])
    self.tNumInsert = int(data[6])
    self.tBaseInsert = int(data[7])
    self.strand = data[7]
    self.qName = data[8]
    self.qSize = int(data[10])
    self.qStart = int(data[11])
    self.qEnd = int(data[12])
    self.tName = data[12]
    self.tSize = int(data[14])
    self.tStart = int(data[15])
    self.tEnd = int(data[16])
    self.blockCount = int(data[17])
    self.blockSizes = data[18]
    self.qStarts = data[19]
    self.tStarts = data[20]


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
  False or NULL (if no strand)
  """
  def __init__(self, chromosome, start, stop, strand):
    self.chromosome = str(chromosome)
    self.start = int(start)
    self.stop = int(stop)
    self.strand = strand

  def __cmp__(self, cI):
    return cmp((self.chromosome, self.start, self.stop, self.strand), (cI.chromosome, cI.start, cI.stop, cI.strand))


class TranscriptAnnotation(object):
  """ Represents an annotation of a transcript, from one of the
  classification bed files
  """
  def __init__(self, chromosomeInterval, name, annotation):
    self.chromosomeInterval = chromosomeInterval
    self.name = str(name)
    self.annotation = list(annotation)

  def bedString(self):
      return "\t".join([self.chromosomeInterval.chromosome,
                        str(self.chromosomeInterval.start),
                        str(self.chromosomeInterval.stop),"/".join(self.annotation + [self.name])])

  def __cmp__(self, annotation):
    """Sort by chromosome interval, then name
    """
    return cmp((self.chromosomeInterval, self.name), (annotation.chromosomeInterval, annotation.name))

class Transcript(object):
  """ Represent a transcript and its annotations
  """
  def __init__(self, chromosomeInterval, name, exons, annotations,
               score, thickStart, thickEnd, itemRgb):
    self.chromosomeInterval = chromosomeInterval
    self.name = str(name)
    self.exons = exons #Is a list of chromosome intervals
    self.annotations = annotations #Is a list of transcript annotations
    #Bed fields
    self.score = score
    self.thickStart = thickStart
    self.thickEnd = thickEnd
    self.itemRgb = itemRgb

  def bedString(self):
    """Write a transcript object to the given file.
    """
    strandChar = "-"
    if self.chromosomeInterval.strand:
        strandChar = "+"
    return "\t".join([self.chromosomeInterval.chromosome,
                      str(self.chromosomeInterval.start),
                      str(self.chromosomeInterval.stop),
                      self.name, str(self.score), strandChar,
                      str(self.thickStart), str(self.thickEnd),
                      self.itemRgb, str(len(self.exons)),
                      ",".join([ str(exon.stop - exon.start) for exon in self.exons]),
                      ",".join([ str(exon.start - self.chromosomeInterval.start) for exon in self.exons])])

  def __cmp__(self, transcript):
    return cmp((self.chromosomeInterval, self.name), (transcript.chromosomeInterval, transcript.name))

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


def tokenizeBedStream(bedStream):
  """ Iterator through bed file, returning lines as list of tokens
  """
  for line in bedStream:
    if line != '':
      tokens = line.split()
      yield tokens


def transcriptIterator(transcriptsBedStream, transcriptDetailsBedStream):
  """ Iterates over the transcripts detailed in the two streams, producing
  Transcript objects. Streams are any iterator that returns bedlines or empty
  strings.
  """
  transcriptsAnnotations = {}
  for tokens in tokenizeBedStream(transcriptDetailsBedStream):
    assert len(tokens) == 4
    tA = TranscriptAnnotation(ChromosomeInterval(
        tokens[0], tokens[1], tokens[2], None), tokens[3].split("/")[-1], tokens[3].split("/")[:-1])
    if tA.name not in transcriptsAnnotations:
      transcriptsAnnotations[tA.name] = []
    transcriptsAnnotations[tA.name].append(tA)

  for tokens in tokenizeBedStream(transcriptsBedStream):
    assert len(tokens) == 12
    # Transcript
    name = tokens[3]
    # Get the chromosome interval
    assert tokens[5] in ('+', '-')
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
                     tokens[10].split(","), tokens[11].split(","))
    # Get the name annotations
    annotations = []
    if name in transcriptsAnnotations:
      annotations = transcriptsAnnotations[name]
    yield Transcript(cI, name, exons, annotations, int(tokens[4]), int(tokens[6]), int(tokens[7]), tokens[8])
    
def writeDetailsBedFile(transcripts, detailsBedFile):
  """Writes out a details bed file for a set of transcripts - that is the set of
  annotations of the transcripts. The bed file must be in chromosome order."""
  annotations = reduce(lambda x, y : x + y, [ transcript.annotations for transcript in transcripts ])
  annotations.sort()
  annotationsFileHandle = open(detailsBedFile, "w")
  for annotation in annotations:
    annotationsFileHandle.write(annotation.bedString())
  annotationsFileHandle.close()

def writeTranscriptBedFile(transcripts, bedFile):
  """Writes out an bed file for a set of transcripts.
  """
  transcripts = transcripts[:]
  transcripts.sort()
  bedFileHandle = open(bedFile, "w")
  for transcript in transcripts:
    bedFileHandle.write(transcript.bedString())
  bedFileHandle.close()

