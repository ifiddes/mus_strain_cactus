""" Test the lib_filter classes and functions
"""
from glob import glob
import os
import shutil
import string
import subprocess
import sys
import unittest
import lib_filter


def makeTempDirParent():
  """ make the parent temp dir directory
  """
  if not os.path.exists(os.path.join(os.curdir, '.tempTestDir')):
    os.mkdir(os.path.join(os.curdir, '.tempTestDir'))


def removeTempDirParent():
  """ remove the parent temp dir directory
  """
  if os.path.exists(os.path.join(os.curdir, '.tempTestDir')):
    shutil.rmtree(os.path.join(os.curdir, '.tempTestDir'))


def makeTempDir(name=None):
  """ make the directory where all temporary test files will be stored.
  """
  makeTempDirParent()
  charSet = string.ascii_lowercase + '123456789'
  if name is None:
    while True:
      name = '%s_%s' % (''.join(random.choice(charSet) for x in xrange(4)),
                        ''.join(random.choice(charSet) for x in xrange(4)))
      if not os.path.exists(os.path.join(os.curdir, '.tempTestDir', name)):
        break
  if not os.path.exists(os.path.join(os.curdir, '.tempTestDir', name)):
    os.mkdir(os.path.join(os.curdir, '.tempTestDir', name))
  return os.path.join(os.curdir, '.tempTestDir', name)


def removeDir(dirpath):
  """ destroy a directory
  """
  if os.path.exists(dirpath):
    shutil.rmtree(dirpath)
  if glob(os.path.join(os.path.dirname(dirpath), '*')) == []:
    # if this is the last tempDir to be destroyed, destroy the parent.
    removeTempDirParent()


def which(program):
  """which() acts like the unix utility which, but is portable between os.
  If the program does not exist in the PATH then 'None' is returned.
  """
  def is_exe(fpath):
    return os.path.exists(fpath) and os.access(fpath, os.X_OK)

  fpath, fname = os.path.split(program)
  if fpath != '':
    if is_exe(program):
      return program
  else:
    for path in os.environ['PATH'].split(os.pathsep):
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return exe_file
  return None


def runCommands(cmds, localTempDir, inPipes=None, outPipes=None, errPipes=None):
  """ Run commands from CMDS list.
  """
  if inPipes is None:
    inPipes = [None] * len(cmds)
  if outPipes is None:
    outPipes = [None] * len(cmds)
  if errPipes is None:
    errPipes = [None] * len(cmds)
  for i, c in enumerate(cmds, 0):
    if inPipes[i] is None:
      sin = None
    else:
      sin = subprocess.PIPE
    if outPipes[i] is None:
      sout = None
    else:
      sout = subprocess.PIPE
    if errPipes[i] is None:
      serr = None
    else:
      serr = subprocess.PIPE
    p = subprocess.Popen(c, cwd=localTempDir, stdin=sin,
                         stdout=sout, stderr=serr)
    if inPipes[i] is None:
      sin = None
    else:
      if not os.path.exists(inPipes[i]):
        raise IOError('Unable to locate inPipe file: %s for command %s'
                      % (inPipes[i], ' '.join(c)))
      sin = open(inPipes[i], 'r').read()
    if outPipes[i] is None:
      pout, perr = p.communicate(sin)
      handleReturnCode(p.returncode, cmds[i])
    else:
      with open(outPipes[i], 'w') as f:
        f.write(p.communicate(sin)[0])
      handleReturnCode(p.returncode, cmds[i])


def handleReturnCode(retcode, cmd):
  """ handle the return codes from runCommands
  """
  if not isinstance(retcode, int):
    raise TypeError('handleReturnCode takes an integer for '
                    'retcode, not a %s.' % retcode.__class__)
  if not isinstance(cmd, list):
    raise TypeError('handleReturnCode takes a list for '
                    'cmd, not a %s.' % cmd.__class__)
  if retcode:
    if retcode < 0:
      raise RuntimeError('Experienced an error while trying to execute: '
                         '%s SIGNAL:%d' %(' '.join(cmd), -retcode))
    else:
      raise RuntimeError('Experienced an error while trying to execute: '
                         '%s retcode:%d' %(' '.join(cmd), retcode))


def createSequenceFile(sequences, tmpDir):
  """ given a dict (key is name, value is sequence) return path to temp file.
  """
  seqfile = os.path.join(tmpDir, 'seq.fa')
  with open(seqfile, 'w') as f:
    for name in sequences:
      f.write('>%s\n%s' % (name, sequences[name]))
  return seqfile


def createAlignmentFile(alignments, tmpDir):
  """ given a list of alignments, return path to a temp file.
  """
  alnfile = os.path.join(tmpDir, 'aln.psl')
  with open(alnfile, 'w') as f:
    for a in alignments:
      f.write('%s\n' % a)
  return alnfile


def bedLine(chrom, chromStart, chromEnd, name, score=None, strand=None,
            thickStart=None, thickEnd=None, itemRgb=None, blockCount=None,
            blockSizes=None, blockStarts=None):
  """ Give the fields, create a bed line string
  """

  s = ('%s %d %d %s'
       % (chrom, chromStart, chromEnd, name))
  if score is not None:
    for v in [strand, thickStart, thickEnd, itemRgb,
              blockCount, blockSizes, blockStarts]:
      assert(v is not None)
    s += (' %d %s %d %d %s %d %s %s'
          % (score, strand, thickStart, thickEnd, itemRgb, blockCount,
          blockSizes, blockStarts))
  return s


class sequenceGetterTests(unittest.TestCase):
  def test_getSequences(self):
    """ getSequences must read a fasta and return a dict of Sequence objects.
    """
    sequences = {'chrA':
                   ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n'
                    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n'
                    'ACGT')
                  ,
                 'bannana.chr0':
                   ('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n'
                    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACTTT\n'
                    'ACGTACGTACGTACGTACGTACGTACGTACG\n')
                  }
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('getSequences'))
    testFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(testFile)
    for name in sequences:
      seq = sequences[name]
      seq = seq.replace('\n', '').strip()
      self.assertTrue(name in seqDict)
      self.assertEqual(seqDict[name].getLength(), len(seq))
      self.assertEqual(seqDict[name].getSequence(), seq)
    for name in seqDict:
      self.assertTrue(name in sequences)
    self.addCleanup(removeDir, tmpDir)


class alignmentGetterTests(unittest.TestCase):
  def test_getAlignment(self):
    """ getAlignment must read a psl file and return a list of PslRow objects.
    """
    alignments = [
      '141 0 0 0 0 0 0 0 - ENSMUST00000178550.1 141 0 141 scaffold-11326 33702 21065 21206 1 141, 0, 21065,',
      '309 0 0 0 0 0 0 0 - ENSMUST00000179623.1 309 0 309 scaffold-1475 11716 9284 9593 1 309, 0, 9284,',
      '700 0 0 0 3 5 2 4 - ENSMUST00000179112.1 705 0 705 scaffold-189833 540197 335509 336213 4 14,129,140,417, 0,15,146,288, 335509,335523,335654,335796,',
                  ]
    fields = ['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert',
              'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName',
              'qSize', 'qStart', 'qEnd', 'tName', 'tSize',
              'tStart', 'tEnd', 'blockCount']
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('getAlignments'))
    testFile = createAlignmentFile(alignments, tmpDir)
    libAlignments = lib_filter.getAlignments(testFile)
    for i, a in enumerate(alignments, 0):
      data = a.split()
      for j, field in enumerate(fields, 0):
        self.assertEqual(data[j], str(getattr(libAlignments[i], field)))
      for j, field in enumerate(['blockSizes', 'qStarts', 'tStarts'], 18):
        self.assertEqual(
          data[j], ','.join(map(str, getattr(libAlignments[i], field))) + ',')
    self.addCleanup(removeDir, tmpDir)


class transcriptIteratorTests(unittest.TestCase):
  def test_transcriptIterator(self):
    """ tests the transcriptIterator function
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113743, 'ENSMUST00000178026.1', 0, '-', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,90,165,105,13,45',
        '0,58,62,698,1209,1305,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113783, 'ENSMUST00000095795.4', 0, '+', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,197,52,105,13,85',
        '0,58,62,698,1209,1418,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        'scaffold-100021', 466, 4248, 'ENSMUST00000034053.5', 0, '-', 466, 4248,
        '128,0,0', 2, '85,152', '0,3630'))
    transcriptBedLines.append(bedLine(
        'scaffold-138877', 4903, 5091, 'ENSMUST00000034053.5', 0, '-', 4903,
        4996, '128,0,0', 1, '188', '0'))
    transcriptBedLines.append(bedLine(
        'scaffold-2051', 13759, 24866, 'ENSMUST00000034053.5', 0, '-', 14291,
        24866, '128,0,0', 4, '722,112,131,188', '0,1977,4316,10919'))

    transcriptDetailsBedLines = []
    transcriptDetailsBedLines.append(bedLine(
        '1', 2812346, 2812349, 'noStop/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        '1', 3113780, 3113783, 'noStart/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 466, 469, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 4245, 4248, 'noStart/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-138877', 4903, 4906, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-2051', 24863, 24866, 'noStart/ENSMUST00000034053.5'))
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    names = ['ENSMUST00000178026.1', 'ENSMUST00000095795.4',
             'ENSMUST00000034053.5', 'ENSMUST00000034053.5',
             'ENSMUST00000034053.5']
    self.assertEqual(len(transcripts), 5)
    for attr, values in [('name', names),
                         ('thickStart', [2812370, 2812370, 466, 4903, 14291]),
                         ('thickEnd', [3038729, 3038729, 4248, 4996, 24866]),
                         ('itemRgb', ['128,0,0', '128,0,0',
                                      '128,0,0', '128,0,0','128,0,0',]),
                         ('chromosomeInterval',
                          [
          lib_filter.ChromosomeInterval('1', 2812346, 3113743, False),
          lib_filter.ChromosomeInterval('1', 2812346, 3113783, True),
          lib_filter.ChromosomeInterval('scaffold-100021', 466, 4248, False),
          lib_filter.ChromosomeInterval('scaffold-138877', 4903, 5091, False),
          lib_filter.ChromosomeInterval('scaffold-2051', 13759, 24866, False),
          ]
                          ),
                         ]:
      for i, n in enumerate(values, 0):
        self.assertEquals(getattr(transcripts[i], attr), n)
    # exons
    self.assertEqual(len(transcripts[0].exons), 9)
    self.assertEqual(transcripts[0].exons[0].chromosome, '1')
    self.assertEqual(transcripts[0].exons[0].start, 2812346)
    self.assertEqual(transcripts[0].exons[0].stop, 2812346 + 54)
    self.assertEqual(transcripts[0].exons[0].strand, False)
    self.assertEqual(transcripts[1].exons[0].chromosome, '1')
    self.assertEqual(transcripts[1].exons[0].start, 2812346)
    self.assertEqual(transcripts[1].exons[0].stop, 2812346 + 54)
    self.assertEqual(transcripts[1].exons[0].strand, True)
    # annotations
    testAnnots = [
      [],
      [lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('1', 2812346, 2812349, None),
        'ENSMUST00000095795.4', ['noStop']),
       lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('1', 3113780, 3113783, None),
        'ENSMUST00000095795.4', ['noStart']),
       ],
      [lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('scaffold-100021', 466, 469, None),
        'ENSMUST00000034053.5', ['noStop']),
       lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('scaffold-100021', 4245, 4248, None),
        'ENSMUST00000034053.5', ['noStart']),
       ],
      [lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('scaffold-138877', 4903, 4906, None),
        'ENSMUST00000034053.5', ['noStop']),
       ],
      [lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('scaffold-2051', 24863, 24866, None),
        'ENSMUST00000034053.5', ['noStart']),
       ],
      ]
    self.assertEqual(len(transcripts), len(testAnnots))
    for i in xrange(0, len(transcripts)):
      self.assertEqual(transcripts[i].annotations, testAnnots[i])
    # Check print functions
    self.assertEquals(transcripts[0].bedString().split(), transcriptBedLines[0].split())
    self.assertEquals(transcripts[1].bedString().split(), transcriptBedLines[1].split())
    self.assertEquals(transcripts[1].annotations[0].bedString().split(), transcriptDetailsBedLines[0].split())
    self.assertEquals(transcripts[1].annotations[1].bedString().split(), transcriptDetailsBedLines[1].split())
    # Check sort function for transcripts
    transcripts.reverse()
    names.reverse()
    for i in xrange(0, len(transcripts)):
      self.assertEquals(transcripts[i].name, names[i])
    transcripts.sort()
    names.reverse()  # ASSUME THAT NAMES AND BED LINES IN TEST ARE ALREADY IN SORTED ORDER
    for i in xrange(0, len(transcripts)):
      self.assertEquals(transcripts[i].name, names[i])

  def test_transcriptWriter(self):
    """ transcripts should be written out correctly.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113743, 'ENSMUST00000178026.1', 0, '-', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,90,165,105,13,45',
        '0,58,62,698,1209,1305,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113783, 'ENSMUST00000095795.4', 0, '+', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,197,52,105,13,85',
        '0,58,62,698,1209,1418,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        'scaffold-100021', 466, 4248, 'ENSMUST00000034053.5', 0, '-', 466, 4248,
        '128,0,0', 2, '85,152', '0,3630'))
    transcriptBedLines.append(bedLine(
        'scaffold-138877', 4903, 5091, 'ENSMUST00000034053.5', 0, '-', 4903,
        4996, '128,0,0', 1, '188', '0'))
    transcriptBedLines.append(bedLine(
        'scaffold-2051', 13759, 24866, 'ENSMUST00000034053.5', 0, '-', 14291,
        24866, '128,0,0', 4, '722,112,131,188', '0,1977,4316,10919'))
    transcriptDetailsBedLines = []
    transcriptDetailsBedLines.append(bedLine(
        '1', 2812346, 2812349, 'noStop/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        '1', 3113780, 3113783, 'noStart/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 466, 469, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 4245, 4248, 'noStart/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-138877', 4903, 4906, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-2051', 24863, 24866, 'noStart/ENSMUST00000034053.5'))
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('transcriptWriter'))
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(outBed, outDetailsBed)
    # test equality.
    self.assertEquals(len(transcripts), len(writtenTranscripts))
    transcripts.sort()
    writtenTranscripts.sort()
    for i in xrange(0, len(transcripts)):
      self.assertEquals(transcripts[i].chromosomeInterval,
                        writtenTranscripts[i].chromosomeInterval)
      self.assertEquals(transcripts[i].annotations,
                        writtenTranscripts[i].annotations)
      self.assertEquals(transcripts[i], writtenTranscripts[i])
    # cleanup
    self.addCleanup(removeDir, tmpDir)

  def test_transcriptObjectManagement(self):
    """ Transcript and TranscriptAnnotation counts should not change.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113743, 'ENSMUST00000178026.1', 0, '-', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,90,165,105,13,45',
        '0,58,62,698,1209,1305,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113783, 'ENSMUST00000095795.4', 0, '+', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,197,52,105,13,85',
        '0,58,62,698,1209,1418,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        'scaffold-100021', 466, 4248, 'ENSMUST00000034053.5', 0, '-', 466, 4248,
        '128,0,0', 2, '85,152', '0,3630'))
    transcriptBedLines.append(bedLine(
        'scaffold-138877', 4903, 5091, 'ENSMUST00000034053.5', 0, '-', 4903,
        4996, '128,0,0', 1, '188', '0'))
    transcriptBedLines.append(bedLine(
        'scaffold-2051', 13759, 24866, 'ENSMUST00000034053.5', 0, '-', 14291,
        24866, '128,0,0', 4, '722,112,131,188', '0,1977,4316,10919'))
    # the following bed lines can cause problems because they represent
    # multilpe instances of a transcript on a single chromosome.
    transcriptBedLines += [
      'X 73726551 74227830 ENSMUST00000101560.3 0 + 73726551 74227830 128,0,0 12 110,22,6,14,4,15,21,18,143,6,21,23 0,500953,500979,500987,501006,501011,501028,501050,501069,501213,501220,501256',
      'X 73726990 74230505 ENSMUST00000101560.3 0 + 73726990 73727560 128,0,0 30 15,10,109,50,54,225,4,44,25,19,1,18,2,56,4,13,20,18,23,199,320,369,88,137,117,67,107,102,348,5 0,19,34,144,195,250,480,487,532,559,580,582,601,501428,501486,501491,501507,501529,501550,501576,501789,502110,502481,502571,502709,502827,502924,503055,503158,503510',
      'X 74211514 74218259 ENSMUST00000101560.3 0 + 74213333 74218259 128,0,0 8 40,15,18,8,148,254,70,133 0,1726,1742,1761,1770,1919,6540,6612',
      ]
    transcriptDetailsBedLines = []
    transcriptDetailsBedLines.append(bedLine(
        '1', 2812346, 2812349, 'noStop/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        '1', 3113780, 3113783, 'noStart/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 466, 469, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 4245, 4248, 'noStart/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-138877', 4903, 4906, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-2051', 24863, 24866, 'noStart/ENSMUST00000034053.5'))
    transcriptDetailsBedLines += [
      'X 73726551 73726554 noStart/ENSMUST00000101560.3',
      'X 73726551 73726661 frameMismatch/ENSMUST00000101560.3',
      'X 73726551 74227830 badFrame/ENSMUST00000101560.3',
      'X 73726990 73726993 noStart/ENSMUST00000101560.3',
      'X 73726990 73727560 badFrame/ENSMUST00000101560.3',
      'X 73727005 73727009 cdsGap/ENSMUST00000101560.3',
      'X 73727009 73727019 frameDiscontig/ENSMUST00000101560.3',
      'X 73727019 73727024 cdsGap/ENSMUST00000101560.3',
      'X 73727133 73727134 cdsGap/ENSMUST00000101560.3',
      'X 73727134 73727184 frameDiscontig/ENSMUST00000101560.3',
      'X 73727184 73727185 cdsGap/ENSMUST00000101560.3',
      'X 73727185 73727239 frameDiscontig/ENSMUST00000101560.3',
      'X 73727239 73727240 cdsGap/ENSMUST00000101560.3',
      'X 73727465 73727470 cdsGap/ENSMUST00000101560.3',
      'X 73727470 73727474 frameDiscontig/ENSMUST00000101560.3',
      'X 73727474 73727477 cdsGap/ENSMUST00000101560.3',
      'X 73727521 73727522 cdsGap/ENSMUST00000101560.3',
      'X 73727522 73727547 frameDiscontig/ENSMUST00000101560.3',
      'X 73727547 73727549 cdsGap/ENSMUST00000101560.3',
      'X 73727559 73727560 noStop/ENSMUST00000101560.3',
      'X 73727568 73727570 utrGap/ENSMUST00000101560.3',
      'X 73727571 73727572 utrGap/ENSMUST00000101560.3',
      'X 73727590 73727591 utrGap/ENSMUST00000101560.3',
      'X 73727593 74228418 unknownUtrSplice/GA..AA/ENSMUST00000101560.3',
      'X 74213255 74213256 utrGap/ENSMUST00000101560.3',
      'X 74213274 74213275 utrGap/ENSMUST00000101560.3',
      'X 74213283 74213284 utrGap/ENSMUST00000101560.3',
      'X 74213333 74218259 badFrame/ENSMUST00000101560.3',
      'X 74213432 74213433 cdsGap/ENSMUST00000101560.3',
      'X 74213433 74213687 frameDiscontig/ENSMUST00000101560.3',
      'X 74218124 74218126 cdsGap/ENSMUST00000101560.3',
      'X 74218258 74218259 noStop/ENSMUST00000101560.3',
      'X 74227504 74227526 frameDiscontig/ENSMUST00000101560.3',
      'X 74227526 74227530 cdsGap/ENSMUST00000101560.3',
      'X 74227530 74227536 frameDiscontig/ENSMUST00000101560.3',
      'X 74227536 74227538 cdsGap/ENSMUST00000101560.3',
      'X 74227552 74227557 cdsGap/ENSMUST00000101560.3',
      'X 74227561 74227562 cdsGap/ENSMUST00000101560.3',
      'X 74227562 74227577 frameDiscontig/ENSMUST00000101560.3',
      'X 74227577 74227579 cdsGap/ENSMUST00000101560.3',
      'X 74227579 74227600 frameDiscontig/ENSMUST00000101560.3',
      'X 74227600 74227601 cdsGap/ENSMUST00000101560.3',
      'X 74227601 74227619 frameDiscontig/ENSMUST00000101560.3',
      'X 74227619 74227620 cdsGap/ENSMUST00000101560.3',
      'X 74227620 74227763 frameDiscontig/ENSMUST00000101560.3',
      'X 74227763 74227764 cdsGap/ENSMUST00000101560.3',
      'X 74227764 74227770 frameDiscontig/ENSMUST00000101560.3',
      'X 74227770 74227771 cdsGap/ENSMUST00000101560.3',
      'X 74227792 74227807 cdsGap/ENSMUST00000101560.3',
      'X 74227829 74227830 noStop/ENSMUST00000101560.3',
      'X 74228474 74228476 utrGap/ENSMUST00000101560.3',
      'X 74228480 74228481 utrGap/ENSMUST00000101560.3',
      'X 74228494 74228497 utrGap/ENSMUST00000101560.3',
      'X 74228517 74228519 utrGap/ENSMUST00000101560.3',
      'X 74228537 74228540 utrGap/ENSMUST00000101560.3',
      'X 74228563 74228566 utrGap/ENSMUST00000101560.3',
      'X 74228765 74228779 utrGap/ENSMUST00000101560.3',
      'X 74229099 74229100 utrGap/ENSMUST00000101560.3',
      'X 74229469 74229471 utrGap/ENSMUST00000101560.3',
      'X 74229559 74229561 utrGap/ENSMUST00000101560.3',
      'X 74229698 74229699 utrGap/ENSMUST00000101560.3',
      'X 74229816 74229817 utrGap/ENSMUST00000101560.3',
      'X 74229884 74229914 unknownUtrSplice/AA..GT/ENSMUST00000101560.3',
      'X 74230021 74230045 unknownUtrSplice/AA..GA/ENSMUST00000101560.3',
      'X 74230147 74230148 utrGap/ENSMUST00000101560.3',
      'X 74230496 74230500 utrGap/ENSMUST00000101560.3',
      ]
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    originalTranscriptCount = len(transcripts)
    originalTranscriptAnnotationCount = sum(
      map(lambda t:len(t.annotations), transcripts))
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('transcriptWriter'))
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(outBed, outDetailsBed)
    newTranscriptCount = len(writtenTranscripts)
    newTranscriptAnnotationCount = sum(
      map(lambda t:len(t.annotations), writtenTranscripts))
    # test equality.
    self.assertEquals(originalTranscriptCount, newTranscriptCount)
    self.assertEquals(originalTranscriptAnnotationCount,
                      newTranscriptAnnotationCount)
    # cleanup
    self.addCleanup(removeDir, tmpDir)

  def test_transcriptReading_0(self):
    """ Transcript objects should be read in correctly
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'scaffold-444', 41415, 87033, 'ENSMUST00000169901.2', 0, '-', 41415,
        45156, '128,0,0', 5, '128,12,219,90,27', '0,131,3590,42232,45591'))
    transcriptBedLines.append(bedLine(
        'scaffold-444', 72633, 82553, 'ENSMUST00000169901.2', 0, '-', 72782,
        82485, '0,128,0', 5, '51,156,104,140,219', '0,129,4370,7482,9701'))
    transcriptDetailsBedLines = []
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-444', 41415, 41418,
        'noStop/alignmentPartialMap/hasOkCopies/count_1/ENSMUST00000169901.2'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-444', 41543, 41546,
        'cdsGap/hasOkCopies/count_1/ENSMUST00000169901.2'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-444', 72633, 82553,
        'hasBadCopies/count_1/ENSMUST00000169901.2'))
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    originalTranscriptCount = len(transcripts)
    originalTranscriptAnnotationCount = sum(
      map(lambda t:len(t.annotations), transcripts))
    # test transcript count
    self.assertEquals(originalTranscriptCount, 2)
    self.assertEquals(len(transcripts[0].annotations), 2)
    self.assertEquals(len(transcripts[1].annotations), 0)
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('transcriptWriter'))
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(outBed, outDetailsBed)
    newTranscriptCount = len(writtenTranscripts)
    newTranscriptAnnotationCount = sum(
      map(lambda t:len(t.annotations), writtenTranscripts))
    # test transcript count
    self.assertEquals(newTranscriptCount, 2)
    self.assertEquals(len(writtenTranscripts[0].annotations), 2)
    self.assertEquals(len(writtenTranscripts[1].annotations), 0)
    # test equality.
    self.assertEquals(originalTranscriptCount, newTranscriptCount)
    self.assertEquals(originalTranscriptAnnotationCount,
                      newTranscriptAnnotationCount)
    # cleanup
    self.addCleanup(removeDir, tmpDir)


if __name__ == '__main__':
  unittest.main()
