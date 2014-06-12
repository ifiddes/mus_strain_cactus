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
    transcriptBedLines = ['1       2812346 3113743 ENSMUST00000178026.1    0       -       2812370 3038729 128,0,0 9       54,2,89,249,90,165,105,13,45    0,58,62,698,1209,1305,226292,301050,301352',
                          '1       2812346 3113783 ENSMUST00000095795.4    0       +       2812370 3038729 128,0,0 9       54,2,89,249,197,52,105,13,85    0,58,62,698,1209,1418,226292,301050,301352',
                          'scaffold-100021  466  4248  ENSMUST00000034053.5  0  -  466  4248  128,0,0  2  85,152  0,3630',
                          'scaffold-138877  4903  5091  ENSMUST00000034053.5  0  -  4903  4996  128,0,0  1  188  0',
                          'scaffold-2051  13759  24866  ENSMUST00000034053.5  0  -  14291  24866  128,0,0  4  722,112,131,188  0,1977,4316,10919',
                          ]
    transcriptDetailsBedLines = ['1       2812346 2812349 noStop/ENSMUST00000095795.4',
                                 '1       3113780 3113783 noStart/ENSMUST00000095795.4',
                                 'scaffold-100021  466  469  noStop/ENSMUST00000034053.5',
                                 'scaffold-100021  4245  4248  noStart/ENSMUST00000034053.5',
                                 'scaffold-138877  4903  4906  noStop/ENSMUST00000034053.5',
                                 'scaffold-2051  24863  24866  noStart/ENSMUST00000034053.5',
                                 ]
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
      print 'test: %d' % i
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
    transcriptBedLines = ['1       2812346 3113743 ENSMUST00000178026.1    0       -       2812370 3038729 128,0,0 9       54,2,89,249,90,165,105,13,45    0,58,62,698,1209,1305,226292,301050,301352',
                          '1       2812346 3113783 ENSMUST00000095795.4    0       +       2812370 3038729 128,0,0 9       54,2,89,249,197,52,105,13,85    0,58,62,698,1209,1418,226292,301050,301352',
                          'scaffold-100021  466  4248  ENSMUST00000034053.5  0  -  466  4248  128,0,0  2  85,152  0,3630',
                          'scaffold-138877  4903  5091  ENSMUST00000034053.5  0  -  4903  4996  128,0,0  1  188  0',
                          'scaffold-2051  13759  24866  ENSMUST00000034053.5  0  -  14291  24866  128,0,0  4  722,112,131,188  0,1977,4316,10919',
                          ]
    transcriptDetailsBedLines = ['1       2812370 2812372 noStop/ENSMUST00000095795.4',
                                 '1       2812370 2812372 noStart/noStop/ENSMUST00000095795.4',
                                 'scaffold-100021  466  469  noStop/ENSMUST00000034053.5',
                                 'scaffold-100021  4245  4248  noStart/ENSMUST00000034053.5',
                                 'scaffold-138877  4903  4906  noStop/ENSMUST00000034053.5',
                                 'scaffold-2051  24863  24866  noStart/ENSMUST00000034053.5',
                                 ]
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


if __name__ == '__main__':
  unittest.main()
