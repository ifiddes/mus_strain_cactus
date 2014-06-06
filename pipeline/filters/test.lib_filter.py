""" Test the lib_filter classes and functions
"""
import os
import shutil
import string
import subprocess
import unittest
import lib_filter


def makeTempDirParent():
  """ make the parent temp dir directory
  """
  if not os.path.exists(os.path.join(os.curdir, 'tempTestDir')):
    os.mkdir(os.path.join(os.curdir, 'tempTestDir'))


def removeTempDirParent():
  """ remove the parent temp dir directory
  """
  if os.path.exists(os.path.join(os.curdir, 'tempTestDir')):
    shutil.rmtree(os.path.join(os.curdir, 'tempTestDir'))


def makeTempDir(name=None):
  """ make the directory where all temporary test files will be stored.
  """
  makeTempDirParent()
  charSet = string.ascii_lowercase + '123456789'
  if name is None:
    while True:
      name = '%s_%s' % (''.join(random.choice(charSet) for x in xrange(4)),
                        ''.join(random.choice(charSet) for x in xrange(4)))
      if not os.path.exists(os.path.join(os.curdir, 'tempTestDir', name)):
        break
  if not os.path.exists(os.path.join(os.curdir, 'tempTestDir', name)):
    os.mkdir(os.path.join(os.curdir, 'tempTestDir', name))
  return os.path.join(os.curdir, 'tempTestDir', name)


def removeDir(dirpath):
  """ destroy a directory
  """
  if os.path.exists(dirpath):
    shutil.rmtree(dirpath)


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
    for path in os.environ["PATH"].split(os.pathsep):
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
  """ given a list of tuples return path to temp file.
  """
  seqfile = os.path.join(tmpDir, 'seq.fa')
  with open(seqfile, 'w') as f:
    for name, seq in sequences:
      f.write('>%s\n%s' % (name, seq))
  return seqfile


class sequenceGetterTests(unittest.TestCase):
  def test_getSequences(self):
    """ getSequences must read a fasta and return a dict of Sequence objects.
    """
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('getSequences'))
    testFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(testFile)
    for name, seq in sequences:
      seq = seq.replace('\n', '').strip()
      self.assertTrue(name in seqDict)
      self.assertTrue(seqDict[name].getLength() == len(seq))
      self.assertTrue(seqDict[name].getSequence() == seq)
    removeDir(tmpDir)


sequences = [('chrA',
              'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n'
              'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n'
              'ACGT'
              ),
             ]

class transcriptIteratorTests(unittest.TestCase):
  def test_getSequences(self):
    """ tests the transcriptIterator function
    """
    transcriptBedLines = ["1       2812346 3113743 ENSMUST00000178026.1    0       -       2812370 3038729 128,0,0 9       54,2,89,249,90,165,105,13,45    0,58,62,698,1209,1305,226292,301050,301352",
                          "1       2812346 3113783 ENSMUST00000095795.4    0       +       2812370 3038729 128,0,0 9       54,2,89,249,197,52,105,13,85    0,58,62,698,1209,1418,226292,301050,301352"]
    transcriptDetailsBedLines = ["1       2812370 2812372 noStop/ENSMUST00000095795.4",
                                 "1       2812370 2812372 noStart/noStop/ENSMUST00000095795.4"]
    transcripts = [ transcript for transcript in lib_filter.transcriptIterator(transcriptBedLines, transcriptDetailsBedLines) ]
    self.assertTrue(len(transcripts) == 2)
    transcript1, transcript2 = transcripts
    self.assertTrue(transcript1.name == "ENSMUST00000178026.1")
    self.assertTrue(transcript2.name == "ENSMUST00000095795.4")
    self.assertTrue(transcript1.chromosomeInterval.chromosome == "1")
    self.assertTrue(transcript1.chromosomeInterval.start == 2812346)
    self.assertTrue(transcript1.chromosomeInterval.stop == 3113743)
    self.assertTrue(transcript1.chromosomeInterval.strand == False)
    self.assertTrue(len(transcript1.exons) == 9)
    self.assertTrue(transcript1.exons[0].chromosome == "1")
    self.assertTrue(transcript1.exons[0].start == 2812346)
    self.assertTrue(transcript1.exons[0].stop == 2812346 + 54)
    self.assertTrue(transcript1.exons[0].strand == False)
    self.assertTrue(transcript2.exons[0].chromosome == "1")
    self.assertTrue(transcript2.exons[0].start == 2812346)
    self.assertTrue(transcript2.exons[0].stop == 2812346 + 54)
    self.assertTrue(transcript2.exons[0].strand == True)
    self.assertTrue(transcript1.annotations == [])
    self.assertTrue(len(transcript2.annotations) == 2)
    self.assertTrue(transcript2.annotations[0].name == "ENSMUST00000095795.4")
    self.assertTrue(transcript2.annotations[1].name == "ENSMUST00000095795.4")
    self.assertTrue(transcript2.annotations[0].annotation == [ "noStop" ])
    self.assertTrue(transcript2.annotations[1].annotation == [ "noStart", "noStop" ])
    self.assertTrue(transcript2.annotations[0].chromosomeInterval.chromosome == "1")
    self.assertTrue(transcript2.annotations[0].chromosomeInterval.start == 2812370)
    self.assertTrue(transcript2.annotations[0].chromosomeInterval.stop == 2812372)
    self.assertTrue(transcript2.annotations[0].chromosomeInterval.strand == None)
    #Bed fields
    self.assertTrue(transcript1.thickStart == 2812370)
    self.assertTrue(transcript1.thickEnd == 3038729)
    self.assertTrue(transcript1.itemRgb == "128,0,0")
    #Check print functions
    self.assertEquals(transcript1.bedString().split(), transcriptBedLines[0].split())
    self.assertEquals(transcript2.bedString().split(), transcriptBedLines[1].split())
    self.assertEquals(transcript2.annotations[0].bedString().split(), transcriptDetailsBedLines[0].split())
    self.assertEquals(transcript2.annotations[1].bedString().split(), transcriptDetailsBedLines[1].split())
    #Check sort function for transcripts
    transcripts.reverse()
    self.assertEquals(transcripts[0].name, "ENSMUST00000095795.4")
    self.assertEquals(transcripts[1].name, "ENSMUST00000178026.1")
    transcripts.sort()
    self.assertEquals(transcripts[1].name, "ENSMUST00000095795.4")
    self.assertEquals(transcripts[0].name, "ENSMUST00000178026.1")


if __name__ == '__main__':
  try:
    unittest.main()
  finally:
    removeTempDirParent()  # cleanup
