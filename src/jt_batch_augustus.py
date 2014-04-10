#!/usr/bin/env python
"""
jobTree based script that runs the augustus gene prediction pipeline.
"""
from argparse import ArgumentParser
from jobTree.scriptTree.stack import Stack
from jobTree.scriptTree.target import Target
from glob import glob
import lib_run
import os
from sonLib.bioio import logger
import subprocess
import sys
##################################################
# Copyright (c) 2014 Dent Earl, Benedict Paten, Mark Diekhans, Craig Hunter
# ... and other members of the Reconstruction Team of David Haussler's
# lab (BME Dept. UCSC)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##################################################


class BatchJob(Target):
  """
  The BatchJob class runs the entire batch. It creates as many
  AugustusCall() instances as required by the input.
  """
  def __init__(self, args):
    Target.__init__(self)
    self.args = args
  def run(self):
    if self.args.tc:
      flavor = 'tc'
    else:
      flavor = 'orig'
    path = os.path.join(self.args.dir, 'predictions', '*.%s.maf' % flavor)
    logger.debug('globbing for %s' % path)
    mafs = glob(path)
    count = 0
    for m in mafs:
      if not os.path.isfile(m):
        continue
      for s in species:
        count += 1
        self.addChildTarget(AugustusCall(m, s, self.args))
    logger.debug('There will be %d AugustusCall children' % count)


class AugustusCall(Target):
  """
  AugustusCall class calls the augustus binary
  """
  def __init__(self, targetMaf, species, args):
    Target.__init__(self)
    self.targetMaf = targetMaf
    self.species = species
    self.out_name = os.path.basename(self.targetMaf) + '.%s.txt' % self.species
    self.out_path = os.path.join(args.out_dir, self.out_name)
    self.args = args
  def run(self):
    cmds = []
    outpipes = []
    # errpipes = []
    followCopy = [] # copy output files from tmp back to the target dir
    followJob = []
    copyPath = lib_run.which('cp')
    followCopy.append([copyPath,
                       os.path.join(self.getLocalTempDir(), self.out_name),
                       self.out_path])
    cmd = [self.args.bin, '--maf=%s' % self.targetMaf,
           '--species=%s' % self.species]
    if self.args.ncoverage:
      cmd.append('--nCoverage')
    cmds.append(cmd)
    outpipes.append(os.path.join(self.getLocalTempDir(), self.out_name))
    # errpipes.append(# put stuff in here...)
    lib_run.runCommandsS(cmds, self.getLocalTempDir(), outPipes=outpipes)
    if (os.path.exists(os.path.join(self.getLocalTempDir(), self.out_name))):
      lib_run.runCommandsS(followCopy, self.getLocalTempDir())


def InitializeArguments(parser):
  logger.debug('Initializing arguments')
  parser.add_argument('--dir', type=str, help=('path to packageXXX/'))
  parser.add_argument('--bin', type=str,
                      help='location of augustus executable.')
  parser.add_argument('--mode', type=str,
                      help='must be primates, mammals, or flies')
  parser.add_argument('--out_dir', type=str,
                      help='location to store output coverage files.')
  parser.add_argument('--ncoverage', default=False, action='store_true',
                      help='turn on ncoverage mode')
  parser.add_argument('--tc', default=False, action='store_true',
                      help='uses transitive closure (tc) mafs instead of orig')


def CheckArguments(args, parser):
  for name, value in [('bin', args.bin), ('dir', args.dir),
                      ('mode', args.mode), ('out_dir', args.out_dir)]:
    if value is None:
      parser.error('Specify --%s' % name)
  if not os.path.exists(args.bin):
    parser.error('--bin %s does not exist' % args.bin)
  if not os.access(args.bin, os.X_OK):
    parser.error('--bin %s is not executable' % args.bin)
  if not os.path.exists(args.dir):
    parser.error('--dir %s does not exist' % args.dir)
  if not os.path.exists(args.out_dir):
    parser.error('--out_dir %s does not exist' % args.out_dir)
  if not os.path.isdir(args.dir):
    parser.error('--dir %s is not a directory' % args.dir)
  if not os.path.isdir(args.out_dir):
    parser.error('--out_dir %s is not a directory' % args.out_dir)
  if args.mode not in ['primates', 'mammals', 'flies']:
    parser.error('Unrecognized --mode %s, must be primates mammals or flies'
                 % args.mode)
  args.bin = os.path.abspath(args.bin)
  args.dir = os.path.abspath(args.dir)
  args.out_dir = os.path.abspath(args.out_dir)
  logger.debug('Arguments checked.\nbin:%s\ndir:%s\n'
               'out_dir:%s\n' % (args.bin, args.dir, args.out_dir))


def LaunchBatch(args):
  jobResult = Stack(BatchJob(args)).startJobTree(args)
  if jobResult:
    raise RuntimeError('The jobTree contained %d failed jobs!\n' % jobResult)


def main():
  description = ('%(prog)s starts a jobTree batch of augustus calls '
                 'for a given input set.')
  parser = ArgumentParser(description=description)
  InitializeArguments(parser)
  Stack.addJobTreeOptions(parser)
  args = parser.parse_args()
  CheckArguments(args, parser)
  LaunchBatch(args)


if __name__ == '__main__':
    from jt_batch_augustus import *
    main()
