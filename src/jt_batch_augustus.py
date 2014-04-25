#!/usr/bin/env python
"""
jobTree based script that runs the augustus gene prediction pipeline.
"""
from argparse import ArgumentParser
from jobTree.scriptTree.stack import Stack
from jobTree.scriptTree.target import Target
from glob import glob
import lib_run
from math import floor
import os
from sonLib.bioio import logger
import subprocess
import sys
import time
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
    count = 0
    for start in xrange(self.args.window_start, self.args.window_end,
                        self.args.window_overlap):
      count += 1
      end = min(self.args.window_end, start + self.args.window_length)
      self.addChildTarget(AugustusCall(self.args.hal_file_path,
                                       self.args.ref_genome,
                                       self.args.ref_sequence, start,
                                       self.args.window_length, end,
                                       count - 1, self.args))
    logger.debug('There will be %d AugustusCall children' % count)
    self.args.batch_start_time = CreateSummaryReport(
      self.args.out_dir, self.args.ref_genome, self.args.ref_sequence,
      self.args.window_start, self.args.window_length, self.args.window_overlap,
      self.args.window_end, count, self.args.batch_start_time,
      self.args.calling_command)


class AugustusCall(Target):
  """
  AugustusCall class calls the augustus binary.
  """
  def __init__(self, hal_file_path, ref_genome, ref_sequence,
               window_start, window_length, window_end,
               window_number, args):
    Target.__init__(self)
    self.hal_file_path = hal_file_path
    self.ref_genome = ref_genome
    self.ref_sequence = ref_sequence
    self.window_start = window_start
    self.window_length = window_length
    self.window_end = window_end
    self.window_number = window_number
    self.out_path = os.path.join(args.out_dir,
                                 'window_%s_%s_%03d'
                                 % (args.ref_genome, args.ref_sequence,
                                    window_number))
    self.args = args
    dbaccess = ReadDBAccess(self.args.dbaccess_file)
    self.aug_parameters = {
      'AUGUSTUS_CONFIG_PATH': os.path.join(self.args.augustus_path, 'config'),
      'treefile': self.args.tree_path,
      'species': self.args.species,
      'dbaccess': dbaccess,
      'temperature': self.args.temperature,
      '/MeaPrediction/x0_E': self.args._MeaPrediction_x0_E,
      '/MeaPrediction/x0_I': self.args._MeaPrediction_x0_I,
      '/MeaPrediction/x1_E': self.args._MeaPrediction_x1_E,
      '/MeaPrediction/x1_I': self.args._MeaPrediction_x1_I,
      '/MeaPrediction/y0_E': self.args._MeaPrediction_y0_E,
      '/MeaPrediction/y0_I': self.args._MeaPrediction_y0_I,
      '/MeaPrediction/alpha_E': self.args._MeaPrediction_alpha_E,
      '/MeaPrediction/alpha_I': self.args._MeaPrediction_alpha_I,
      '/MeaPrediction/i1_E': self.args._MeaPrediction_i1_E,
      '/MeaPrediction/i1_I': self.args._MeaPrediction_i1_I,
      '/MeaPrediction/i2_E': self.args._MeaPrediction_i2_E,
      '/MeaPrediction/i2_I': self.args._MeaPrediction_i2_I,
      '/MeaPrediction/j1_E': self.args._MeaPrediction_j1_E,
      '/MeaPrediction/j1_I': self.args._MeaPrediction_j1_I,
      '/MeaPrediction/j2_E': self.args._MeaPrediction_j2_E,
      '/MeaPrediction/j2_I': self.args._MeaPrediction_j2_I,
      '/CompPred/ec_score': self.args._CompPred_ec_score,
      '/CompPred/ec_addend': self.args._CompPred_ec_addend,
      '/CompPred/ec_factor': self.args._CompPred_ec_factor,
      '/CompPred/dd_factor': self.args._CompPred_dd_factor,
      '/CompPred/exon_gain': self.args._CompPred_exon_gain,
      '/CompPred/exon_loss': self.args._CompPred_exon_loss,
      '/CompPred/phylo_factor': self.args._CompPred_phylo_factor,
      }
    if self.args._CompPred_only_species is not None:
      self.aug_parameters['/CompPred/only_species'] = (
        self.args._CompPred_only_species)

  def run(self):
    if not os.path.exists(self.out_path):
      os.mkdir(self.out_path)
    if self.args.maf_file_path is None:
      self.maf_file = os.path.join(self.getLocalTempDir(),
                                   'window_%s_%s_%03d.maf'
                                   % (self.args.ref_genome,
                                      self.args.ref_sequence,
                                      self.window_number))
    else:
      self.maf_file = self.args.maf_file_path
    self.aug_parameters['alnfile'] = self.maf_file
    # extract the region needed as maf
    hal2maf_cmd = [os.path.join(self.args.hal_path, 'bin', 'hal2maf')]
    hal2maf_cmd.append('--refGenome')
    hal2maf_cmd.append(self.args.ref_genome)
    hal2maf_cmd.append('--refSequence')
    hal2maf_cmd.append(self.args.ref_sequence)
    hal2maf_cmd.append('--ucscNames')
    hal2maf_cmd.append('--start')
    hal2maf_cmd.append(str(self.window_start))
    hal2maf_cmd.append('--length')
    hal2maf_cmd.append(str(self.window_length))
    hal2maf_cmd.append(self.args.hal_file_path)
    hal2maf_cmd.append(self.maf_file)
    hal2maf_cmds = [hal2maf_cmd]
    if self.args.maf_file_path is None:
      time_start = TimeStamp(self.out_path)
      LogCommand(self.out_path, hal2maf_cmds)
      if not self.args.debug:
        lib_run.RunCommandsS(hal2maf_cmds, self.getLocalTempDir())
        TimeStamp(self.out_path, time_start)
    # run augustus on the maf
    err_pipe = [os.path.join(self.getLocalTempDir(), 'stderr.out')]
    out_pipe = [os.path.join(self.getLocalTempDir(), 'stdout.out')]
    aug_cmd = [os.path.join(self.args.augustus_path, 'bin', 'augustus')]
    for key in self.aug_parameters:
      aug_cmd.append('--%s=%s' % (key, str(self.aug_parameters[key])))
    aug_cmds = [aug_cmd]
    LogCommand(self.out_path, aug_cmds, out_pipe=out_pipe,
               err_pipe=err_pipe)
    if not self.args.debug:
      time_start = TimeStamp(self.out_path)
      lib_run.RunCommandsS(aug_cmds, self.getLocalTempDir(),
                           out_pipes=out_pipe, err_pipes=err_pipe)
      TimeStamp(self.out_path, time_start)
    # copy output files from tmp back to the target dir
    copy_cmds = []
    # todo: copy out actual results
    for suffix in ['dot', 'gff', 'gff3', 'out', 'wig']:
      files = glob(os.path.join(self.getLocalTempDir(), '*.%s' % suffix))
      for f in files:
        copy_cmds.append([lib_run.Which('cp'), f, os.path.join(self.out_path)])
    LogCommand(self.out_path, copy_cmds)
    if not self.args.debug:
      time_start = TimeStamp(self.out_path)
      # we could use RunCommandsP here, but we might hammer the disk if we did.
      lib_run.RunCommandsS(copy_cmds, self.getLocalTempDir())
      TimeStamp(self.out_path, time_start)


def TimeStamp(out_path, time_start=None):
  """ Open up the log file and make a timestamp.
  """
  now = time.time()
  f = open(os.path.join(out_path, 'jt_issued_commands.log'), 'a')
  if time_start is not None:
    elapsed_time = now - time_start
    f.write('[%s] End (elapsed: %s)\n' %
            (time.strftime("%a, %d %b %Y %H:%M:%S (%Z)", time.localtime(now)),
             PrettyTime(elapsed_time)))
  else:
    f.write('[%s] Start\n' % (time.strftime("%a, %d %b %Y %H:%M:%S (%Z)",
                                            time.localtime(now))))
  f.close()


def ReadDBAccess(dbaccess_file):
  """ Open the file and read the first line, report it back.
  """
  f = open(dbaccess_file, 'r')
  line = f.read()
  line = line.strip()
  return line


def LogCommand(out_path, cmds, out_pipe=None, err_pipe=None):
  """ Write out the commands that will be executed for this run.
  """
  f = open(os.path.join(out_path, 'jt_issued_commands.log'), 'a')
  if out_pipe is None:
    out_str = ''
  else:
    out_str = ' 1>%s' % ' '.join(out_pipe)
  if err_pipe is None:
    err_str = ''
  else:
    err_str = ' 2>%s' % ' '.join(err_pipe)
  for c in cmds:
    f.write('[%s] %s%s%s\n' % (time.strftime("%a, %d %b %Y %H:%M:%S (%Z)",
                                             time.localtime(time.time())),
                               ' '.join(c),
                               err_str, out_str))
  f.close()


def CreateSummaryReport(out_dir, ref_genome, ref_sequence, window_start,
                        window_length, window_overlap, window_end, count,
                        now, command):
  """ Create a summary report in the root output directory.
  """
  f = open(os.path.join(out_dir, 'summary_report.txt'), 'w')
  f.write('run started: %s\n' % time.strftime("%a, %d %b %Y %H:%M:%S (%Z)",
                                              time.localtime(now)))
  f.write('command: %s\n' % command)
  f.write('window start:   %d\n' % window_start)
  f.write('window length:  %d\n' % window_length)
  f.write('window overlap: %d\n' % window_overlap)
  f.write('window end:     %d\n' % window_end)
  f.write('region length:  %d\n' % (window_end - window_start))
  f.write('num windows:    %d\n' % count)
  f.close()
  return now


def InitializeArguments(parser):
  logger.debug('Initializing arguments')
  parser.add_argument('--augustus_path', type=str,
                      help='location of augustus directory.')
  parser.add_argument('--hal_path', type=str,
                      help='location of hal tools directory.')
  parser.add_argument('--hal_file_path', type=str,
                      help='location hal file.')
  parser.add_argument('--tree_path', type=str,
                      help='location newick tree file.')
  parser.add_argument('--out_dir', type=str,
                      help='location to store output files.')
  parser.add_argument('--dbaccess_file', type=str,
                      help='location of dbaccess file containing login info.')
  parser.add_argument('--maf_file_path', type=str,
                      help=('location maf file. Overrides all hal window '
                            'extraction'))
  parser.add_argument('--debug', default=False, action='store_true',
                      help='turns off execution of commands, just writes logs.')
  window = parser.add_argument_group('Window options')
  window.add_argument('--ref_genome', type=str,
                      help='reference genome to use for region extraction.')
  window.add_argument('--ref_sequence', type=str,
                      help=('reference sequence (chr) to use for '
                            'region extraction.'))
  window.add_argument('--window_start', type=int, default=0,
                      help='start of windowing region.')
  window.add_argument('--window_end', type=int, default=None,
                      help='end of windowing region.')
  window.add_argument('--window_length', type=int, default=2000000,
                      help='length of each window. default=%(default)s')
  window.add_argument('--window_overlap', type=int, default=1000000,
                      help='overlap of each window. default=%(default)s')
  augustus = parser.add_argument_group('Augustus options')
  augustus.add_argument('--species', default='human', type=str,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--temperature', default=3, type=int,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/x0_E', default=-1.25,
                        dest='_MeaPrediction_x0_E', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/x0_I', default=-0.78125,
                        dest='_MeaPrediction_x0_I', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/x1_E', default=5,
                        dest='_MeaPrediction_x1_E', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/x1_I', default=10,
                        dest='_MeaPrediction_x1_I', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/y0_E', default=0.5,
                        dest='_MeaPrediction_y0_E', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/y0_I', default=0.9,
                        dest='_MeaPrediction_y0_I', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/alpha_E', default=9.375,
                        dest='_MeaPrediction_alpha_E', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/alpha_I', default=2.5075,
                        dest='_MeaPrediction_alpha_I', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/i1_E', default=0.,
                        dest='_MeaPrediction_i1_E', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/i1_I', default=0.,
                        dest='_MeaPrediction_i1_I', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/i2_E', default=0.5,
                        dest='_MeaPrediction_i2_E', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/i2_I', default=0.9,
                        dest='_MeaPrediction_i2_I', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/j1_E', default=-1.25,
                        dest='_MeaPrediction_j1_E', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/j1_I', default=-0.78125,
                        dest='_MeaPrediction_j1_I', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/j2_E', default=0.,
                        dest='_MeaPrediction_j2_E', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/MeaPrediction/j2_I', default=0.,
                        dest='_MeaPrediction_j2_I', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/CompPred/ec_score', default=-13,
                        dest='_CompPred_ec_score', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/CompPred/ec_addend', default=-20.6,
                        dest='_CompPred_ec_addend', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/CompPred/ec_factor', default=6,
                        dest='_CompPred_ec_factor', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/CompPred/dd_factor', default=20,
                        dest='_CompPred_dd_factor', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/CompPred/exon_gain', default=0.0001,
                        dest='_CompPred_exon_gain', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/CompPred/exon_loss', default=0.0001,
                        dest='_CompPred_exon_loss', type=float,
                        help='Augustus option. default=%(default)s')
  augustus.add_argument('--/CompPred/phylo_factor', default=50.,
                        dest='_CompPred_phylo_factor', type=float,
                        help='Augustus option. default=%(default)s')
  #####
  augustus.add_argument('--/CompPred/only_species', default=None,
                        dest='_CompPred_only_species', type=str,
                        help='Augustus option. default=%(default)s')


def CheckArguments(args, parser):
  # check for setting
  for name, value in [('augustus_path', args.augustus_path),
                      ('hal_path', args.hal_path),
                      ('hal_file_path', args.hal_path),
                      ('tree_path', args.tree_path),
                      ('out_dir', args.out_dir),
                      ('dbaccess_file', args.dbaccess_file)]:
    if value is None:
      parser.error('Specify --%s' % name)
    else:
      value = os.path.abspath(value)
  # check for existence
  for name, value in [('augustus_path', args.augustus_path),
                      ('hal_path', args.hal_path),
                      ('hal_file_path', args.hal_file_path),
                      ('tree_path', args.tree_path),
                      ('dbaccess_file', args.dbaccess_file),
                      ]:
    if not os.path.exists(value):
      parser.error('--%s %s does not exist' % (name, value))
  if not os.path.exists(args.out_dir):
    os.mkdir(args.out_dir)
  # check for directories
  for name, value in [('augustus_path', args.augustus_path),
                      ('hal_path', args.hal_path),
                      ('out_dir', args.out_dir)]:
    if not os.path.isdir(value):
      parser.error('--%s %s is not a directory' % (name, value))
  # check for files
  for value in [os.path.join(args.augustus_path, 'bin', 'augustus'),
                os.path.join(args.hal_path, 'bin', 'hal2maf'),
                args.hal_file_path,
                args.tree_path,
                args.dbaccess_file,
                ]:
    if not os.path.isfile(value):
      parser.error('%s is not a file' % value)
  # check for executability
  for value in [(os.path.join(args.augustus_path, 'bin', 'augustus')),
                      (os.path.join(args.hal_path, 'bin', 'hal2maf'))]:
    if not os.access(value, os.X_OK):
      parser.error('%s is not executable' % value)
  args.augustus_path = os.path.abspath(args.augustus_path)
  args.hal_path = os.path.abspath(args.hal_path)
  args.hal_file_path = os.path.abspath(args.hal_file_path)
  args.tree_path = os.path.abspath(args.tree_path)
  args.out_dir = os.path.abspath(args.out_dir)
  logger.debug('Arguments checked.\n'
               'augustus_path:%s\n'
               'hal_path:%s\n'
               'hal_file_path:%s\n'
               'tree_path:%s\n'
               'out_dir:%s\n'
               % (args.augustus_path, args.hal_path, args.hal_file_path,
                  args.tree_path, args.out_dir))
  args.calling_command = '%s' % ' '.join(sys.argv[0:])


def PrettyTime(t):
  """ Given input t as seconds, return a nicely formated string.
  """
  plural_dict = {True: 's', False: ''}
  if t < 120:
    return '%ds' % t
  if t < 120 * 60:
    m = floor(t / 60.)
    s = t % 60
    return '%dm %ds' % (m, s)
  if t < 25 * 60 * 60:
    h = floor(t / 60. / 60.)
    m = floor((t - (h * 60. * 60.)) / 60.)
    s = t % 60
    return '%dh %.0fm %ds' % (h, m, s)
  if t < 7 * 24 * 60 * 60:
    d = floor(t / 24. / 60. / 60.)
    h = floor((t - (d * 24. * 60. * 60.)) / 60. / 60.)
    m = floor((t
               - (d * 24. * 60. * 60.)
               - (h * 60. * 60.))
              / 60.)
    s = t % 60
    d_plural = plural_dict[d > 1]
    return '%d day%s %dh %dm %ds' % (d, d_plural, h, m, s)
  w = floor(t / 7. / 24. / 60. / 60.)
  d = floor((t - (w * 7 * 24 * 60 * 60)) / 24. / 60. / 60.)
  h = floor((t
             - (w * 7. * 24. * 60. * 60.)
             - (d * 24. * 60. * 60.))
            / 60. / 60.)
  m = floor((t
             - (w * 7. * 24. * 60. * 60.)
             - (d * 24. * 60. * 60.)
             - (h * 60. * 60.))
            / 60.)
  s = t % 60
  w_plural = plural_dict[w > 1]
  d_plural = plural_dict[d > 1]
  return '%d week%s %d day%s %dh %dm %ds' % (w, w_plural, d,
                                             d_plural, h, m, s)


def LaunchBatch(args):
  args.batch_start_time = time.time()
  jobResult = Stack(BatchJob(args)).startJobTree(args)
  if jobResult:
    raise RuntimeError('The jobTree contained %d failed jobs!\n' % jobResult)
  f = open(os.path.join(args.out_dir, 'summary_report.txt'), 'a')
  now = time.time()
  f.write('run finished: %s\n' %
          time.strftime("%a, %d %b %Y %H:%M:%S (%Z)", time.localtime(now)))
  f.write('elapsed time: %s\n' % PrettyTime(now - args.batch_start_time))
  f.close()


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
