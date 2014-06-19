#!/usr/bin/env python
"""
result2bedDirs
19 June 2014
dent earl, dearl a soe ucsc edu

Script to take in locations of many results of
msca project filtering pipeline and kick out a
bed directory for use with the hal2assemblyHub
script.

"""
from argparse import ArgumentParser, FileType
import os
import shutil
import sys
sys.path.append(
  os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))), 'filters'))
import lib_filter


class Result(object):
  def __init__(self, d):
    self.original = d
    self.track = dirToTrack(d)
    self.genome = dirToGenome(d)
    self.bed = os.path.join(d, 'out.bed')
    self.bedDetails = os.path.join(d, 'out_details.bed')


def initializeArguments(parser):
  parser.add_argument('resultDirs', nargs='*', type=lib_filter.DirType)
  parser.add_argument('--outDir', type=lib_filter.DirType)


def checkArguments(args, parser):
  # setting
  pairs = tuple((item, getattr(args, item)) for item in ['outDir'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)


def dirToTrack(d):
  """ Given the path to a pipeline result directory, return the track name.
  """
  return os.path.basename(d).split('.')[0]


def dirToGenome(d):
  """ Given the path to a pipeline result directory, return the genome name.
  """
  return os.path.basename(d).split('.')[-1]


def getResults(args):
  """ create a Result for each result directory
  """
  results = []
  for d in args.resultDirs:
    results.append(Result(d))
  return results


def copyResults(results, args):
  """ create all the necessary directories and make copies of data.
  """
  if results == []:
    return
  for r in results:
    destDirs = [os.path.join(args.outDir, r.track),
                os.path.join(args.outDir, r.track + '_detail'),]
    for d in destDirs:
      if not os.path.exists(d):
        os.mkdir(d)
      if not os.path.exists(os.path.join(d, r.genome)):
        os.mkdir(os.path.join(d, r.genome))
    bedDir = os.path.join(destDirs[0], r.genome)
    bedDetailsDir = os.path.join(destDirs[1], r.genome)
    shutil.copy2(r.bed, os.path.join(bedDir, '%s.%s.bed'
                                     % (r.genome, r.track)))
    shutil.copy2(r.bedDetails, os.path.join(bedDetailsDir, '%s.%s.bed'
                                            % (r.genome, r.track)))


def main():
  parser = ArgumentParser()
  initializeArguments(parser)
  args = parser.parse_args()
  checkArguments(args, parser)
  results = getResults(args)
  copyResults(results, args)


if __name__ == '__main__':
  main()
