#!/usr/bin/env python
"""
fruitPicker
1 Oct 2014
dent earl, dearl a soe ucsc edu

Script to scan through the mouse annotation project gene check pipeline
and pull out transcripts that may be interesting.
"""
from argparse import ArgumentParser
import itertools
import math
import sys
import os
sys.path.append(
  os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))), 'filters'))
import lib_filter


def initializeArguments(parser):
  parser.add_argument('--geneCheck', type=lib_filter.FileType)
  parser.add_argument('--geneCheckDetails', type=lib_filter.FileType)
  parser.add_argument('--label', choices=['nonsynon', 'synon', 'nonsense'],
                      help='single label to search for.')


def checkArguments(args, parser):
  # setting
  pairs = tuple((item, getattr(args, item)) for item in
                ['geneCheck', 'geneCheckDetails', 'label'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)


def containsOnlyLabel(transAnns, args):
  """ For a list of TranscriptAnnotation objects, return true if label
  only include the target label.
  """
  success = False
  for a in transAnns:
    if args.label not in a.labels:
      return False
    if args.label == a.labels[0] and len(a.labels) == 2:
      success = True
    elif args.label == 'nonsense' and args.label == a.labels[0]:
      success = True
    else:
      return False
  return success


def processTranscripts(args):
  """ Read the transcript bed files and return a dict of categories and counts.
  """
  transcripts = lib_filter.getTranscripts(args.geneCheck, args.geneCheckDetails)
  matches = []
  for t in transcripts:
    if t.annotations == []:
      continue
    if containsOnlyLabel(t.annotations, args):
      matches.append(t)
  return matches


def reportTranscripts(transcripts):
  """ Given list of transcripts, print them out
  """
  for t in transcripts:
    if t.chromosomeInterval.strand:
      strand = '+'
    else:
      strand = '-'
    try:
      print('%s %s %d %d %s %s' %
            (t.name, t.chromosomeInterval.chromosome,
             t.chromosomeInterval.start, t.chromosomeInterval.stop,
             strand,
             ', '.join(map(lambda a: str(list(a.labels)), t.annotations))))
    except IOError:
      # allow breaks for linux tools like head
      sys.exit(0)


def main():
  parser = ArgumentParser()
  initializeArguments(parser)
  args = parser.parse_args()
  checkArguments(args, parser)
  transcripts = processTranscripts(args)
  reportTranscripts(transcripts)


if __name__ == '__main__':
  main()
