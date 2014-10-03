#!/usr/bin/env python
"""
getTranscriptsByLabel
20 June 2014
dent earl, dearl a soe ucsc edu

Script to scan through the mouse annotation project gene check pipeline
and pull out transcripts by their labels.
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
  parser.add_argument('--labels',
                      help='comma separated list of labels to search for')
  parser.add_argument('--dontSplit', default=False, action='store_true',
                      help=('do not itemize labels, treat them as-is. '
                            'default=%(default)s'))


def checkArguments(args, parser):
  # setting
  pairs = tuple((item, getattr(args, item)) for item in
                ['geneCheck', 'geneCheckDetails', 'labels'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)
  args.labels = args.labels.split(',')


def getAnnotationSet(transAnns, args):
  """ For a list of TranscriptAnnotation objects, return a set that contains
  every single label from the TranscriptAnnotation.labels lists along with all
  suffixes for multi-item lists.
  """
  labels = set()
  for a in transAnns:
    if args.dontSplit:
      labels.add('_'.join(a.labels))
      continue
    if len(a.labels) > 1:
      # all suffixes
      for i in xrange(1, len(a.labels) + 1):
        if i < len(a.labels):
          suffix = '+'
        else:
          suffix = ''
        labels.add('%s%s' % ('_'.join(list(itertools.islice(a.labels, 0, i))), suffix))
    else:
      labels.add(a.labels[0])
    printed = False
    for l in a.labels:
      if l.startswith('nonsynon') and not printed:
        print a.labels
        printed = True
  return labels


def processTranscripts(args):
  """ Read the transcript bed files and return a dict of categories and counts.
  """
  transcripts = lib_filter.getTranscripts(args.geneCheck, args.geneCheckDetails)
  matches = []
  print 'labels: %s' % str(args.labels)
  for t in transcripts:
    if t.annotations == []:
      continue
    labelSet = getAnnotationSet(t.annotations, args)
    # check to see if there is any overlap between the two lists
    match = False
    match = reduce(lambda x, y: match or x or y,
                   (label in labelSet for label in args.labels))
    if match:
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
