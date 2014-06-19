#!/usr/bin/env python
"""
geneCheckStats
11 June 2014
dent earl, dearl a soe ucsc edu

Script to scan through geneCheck formatted bed files and output
statistics on the labels contained in annotation rows.
"""
import sys
import os
sys.path.append(
  os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))), 'filters'))
from argparse import ArgumentParser
import math
import lib_filter


def initializeArguments(parser):
  parser.add_argument('--geneCheck', type=lib_filter.FileType)
  parser.add_argument('--geneCheckDetails', type=lib_filter.FileType)
  parser.add_argument('--dontSplit', default=False, action='store_true',
                      help=('do not itemize labels, treat them as-is. '
                            'default=%(default)s'))


def checkArguments(args, parser):
  # setting
  pairs = tuple((item, getattr(args, item)) for item in
                ['geneCheck', 'geneCheckDetails'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)


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
          suffix = '*'
        else:
          suffix = ''
        labels.add('%s%s' % ('_'.join(a.labels[0:i]), suffix))
    else:
      labels.add(a.labels[0])
  return labels


def processTranscripts(args):
  """ Read the transcript bed files and return a dict of categories and counts.
  """
  transcripts = lib_filter.getTranscripts(args.geneCheck, args.geneCheckDetails)
  categories = ['ok', 'any', 'badFrame', 'cdsGap', 'frameDiscontig',
                'frameMismatch', 'noStart', 'noStop', 'orfStop',
                'unknownCdsSplice', 'unknownUtrSplice', 'utrGap',
                'containsNs', 'alignmentPartialMap', 'alignmentAbutsEdge',
                'total']
  counts = {}
  for cat in categories:
    counts[cat] = 0
  for t in transcripts:
    counts['total'] += 1
    if t.annotations == []:
      counts['ok'] += 1
      continue
    counts['any'] += 1
    labelSet = getAnnotationSet(t.annotations, args)
    for label in labelSet:
      if label not in counts:
        counts[label] = 0
      counts[label] += 1
  return counts


def reportCounts(counts):
  """ Given a dict of counts (key: category, value: integer count, print it out.
  """
  max_str = 0
  max_digit = 0
  for c in counts:
    if max_str <= len(c):
      max_str = len(c) + 1
    try:
      v = round(math.log10(counts[c]))
    except ValueError:
      v = 1
    if max_digit < v:
      max_digit = int(v)
  print '%*s %*d (%.3f)' % (max_str, 'total', max_digit, counts['total'], 1.0)
  print '%*s %*d (%.3f)' % (max_str, 'ok', max_digit, counts['ok'],
                            float(counts['ok']) / counts['total'])
  print '%*s %*d (%.3f)' % (max_str, 'any', max_digit, counts['any'],
                            float(counts['any']) / counts['total'])
  order = sorted(counts, key=lambda x: counts[x], reverse=True)
  for cat in order:
    if cat in ['total', 'ok', 'any']:
      continue
    print '%*s %*d (%.3f)' % (max_str, cat, max_digit, counts[cat],
                               float(counts[cat]) / counts['total'])


def main():
  parser = ArgumentParser()
  initializeArguments(parser)
  args = parser.parse_args()
  checkArguments(args, parser)
  counts = processTranscripts(args)
  reportCounts(counts)


if __name__ == '__main__':
  main()
