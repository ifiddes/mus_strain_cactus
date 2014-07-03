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
import numpy
import lib_filter


def initializeArguments(parser):
  parser.add_argument('--geneCheck', type=lib_filter.FileType)
  parser.add_argument('--geneCheckDetails', type=lib_filter.FileType)
  parser.add_argument('--dontSplit', default=False, action='store_true',
                      help=('do not itemize labels, treat them as-is. '
                            'default=%(default)s'))
  parser.add_argument('--summaryCounts', default=False, action='store_true',
                      help=('put out a summary table of object counts at '
                            'the end.'))


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
          suffix = '+'
        else:
          suffix = ''
        labels.add('%s%s' % ('_'.join(a.labels[0:i]), suffix))
    else:
      labels.add(a.labels[0])
  return labels


def isOk(annots):
  """ given a list of TranscriptAnnotations, return True if data is OK.
  """
  if annots == []:
    return True
  for a in annots:
    for label in a.labels:
      if label not in ['hasOkCopy', 'hasBadCopy']:
        return False
  return True


def processTranscripts(args):
  """ Read the transcript bed files and return a dict of categories and counts.
  """
  transcripts = lib_filter.getTranscripts(args.geneCheck, args.geneCheckDetails)
  summary = {}
  summary['Transcripts'] = len(transcripts)
  summary['TranscriptAnnotations'] = 0
  summary['TA_per_T'] = []
  summary['L_per_TA'] = []
  categories = ['ok', 'not ok', 'badFrame', 'cdsGap', 'frameDiscontig',
                'frameMismatch', 'noStart', 'noStop', 'orfStop',
                'unknownCdsSplice', 'unknownUtrSplice', 'utrGap',
                'containsNs', 'alignmentPartialMap', 'alignmentAbutsEdge',
                'total']
  counts = {}
  for cat in categories:
    counts[cat] = 0
  for t in transcripts:
    counts['total'] += 1
    summary['TranscriptAnnotations'] += len(t.annotations)
    summary['TA_per_T'].append(len(t.annotations))
    if isOk(t.annotations):
      counts['ok'] += 1
      continue
    counts['not ok'] += 1
    labelSet = getAnnotationSet(t.annotations, args)
    summary['L_per_TA'].append(len(labelSet))
    for label in labelSet:
      if label not in counts:
        counts[label] = 0
      counts[label] += 1
  return counts, summary


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
  print '%*s %*d (%.3f)' % (max_str, 'not ok', max_digit, counts['not ok'],
                            float(counts['not ok']) / counts['total'])
  order = sorted(counts, key=lambda x: counts[x], reverse=True)
  for cat in order:
    if cat in ['total', 'ok', 'not ok']:
      continue
    print '%*s %*d (%.3f)' % (max_str, cat, max_digit, counts[cat],
                               float(counts[cat]) / counts['total'])


def reportSummary(summaryCounts, args):
  """ Given the summaryCounts object, print out a table of data.
  """
  if not args.summaryCounts:
    return
  print '##############################'
  print 'Total number of Transcripts: %d' % summaryCounts['Transcripts']
  print('Total number of TranscriptAnnotations: %d'
        % summaryCounts['TranscriptAnnotations'])
  print('%12s %10s %10s %10s %10s %10s'
        % ('Category', 'min', 'med', 'ave', 'max', 'std'))
  print('%12s %10d %10.1f %10.1f %10d %10.1f'
        % ('TA per T',
           numpy.min(summaryCounts['TA_per_T']),
           numpy.median(summaryCounts['TA_per_T']),
           numpy.average(summaryCounts['TA_per_T']),
           numpy.max(summaryCounts['TA_per_T']),
           numpy.std(summaryCounts['TA_per_T']),
           ))
  print('%12s %10d %10.1f %10.1f %10d %10.1f'
        % ('Label per TA',
           numpy.min(summaryCounts['L_per_TA']),
           numpy.median(summaryCounts['L_per_TA']),
           numpy.average(summaryCounts['L_per_TA']),
           numpy.max(summaryCounts['L_per_TA']),
           numpy.std(summaryCounts['L_per_TA']),
           ))


def main():
  parser = ArgumentParser()
  initializeArguments(parser)
  args = parser.parse_args()
  checkArguments(args, parser)
  counts, summaryCounts = processTranscripts(args)
  reportCounts(counts)
  reportSummary(summaryCounts, args)


if __name__ == '__main__':
  main()
