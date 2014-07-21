#!/usr/bin/env python
"""
validateParalogs
dent earl, dearl a soe soe edu
21 july 2014

script to read through the name fields in bed files to collect information
about names that end in hyphen number (i.e. ENSMUST00000161846.1-3).
"""
from argparse import ArgumentParser
import os
import sys
sys.path.append(
  os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))), 'filters'))
import lib_filter


def initializeArguments(parser):
  parser.add_argument('--bed', type=lib_filter.FileType,
                      help='input bedfile to read')


def checkArguments(args, parser):
  # check for setting
  pairs = tuple((item, getattr(args, item)) for item in
                ['bed'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)


def bedIterator(infile):
  """ Provide an iterator for bed files.
  """
  while True:
    line = infile.readline().strip()
    if line == '':
      return
    data = line.split()
    yield data[3].split('/')[-1]


def readBed(filename):
  """ Given a FILENAME, open the bed file and return a list of names.
  """
  names = []
  with open(filename, 'r') as f:
    for line in bedIterator(f):
      names.append(line)
  return names


def validateData(names):
  """ Given a list of names DATA, ensure the data makes sense.
  """
  namesCount = {}
  for name in names:
    n = name.split('-')[0]
    c = namesCount.get(n, 0)
    namesCount[n] = c + 1
  for name in names:
    d = name.split('-')
    n = d[0]
    if len(d) > 1:
      if namesCount[n] <= 1:
        print n, namesCount[n]
      assert(namesCount[n] > 1)
  print 'names: %d' % len(names)
  print 'namesCount: %d' % len(namesCount)


def main():
  parser = ArgumentParser()
  initializeArguments(parser)
  args = parser.parse_args()
  checkArguments(args, parser)
  names = readBed(args.bed)
  validateData(names)


if __name__ == '__main__':
  main()
