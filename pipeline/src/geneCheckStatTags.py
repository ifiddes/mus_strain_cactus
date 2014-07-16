#!/usr/bin/env python
"""
geneCheckStatCreator
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
import lib_filter
import lib_stat_graph as lsg


def initializeArguments(parser):
  parser.add_argument('--xml', type=lib_filter.FileType)
  parser.add_argument('--tag', type=str,
                      help='gather stats on a specific tag')
  parser.add_argument('--tagLowerBound', type=int,
                      help='filter out all counts below this value')


def checkArguments(args, parser):
  # check for setting
  pairs = tuple((item, getattr(args, item)) for item in
                ['xml', 'tag'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)


def main():
  parser = ArgumentParser()
  initializeArguments(parser)
  args = parser.parse_args()
  checkArguments(args, parser)
  graph = lsg.readStatGraph(args.xml)
  stats = lsg.getTagStats(graph, args.tag)
  lsg.reportTagStats(stats, args.tagLowerBound)


if __name__ == '__main__':
  main()
