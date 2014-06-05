"""
convenience library for assisting filters.
"""
import os


def InitializeArguments(parser):
  parser.add_argument('--refGenome')
  parser.add_argument('--genome')
  parser.add_argument('--geneCheckBed')
  parser.add_argument('--geneCheckBedDetails')
  parser.add_argument('--alignment')
  parser.add_argument('--sequence')
  parser.add_argument('--chromSizes')
  parser.add_argument('--outDir')


def CheckArguments(args, parser):
  # setting
  pairs = tuple((item, getattr(args, item)) for item in
                ['refGenome', 'genome',
                 'geneCheckBed', 'geneCheckBedDetails',
                 'alignment', 'sequence', 'chromSizes',
                 'outDir'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)
  # existence
  pairs = tuple((item, getattr(args, item)) for item in
                ['geneCheckBed', 'geneCheckBedDetails',
                 'alignment', 'sequence', 'chromSizes',
                 'outDir'])
  for name, value in pairs:
    if not os.path.exists(value):
      parser.error('--%s=%s does not exist' % (name, value))
  # regular file
  pairs = tuple((item, getattr(args, item)) for item in
                ['geneCheckBed', 'geneCheckBedDetails',
                 'alignment', 'sequence', 'chromSizes',
                 ])
  for name, value in pairs:
    if not os.path.isfile(value):
      parser.error('--%s=%s is not a file' % (name, value))
  # directory
  if not os.path.isdir(args.outDir):
    parser.error('--outDir=%s is not a directory' % args.outDir)
