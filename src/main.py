import os
from argparse import ArgumentParser
import pysqlite as sql
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions
from src.sqlite_lib import initalizeTable as initalizeTable
from src.unknown_bases import UnknownBases

#classifiers we are currently working with
classifiers = [UnknownBases]

#hard coded file extension types that we are looking for
gene_check_files = {"bed":".coding.gene-check.bed", "details":".coding.gene-check-details.bed"]
sequence_files = {"fasta":".fa", "sizes":".sizes"}
alignment_files = {"psl":".chained.psl"}

def build_parser(args):
    """
    Builds an argument parser for this run
    """
    parser = ArgumentParser()
    parser.add_argument('--refGenome', type=str)
    parser.add_argument('--genomes', nargs="+")
    parser.add_argument('--geneCheckDir', type=DirType)
    parser.add_argument('--alignmentDir', type=DirType)
    parser.add_argument('--sequenceDir', type=DirType)    
    parser.add_argument('--originalGeneCheckBed', type=FileType)
    parser.add_argument('--originalGeneCheckBedDetails', type=FileType)
    parser.add_argument('--outDb', type=dirType, default="results.db")
    return parser


def parse_dir(genomes, targetDir, fileTypes, pathDict={}):
    """
    Given a directory, a dict of filetypes (mapping extension to name), and a list of genomes,
    returns a dict mapping all files matching those fileTypes and checks
    for existence of these files
    """
    for g in genomes:
        for t, ext in fileTypes.iteritems():
            path = os.path.join(targetDir, g + ext)
            if not os.path.exists(path):
                raise RuntimeError("{} does not exist in {}".format(f, targetDir))
            if g not in pathDict:
                pathDict[g] = {}
            pathDict[g][t] = path
    return pathDict


def build_genome_analysis(target, pathDict, args):
    for genome in args.genomes:
        initalizeTable(args.outDb, genome, classifiers)
        paths = pathDict[genome]
        for classifier in classifiers:
            target.addChildTarget(classifier(genome, paths, args.originalGeneCheckBed, 
                    args.originalGeneCheckBedDetails, args.outDb, args.refGenome))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    if not os.path.exists(args.outDb):
        os.mkdir(args.outDb)

    logger.info("Building paths to the required files")
    pathDict = parse_dir(args.genomes, args.geneCheckDir, gene_check_files)
    pathDict = parse_dir(args.genomes, args.alignmentDir, alignment_files, pathDict=pathDict)
    pathDict = parse_dir(args.genomes, args.sequenceDir, sequence_files, pathDict=pathDict)

    refSequence = os.path.join(args.sequenceDir, args.refGenome + ".fa")
    if not os.path.exists(refSequence):
        raise RuntimeError("Reference genome fasta not present at {}".format(refSequence))
    args.refSequence = refSequence

    i = Stack(Target.makeTargetFn(build_analysis, args=(pathDict, args))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")



if __name__ == '__main__':
        main()