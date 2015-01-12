import os
import argparse
import sqlite3 as sql
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions
from lib.sqlite_lib import initializeTable, insertRow
from lib.general_lib import FileType, DirType, FullPaths

#basic_attributes has many basic classes stuck together that are not really classifiers
from src.basic_attributes import *
#psl_attributes takes attributes from the psl and makes columns, mainly genomic positions
from src.psl_attributes import *

from src.unknown_bases import UnknownBases
from src.end_stop import EndStop
from src.begin_start import BeginStart
from src.in_frame_stop import InFrameStop

#classifiers we are currently working with
classifiers = [EndStop, UnknownBases, BeginStart, InFrameStop]

#add in all of the basic attribute columns
#classifiers = classifiers + [TranscriptID, GeneID, GeneName, GeneType, TranscriptType]
#add in all of the psl attribute columns
#classifiers = classifiers + [SourceChrom, SourceStart, SourceStop, SourceStrand,
#                            DestChrom, DestStart, DestStop, DestStrand]


#hard coded file extension types that we are looking for
alignment_ext = ".filtered.psl"
sequence_ext = ".2bit"
gene_check_ext = ".bed"
#gene_check_details_ext = ".coding-gene-check-details.bed"

def build_parser():
    """
    Builds an argument parser for this run
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--refGenome', type=str)
    parser.add_argument('--genomes', nargs="+")
    parser.add_argument('--annotationBed', type=FileType)
    parser.add_argument('--dataDir', type=DirType, action=FullPaths)
    parser.add_argument('--gencodeAttributeMap', type=FileType)
    parser.add_argument('--outDb', type=str, default="results.db")
    parser.add_argument('--primaryKey', type=str, default="AlignmentID")
    parser.add_argument('--overwriteDb', action="store_true")
    return parser


def parse_dir(genomes, targetDir, ext):
    pathDict = {}
    for g in genomes:
        path = os.path.join(targetDir, g + ext)
        if not os.path.exists(path):
            raise RuntimeError("{} does not exist".format(path))
        pathDict[g] = path
    return pathDict


def build_analysis(target, alnPslDict, seqFastaDict, geneCheckBedDict, gencodeAttributeMap,
            genomes, annotationBed, outDb, primaryKeyColumn, refGenome):
    for genome in genomes:
        alnPsl, seqFasta = alnPslDict[genome], seqFastaDict[genome]
        geneCheckBed = geneCheckBedDict[genome]
        initialize_sql_columns(genome, outDb, primaryKeyColumn)
        for classifier in classifiers:
            target.addChildTarget(classifier(genome, alnPsl, seqFasta, annotationBed,
                    gencodeAttributeMap, geneCheckBed, outDb, refGenome, primaryKeyColumn))


def initialize_sql_columns(genome, outDb, primaryKeyColumn):
    con = sql.connect(outDb)
    columns = [[x.__name__, x.__type__()] for x in classifiers]
    with con:
        initializeTable(con.cursor(), genome, columns, primaryKeyColumn)


def initialize_sql_rows(genome, outDb, alnPsl, primaryKeyColumn):
    con = sql.connect(outDb)
    alnIds = set(x.split()[9] for x in open(alnPsl))
    for alnId in alnIds:
        insertRow(con.cursor(), genome, primaryKeyColumn, alnId)


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    if args.overwriteDb is True and os.path.exists(args.outDb):
        os.remove(args.outDb)

    logger.info("Building paths to the required files")
    alnPslDict = parse_dir(args.genomes, args.dataDir, alignment_ext)
    seqFastaDict = parse_dir(args.genomes, args.dataDir, sequence_ext)
    geneCheckBedDict = parse_dir(args.genomes, args.dataDir, gene_check_ext)
    #geneCheckBedDetailsDict = parse_dir(args.genomes, args.geneCheckDir, gene_check_details_ext)

    refSequence = os.path.join(args.dataDir, args.refGenome + ".fa")
    if not os.path.exists(refSequence):
        raise RuntimeError("Reference genome fasta not present at {}".format(refSequence))
    args.refSequence = refSequence

    i = Stack(Target.makeTargetFn(build_analysis, args=(alnPslDict, seqFastaDict, geneCheckBedDict, 
            args.gencodeAttributeMap, args.genomes, args.annotationBed, args.outDb, args.primaryKey, 
            args.refGenome))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
        main()