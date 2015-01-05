import os
from jobTree.scriptTree.target import Target
from jobTree.src.bioio import logger
import src.psl_genecheck_lib as lib
from collections import defaultdict

class AbstractClassifier(Target):
    def __init__(self, genome, filepaths, originalGeneCheckBed, originalGeneCheckBedDetails, outDir, refGenome):
        #initialize the Target
        Target.__init__(self)

        #store basic information
        self.genome = genome
        self.refGenome = refGenome
        self.orig_bed = originalGeneCheckBed
        self.orig_bed_details = originalGeneCheckBedDetails
        self.output = outDatabase
        self.bed_file = filepaths["bed"]
        self.bed_details_file = filepaths["details"]

        #use dent's library to pull down dicts of alignments, sequences, etc
        self.sequences = lib.getSequences(filepaths["fasta"])
        self.sizes = lib.getChromSizes(filepaths["sizes"])
        self.alignments = lib.getAlignments(filepaths["psl"])

        #really don't need these at all, as orig+psls has all the information
        self.transcripts = lib.getUniqueTranscripts(self.bed_file, self.bed_details_file)
        self.transcript_dict = lib.transcriptListToDict(self.transcripts, noDuplicates=True)

        self.orig_transcripts = lib.getUniqueTranscripts(self.orig_bed, self.orig_bed_details)
        #hope there isn't duplicates here...
        self.orig_transcript_dict = lib.transcriptListToDict(self.orig_transcripts, noDuplicates=True)

        #this dict is keyed on (psl.qName, psl.tName) to ensure uniqueness
        self.alignment_dict = lib.alignmentListToDict(self.alignments, noduplicates=True)

        #the path of the sqlite3 db to be built
        self.db = os.path.join(outDir, refGenome + "_" + genome + ".db")

    def run(self, classifier):
        logger.info("Running classifier {} on {}".format(classifier, self.genome))


    def finish(self, classifier):
        logger.info("Finished classifier {} on {}".format(classifier, self.genome))