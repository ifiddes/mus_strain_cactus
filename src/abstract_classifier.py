from jobTree.scriptTree.target import Target
from jobTree.src.bioio import logger
import lib.psl_genecheck_lib as psl_lib

class AbstractClassifier(Target):
    def __init__(self, genome, filepaths, originalGeneCheckBed, originalGeneCheckBedDetails, outDb, refGenome, primaryKey):
        #initialize the Target
        Target.__init__(self)

        #store basic information
        self.genome = genome
        self.refGenome = refGenome
        self.orig_bed = originalGeneCheckBed
        self.orig_bed_details = originalGeneCheckBedDetails
        self.output = outDb
        self.bed_file = filepaths["bed"]
        self.bed_details_file = filepaths["details"]
        self.primary_key = primaryKey
        self.db = outDb
        self.filepaths = filepaths

    def get_sequences(self):
        self.sequences = psl_lib.getSequences(self.filepaths["fasta"])

    def get_sizes(self):
        self.sizes = psl_lib.getChromSizes(self.filepaths["sizes"])

    def get_alignments(self):
        self.alignments = psl_lib.getAlignments(self.filepaths["psl"])

    def make_alignment_dict(self):
        if not hasattr(self, 'alignments'):
            self.get_alignments()
        self.alignment_dict = psl_lib.alignmentListToDict(self.alignments)

    def get_transcripts(self):
        self.transcripts = psl_lib.getUniqueTranscripts(self.bed_file, self.bed_details_file)

    def make_transcript_dict(self, noDuplicates=True):
        if not hasattr(self, 'transcripts'):
            self.get_transcripts()
        self.transcript_dict = psl_lib.getUniqueTranscripts(self.transcripts, noDuplicates)

    def get_orig_transcripts(self):
        self.orig_transcripts = psl_lib.getUniqueTranscripts(self.orig_bed, self.orig_bed_details)
        
    def make_orig_transcript_dict(self, noDuplicates=True):
        if not hasattr(self, 'orig_transcripts'):
            self.get_orig_transcripts()
        self.orig_transcript_dict = psl_lib.transcriptListToDict(self.orig_transcripts, noDuplicates)