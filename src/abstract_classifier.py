from jobTree.scriptTree.target import Target
from jobTree.src.bioio import logger
import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib

class AbstractClassifier(Target):
    def __init__(self, genome, alnPsl, seqFasta, annotationBed, gencodeAttributeMap,  
                geneCheckBed, outDb, refGenome, primaryKey):
        #initialize the Target
        Target.__init__(self)

        #store basic information
        self.genome = genome
        self.refGenome = refGenome
        self.output = outDb
        self.alnPsl = alnPsl
        self.seqFasta = seqFasta
        self.annotationBed = annotationBed
        self.gencodeAttributeMap = gencodeAttributeMap
        self.geneCheckBed = geneCheckBed
        self.primary_key = primaryKey
        self.db = outDb

    def get_original_transcripts(self):
        self.original_transcripts = seq_lib.getTranscripts(self.annotationBed)

    def get_transcript_attributes(self):
        self.attribute_dict = seq_lib.getTranscriptAttributeDict(self.gencodeAttributeMap)

    def get_original_transcript_dict(self):
        if not hasattr(self, 'original_transcripts'):
            self.original_transcripts = seq_lib.getTranscripts(self.annotationBed)
        self.original_transcript_dict = seq_lib.transcriptListToDict(self.original_transcripts, noDuplicates=True)

    def get_transcripts(self):
        self.transcripts = seq_lib.getTranscripts(self.geneCheckBed)

    def get_transcript_dict(self):
        if not hasattr(self, 'transcripts'):
            self.transcripts = seq_lib.getTranscripts(self.geneCheckBed)
        self.transcript_dict = seq_lib.transcriptListToDict(self.transcripts, noDuplicates=True)

    def get_seq_dict(self):
        self.seq_dict = seq_lib.getFastaDict(self.seqFasta, upper=True)

    def get_alignments(self):
        self.alignments = psl_lib.readPsl(self.alnPsl, uniqify=True)

    def make_alignment_dict(self):
        if not hasattr(self, 'alignments'):
            self.get_alignments()
        self.alignment_dict = psl_lib.getPslDict(self.alignments, noDuplicates=True)

    def upsert_wrapper(self, cur, alignmentName, value):
        """convenience wrapper for upserting into a column in the sql lib.
        So you don't have to call __name__, self.primaryKey, etc each time"""
        sql_lib.upsert(cur, self.genome, self.primary_key, alignmentName, 
                self.__class__.__name__, value)