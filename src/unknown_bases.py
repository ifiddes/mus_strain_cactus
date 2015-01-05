import os
from src.abstract_classifier import AbstractClassifier
import src.psl_genecheck_lib as lib
import src.sqlite_lib as sql_lib

class UnknownBases(AbstractClassifier):
    def run(self):
        AbstractClassifier.run(self, "UnknownBases")
        n_counts = {}
        for transcript in self.transcripts:
            orig_transcript = self.orig_transcripts(lib.removeAlignmentNumber(transcript.name))
            alignments = self.alignments_dict(lib.removeAlignmentNumber(transcript.name),
                    transcript.chromosomeInterval.chromosome)
             n_counts[transcript.name] = sum(a.nCount for a in alignments)

        con = sql_lib.ExclusiveSqlConnection(self.output)
        if not con.hasTable(self.genome):
            