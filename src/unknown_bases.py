import os
from src.abstract_classifier import AbstractClassifier
import lib.psl_genecheck_lib as psl_lib
import lib.sqlite_lib as sql_lib

class UnknownBases(AbstractClassifier):
    @staticmethod
    def __type__():
        return "INT"

    @staticmethod
    def __name__():
        return "UnknownBases"

    def run(self):
        #call methods from abstract_classifier that pull data needed for this classifier
        self.get_transcripts()
        self.make_alignment_dict()

        n_counts = {}
        for transcript in self.transcripts:
            alignments = self.alignment_dict[psl_lib.removeAlignmentNumber(transcript.name),
                    transcript.chromosomeInterval.chromosome]
            n_counts[transcript.name] = sum(a.nCount for a in alignments)

        with sql_lib.ExclusiveSqlConnection(self.output) as cur:
            for transcript, count in n_counts.iteritems():
                sql_lib.upsert(cur, self.genome, self.primary_key, transcript, self.__name__(), count)