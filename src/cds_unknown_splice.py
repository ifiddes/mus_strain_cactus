from collections import defaultdict

from src.abstract_classifier import AbstractClassifier
import lib.sequence_lib as seq_lib


class CdsUnknownSplice(AbstractClassifier):
    """

    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.

    Since sqlite3 lacks a BOOL type, reports 1 if TRUE and 0 if FALSE

    """
    @staticmethod
    def __type__():
        return "INTEGER"

    def bad_splice(self, donor, acceptor):
        m = {"GT":"AG", "GC":"AG", "AT":"AC"}
        d = donor.upper()
        a = acceptor.upper()
        if d in m and m[d] != a:
            return True
        else:
            return False

    def run(self, minimum_intron_size=30):
        self.get_transcript_dict()
        self.get_seq_dict()

        s_dict = defaultdict(int)
        for a, t in self.transcript_dict.iteritems():
            chrom_seq = self.seq_dict[t.chromosomeInterval.chromosome]
            for i in xrange(len(t.intronIntervals)):
                if t.exons[i].containsCds() is True and t.exons[i+1].containsCds() is True:
                    if len(t.intronIntervals[i]) >= minimum_intron_size:
                        donor = chrom_seq[t.intronIntervals[i].start : t.intronIntervals[i].start + 2]
                        acceptor = chrom_seq[t.intronIntervals[i].stop - 2 : t.intronIntervals[i].start]
                        if self.bad_splice(donor, acceptor) is True:
                            s_dict[a] = 1
                            break

        self.upsert_dict_wrapper(s_dict)