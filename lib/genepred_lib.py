"""
UNTESTED genepred lib if you don't want to use BED files
original author: Dent Earl
modified by: Ian Fiddes
"""

from sequence_lib import ChromosomeInterval, Sequence


class genePredRow(object):
    """
    Represents one transcript from a genePred (.gp) file. Coordinates are translated
    into ChromosomeInterval objects representing the mRNA, CDS and each Exon
    """
    __slots__ = ("name", "transcript", "cds", "exons", "thickStart", "thickEnd")
    
    def __init__(self, line):
        data = line.split()
        assert len(data) == 14
        
        self.name = data[0]
        
        chrom, strand = data[1:3]
        txStart, txEnd, cdsStart, cdsEnd, exonCount = map(int, data[3:8])
        exonStarts = map(int, data[8].split(",")[:-1]) #extra comma on end
        exonEnds = map(int, data[9].split(",")[:-1])

        self.transcript = ChromosomeInterval(chrom, strand, txStart, txEnd)
        self.cds = ChromsomeInterval(chrom, strand, cdsStart, cdsEnd)

        #for legacy reasons (not messing with Dent's code), we store
        #cdsStart also as thickStart and cdsEnd also as thickEnd
        self.thickStart = cdsStart
        self.thickEnd = cdsEnd

        self.exons = []
        for start, stop in zip(exonStarts, exonEnds):
            self.exons.append(ChromosomeInterval(chrom, strand, start, stop))

    def __eq__(self, other):
        return self.hashkey() == other.hashkey()

    def hashkey(self):
        """
        Returns a unique* hashable key for this row.

        *In testing this appears to almost always be true. In (rare) cases where it is not true, 
        the genepred lines are identical."""
        return (self.name, self.chrom, self.txStart, self.txEnd, self.cdsStart, self.cdsEnd)


    def hashkey(self, deuniquify=False):
        """ return a string to use as dict key.
        """
        if deuniquify:
            return '%s_%s_%d_%d' % (removeAlignmentNumber(self.name),
                                        self.chromosomeInterval.chromosome,
                                        self.chromosomeInterval.start,
                                        self.chromosomeInterval.stop)
        else:
            return '%s_%s_%d_%d' % (self.name, self.chromosomeInterval.chromosome,
                                        self.chromosomeInterval.start,
                                        self.chromosomeInterval.stop)

    def mRnaCoordinateToCodon(self, p):
        """ Take position P with 0-based mRNA-relative position and convert it
        to 0-based (codon, codon position) tuple.
        """
        # (strand irrelevant)
        #       0    5    10
        #       |    |    |
        # mRNA  +++++++++++
        # codon +++++++++++
        #       |  |  |  |
        #       0,0   2,0
        if p is None: return None
        assert(p >= 0)
        assert(p < sum([(e.stop - e.start) for e in self.exons]))    # could be a tighter bound
        return (int(p / 3), p % 3)

    def codonCoordinateToMRna(self, p):
        """ Take (codon, codon position) tuple P and convert it
        to 0-based mRNA-relative position.
        """
        # (strand irrelevant)
        #       0    5    10
        #       |    |    |
        # mRNA  +++++++++++
        # codon +++++++++++
        #       |  |  |  |
        #       0,0   2,0
        if p is None: return None
        assert(p[0] >= 0)
        assert(p[1] >= 0)
        return p[0] * 3 + p[1]

    def codonCoordinateToChromosome(self, p):
        """ Take (codon, codon position) tuple P and convert it
        to 0-based chromosome-relative position.
        """
        m = self.codonCoordinateToMRna(p)
        return self.mRnaCoordinateToChromosome(m)

    def chromosomeCoordinateToCodon(self, p):
        """ Take 0-based chromosome-relative P and convert it
        to (codon, codon position) tuple position.
        """
        m = self.chromosomeCoordinateToMRna(p)
        return self.mRnaCoordinateToCodon(m)

    def mRnaCoordinateToExon(self, p):
        """ Take position P with 0-based mRNA-relative position and convert it
        to 0-based exon-relative position.
        """
        assert(len(self.exons))
        if p is None: return None
        if p < 0: return None
        if p >= sum([(e.stop - e.start) for e in self.exons]): return None
        # positive strand
        # chromosome +++++++++++
        #            |    |    |
        # exon        ..++ ++++.  two exons (thick and thin parts)
        #             |     |  |
        # mrna          ++ ++++
        #               |     |
        # so to go from mrna to exon, we must add on the difference
        # between the thick start and thin start from the "start".
        if self.chromosomeInterval.strand:
            # positive strand, offset is first exon start to thickStart
            for e in self.exons:
                if e.start < self.thickStart and e.stop <= self.thickStart:
                    # add the whole exon to the offset
                    p += e.stop - e.start
                elif e.start < self.thickStart and self.thickStart <= e.stop:
                    # only add the thin part of this exon
                    p += self.thickStart - e.start
                    break
        else:
            for e in reversed(self.exons):
                if self.thickEnd < e.start and self.thickEnd < e.stop:
                    # add the whole exon to the offset
                    p += e.stop - e.start
                elif e.start < self.thickEnd and self.thickEnd < e.stop:
                    # only add the thin part of this exon
                    p += e.stop -    self.thickEnd
                    break
        return p

    def exonCoordinateToMRna(self, p):
        """ Take position P with 0-based exon-relative position and convert it
        to 0-based mRNA-relative position. If position does not exist in mRNA,
        return None.
        """
        if p is None: return None
        # find the thickStart, thickEnd offsets in exon coordinates
        exonThickStart, exonThickEnd = None, None
        x = 0    # exon coordinate
        for e in self.exons:
            length = e.stop - e.start
            if exonThickStart is None and e.start >= self.thickStart:
                # thickStart fell between exons
                exonThickStart = x
            if exonThickStart is None and e.stop > self.thickStart:
                # exon contains thickStart
                exonThickStart = x + self.thickStart - e.start
            if exonThickEnd is None and e.start >= self.thickEnd:
                # thickEnd fell between exons
                exonThickEnd = x
            if exonThickEnd is None and e.stop >= self.thickEnd:
                # exon contains thickEnd
                exonThickEnd = x + self.thickEnd - e.start
            x += length
        if not self.chromosomeInterval.strand:
            exonThickStart, exonThickEnd = exonThickEnd, exonThickStart
            exonThickStart = x - exonThickStart
            exonThickEnd = x - exonThickEnd
        if p < exonThickStart:
            return None
        if p >= exonThickEnd:
            return None
        return p - exonThickStart

    def mRnaCoordinateToChromosome(self, p):
        """ Take position P with 0-based mRNA-relative position and convert it
        to 0-based chromosome-relative position.
        """
        assert(len(self.exons))
        if p is None: return None
        if p < 0: return None
        limit = sum([(e.stop - e.start) for e in self.exons])
        if p >= limit: return None
        p = self.mRnaCoordinateToExon(p)
        if p >= limit: return None
        return self.exonCoordinateToChromosome(p)

    def exonCoordinateToChromosome(self, p):
        """ Take position P with 0-based exon-relative position and convert it
        to 0-based chromosome-relative position.
        """
        if p is None: return None
        if p < 0:
            return None
        if p >= sum([(e.stop - e.start) for e in self.exons]):
            return None
        assert(len(self.exons))
        c = 0    # cumulative position through exon space
        if not self.chromosomeInterval.strand:
            p = sum([(e.stop - e.start) for e in self.exons]) - 1 - p
        e_start = self.exons[0].start
        for e in self.exons:
            if p < c + e.stop - e.start:
                # the position is within this exon
                return p - c + e.start
            else:
                # sorry mario, your position is in another exon
                c += e.stop - e.start
        assert(False)    # we should never get here

    def chromosomeCoordinateToExon(self, p):
        """ Take position P with 0-based chromosome-relative position and convert it
        to 0-based exon-relative position. If position does not exist in
        exon, return None.
        """
        if p is None: return None
        if self.chromosomeInterval.strand:
            def _stranded(v): return v
        else:
            def _stranded(v):
                return sum([(e.stop - e.start) for e in self.exons]) - 1 - v
        c = 0    # cumulative position through exon space
        e_start = self.exons[0].start
        for e in self.exons:
            if p < e.start:
                # p is not in an exon
                return None
            if p < e.stop:
                # the position is within this exon
                return _stranded(c + p - e.start)
            else:
                # sorry mario, your position is in another exon
                c += e.stop - e.start
        return None

    def chromosomeCoordinateToMRna(self, p):
        """ Take position P with 0-based chromosome-relative position and convert it
        to 0-based mRNA-relative position. If position does not exist in mRNA,
        return None.
        """
        if p is None: return None
        if p < 0:
            return None
        if p >= self.chromosomeInterval.stop:
            return None
        q = self.chromosomeCoordinateToExon(p)
        if q is None:
            return None
        if q < 0:
            return None
        if q >= sum([(e.stop - e.start) for e in self.exons]):
            return None
        return self.exonCoordinateToMRna(q)

    def getMRna(self, sequence):
        """ Return the mRNA sequence for the transcript (based on the exons) using
        a SEQUENCE object as the source for dna sequence.
        The returned sequence is in the correct 5'-3' orientation (i.e. it has
        been reverse complemented if necessary).
        """
        assert(self.chromosomeInterval.chromosome == sequence.name)
        assert(self.chromosomeInterval.stop <= sequence.getLength())
        s = ''
        # chromosome    ttttTTTTTTTTTTTtttt  t: thin T: THICK
        # exon            eeeeee eeee eee
        # mrna              mmmm mmmm m
        for e in self.exons:
            if self.thickStart < e.start and e.stop < self.thickEnd:
                # squarely in the CDS
                s += sequence.sliceSequence(e.start, e.stop)
            elif (e.start <= self.thickStart and e.stop < self.thickEnd
                        and self.thickStart < e.stop):
                # thickStart marks the start of the mRNA
                s += sequence.sliceSequence(self.thickStart, e.stop)
            elif e.start <= self.thickStart and self.thickEnd <= e.stop:
                # thickStart and thickEnd mark the whole mRNA
                s += sequence.sliceSequence(self.thickStart, self.thickEnd)
            elif (self.thickStart < e.start and self.thickEnd <= e.stop
                        and e.start < self.thickEnd):
                # thickEnd marks the end of the mRNA
                s += sequence.sliceSequence(e.start, self.thickEnd)
        if not self.chromosomeInterval.strand:
            s = reverseComplement(s)
        return s

    def getMRnaCodon(i, sequence):
        """ Given a nonnegative number I, and a SEQUENCE object,
        return the three nucleotide codon from the mRNA.
        """
        return self.getMRna(sequence)[i * 3:i * 3 + 3]

    def bedString(self):
        """ Write a transcript object to the given file.
        """
        strandChar = '-'
        if self.chromosomeInterval.strand:
                strandChar = '+'
        return '\t'.join(
            [self.chromosomeInterval.chromosome,
             str(self.chromosomeInterval.start),
             str(self.chromosomeInterval.stop),
             self.name, str(self.score), strandChar,
             str(self.thickStart), str(self.thickEnd),
             self.itemRgb, str(len(self.exons)),
             ','.join([str(exon.stop - exon.start) for exon in self.exons]),
             ','.join([str(exon.start - self.chromosomeInterval.start)
                                 for exon in self.exons])])

    def getIntrons(self):
        """Get a list of ChromosomeIntervals representing the introns for this
        transcript. The introns are in *+ strand of CHROMOSOME* ordering,
        not the order that they appear in the transcript!"""
        introns = []
        prevExon = None
        for exon in self.exons:
            if prevExon is not None:
                assert exon.start > prevExon.stop
                assert exon.strand == prevExon.strand
                intron = ChromosomeInterval(exon.chromosome,
                                            prevExon.stop,
                                            exon.start,
                                            exon.strand)
                introns.append(intron)
            prevExon = exon
        return introns

    def __cmp__(self, transcript):
        return cmp((self.chromosomeInterval, self.name),
                    (transcript.chromosomeInterval, transcript.name))


_nuc_pairs = [('a', 't'), ('g', 'c'), ('n', 'n')]
_complement = {}
for a, b in _nuc_pairs:
    _complement[a], _complement[b] = b, a
    _complement[a.upper()], _complement[b.upper()] = b.upper(), a.upper()
_complement['-'] = '-'


def complement(seq):
    """ given a sequence, return the complement.
    """
    seq = ''.join([_complement[s] for s in seq])
    return seq


def reverseComplement(seq):
    """ Given a sequence, return the reverse complement.
    """
    seq = seq[::-1]    # reverse
    seq = complement(seq)
    return seq

