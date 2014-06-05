
"""The following data types are used for iterating over gene-check-detail and gene-check bed files.
An example of entries from such files.

[benedict@hgwdev tracks]$ more Rattus.coding.gene-check-details.bed
1       2812370 2812372 noStop/ENSMUST00000065527.4
1       2812370 2812372 noStop/ENSMUST00000095795.4

[benedict@hgwdev tracks]$ more Rattus.coding.gene-check.bed 
1       2812346 3113743 ENSMUST00000178026.1    0       -       2812370 3038729
128,0,0 9       54,2,89,249,90,165,105,13,45    0,58,62,698,1209,1305,226292,301
050,301352
1       2812346 3113783 ENSMUST00000095795.4    0       -       2812370 3038729
128,0,0 9       54,2,89,249,197,52,105,13,85    0,58,62,698,1209,1418,226292,301
050,301352
"""

class ChromosomeInterval:
    """Represents an interval of a chromosome. BED coordinates, strand is True, False or NULL (if no strand)
    """
    def __init__(self, chromosome, start, stop, strand):
        self.chromosome = str(chromosome)
        self.start = int(start)
        self.stop = int(stop)
        self.strand = strand
    
class TranscriptAnnotation: #Maps to bed 4 + annotation field
    """Represents an annotation of a transcript, from one of the classification bed files
    """
    def __init__(self, chromosomeInterval, transcript, annotation):
        self.chromosomeInterval = chromosomeInterval
        self.transcript = str(transcript)
        self.annotation = annotation

class Transcript: 
    """Represent a transcript.
    """
    def __init__(self, chromosomeInterval, transcript, exons, annotations):
        self.chromosomeInterval = chromosomeInterval
        self.transcript = str(transcript)
        self.exons = exons #Is a list of chromosome intervals
        self.annotations = annotations #Is a list of transcript annotations

def tokenizeBedFile(bedFile):
    """Iterator through bed file, returning lines as list of tokens
    """
    fileHandle = open(bedFile, 'r')
    for line in fileHandle:
        if line != '':
            tokens = line.split()
            yield tokens
    fileHandle.close()

def transcriptIterator(transcriptsBedFile, transcriptClassificationBedFile):
    """Iterates over the transcripts detailed in the two files, producing Transcript objects.
    """
    transcriptsAnnotations = {}
    for bedTokens in tokenizeBedFile(transcriptClassificationBedFile):
        assert len(bedTokens) == 5 #We expect to be able to get 5 fields out of this bed.
        tA = TranscriptAnnotation(ChromosomeInterval(tokens[0], tokens[1], tokens[2], None), tokens[3], tokens[4])
        if tA.transcript not in transcriptsClassifications:
            transcriptsAnnotations[tA.transcript] = []
        transcriptsAnnotations[tA.transcript].append(tA)
    
    for bedTokens in tokenizeBedFile(transcriptsBedFile):
        assert len(bedTokens) == 12 #We expect to be able to get 12 fields out of this bed.
        #Transcript
        transcript = tokens[3]
        #Get the chromosome interval
        assert tokens[5] in ('+', '-')
        cI = ChromosomeInterval(tokens[0], tokens[1], tokens[2], tokens[5] == '+')
        #Get the exons
        def getExons(exonNumber, blockSizes, blockStarts):
            assert exonNumber == len(blockSizes)
            assert exonNumber == len(blockStarts)
            return [ ChromosomeInterval(cI.chromosome, cI.start + int(blockStarts[i]), cI.start + int(blockStarts[i]) + int(blockSizes[i]), cI.strand) \
                    for i in range(exonNumber) ]
        exons = getExons(int(tokens[9]), ",".split(tokens[10]), ",".split(tokens[11]))
        #Get the transcript annotations
        annotations = []
        if transcript in transcriptsAnnotations:
            annotations = transcriptsAnnotations[transcript]
        yield Transcript(chromosomeInterval, transcript, exons, annotations)
        
    