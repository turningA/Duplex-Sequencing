#by Mike Schmitt and Scott Kennedy
#version 1.1
#June 14, 2012

import sys
import pysam
from optparse import OptionParser

"""

    This program takes as input a BAM file with Duplex Tags in the header, and searches for the partner tag. Reads with paired tags are kept. Non-agreeing positions within a read are replaced with N's.
    
    Usage: python DCSfilter.py --infile seqfile.bam --outfile seqfile.bam.DCSfilter
    
"""
    
parser=OptionParser()
parser.add_option("--infile", action="store", type='string', dest="infile", help="input BAM file", default='sys.stdin')
parser.add_option("--outfile",  action="store", type='string', dest="outfile", help="output BAM file",  default='/dev/stdout')
parser.add_option("--readnumloc",  action="store", type='int', dest="readnumloc", help="header field containing read number",  default='2')
parser.add_option("--tagloc",  action="store", type='int', dest="tagloc", help="header field containing duplex tag",  default='1')
parser.add_option("--index", action="store", type='string', dest="indexSeq", help="Sequence of sample index.  Needed for GATK support")
o, args = parser.parse_args()

inBam = pysam.Samfile( o.infile, "rb" )
readDict = {}

dictctr = 0
seqctr = 0
tagmatchctr = 0
partialmatchctr = 0
seqreplacectr = 0
Nctr = 0

#first, build a dictionary with read 1 duplex tags as key, and the corresponding sequence as an entry.
for line in inBam :

    lineSplit = line.qname
    read = lineSplit.split(":")[o.readnumloc]
    tag = lineSplit.split(":")[o.tagloc]
    if read == '1' and tag not in readDict :
 
        readDict[tag] = [line.seq, '']
        dictctr += 1

    if dictctr % 1000000 == 0 :
        print >> sys.stderr, "sequences added to dictionary:", dictctr
        dictctr += 1

inBam.close()
inBam = pysam.Samfile( o.infile, "rb" )

#next, evaluate every read 2 tag for a match in the dictionary

for line in inBam :

    seqctr += 1
    lineSplit = line.qname
    read = lineSplit.split(":")[o.readnumloc]     
    tag = lineSplit.split(":")[o.tagloc]
    switchtag = tag[12:24] + tag[:12]
    
    if read == '2' and switchtag in readDict :
        tagmatchctr += 1

        if len(line.seq) == len(readDict[switchtag][0]) :
            newSeq = ''
            for i in xrange (len(line.seq) ) :
                
                if line.seq[i] == readDict[switchtag][0][i] :
                    newSeq = newSeq + line.seq[i]
                else :
                    newSeq = newSeq + 'N'
                   
        
            if line.seq != readDict[switchtag][0] and newSeq.count('N') < 20 :
                partialmatchctr += 1
                Nctr += newSeq.count('N')
                    
            if newSeq.count('N') < (readDict[switchtag][1]).count('N') or ( readDict[switchtag][1] == '' and newSeq.count('N') < 20 ) :
                readDict[switchtag][1] = newSeq
                seqreplacectr += 1

    if seqctr % 1000000 == 0 :
        print >> sys.stderr, "tags processed for matches:", seqctr
        print >> sys.stderr, "tag matches:", tagmatchctr
        print >> sys.stderr, "total sequence matches:", seqreplacectr
        print >> sys.stderr, "reads containing disagreeing bases (replaced with N's):", partialmatchctr
        print >> sys.stderr, "number of N's added:", Nctr


inBam.close()

# Done generating tag dictionary.  Reinterate over bamfile and write entries that have a sequence match.

inBam = pysam.Samfile( o.infile, "rb" )
outBam = pysam.Samfile ( o.outfile, "wb", template=inBam)

printlinectr = 0
printlinematch = 0

for line in inBam :
    
    printlinectr += 1
    
    lineSplit = line.qname

    tag = lineSplit.split(":")[o.tagloc]
    read = lineSplit.split(":")[o.readnumloc]     
            
    if tag in readDict and read == '1' and len (readDict[tag][1]) > 0 :
        
        line.seq = readDict[tag][1]
        line.qual = 'i' * line.rlen
        line.qname = line.qname + '#' + o.indexSeq
        readDict[tag][1] = ''
        printlinematch += 1
        outBam.write(line)

    if printlinectr % 1000000 == 0:
        print >> sys.stderr, "Lines evaluated for printing:", printlinectr
        print >> sys.stderr, "Matching sequences printed:", printlinematch


print >> sys.stderr, "Total tags processed for matches:", seqctr
print >> sys.stderr, "Total tag matches:", tagmatchctr
print >> sys.stderr, "Total sequence matches:", seqreplacectr
print >> sys.stderr, "Total reads containing disagreeing bases (replaced with N's):", partialmatchctr
print >> sys.stderr, "total number of N's added:", Nctr

inBam.close()
outBam.close()
