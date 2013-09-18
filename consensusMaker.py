#By Scott Kennedy
#Version 1.1
#June 14, 2012

import sys
import pysam
import re
from collections import defaultdict
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", default='sys.stdin')
parser.add_argument("--tagfile",  action="store",  dest="tagfile", help="output tagcounts file",  default='sys.stdout')
parser.add_argument("--fastqfile",  action="store", dest="fastqfile", help="output FASTQ file",  default='sys.stdout')
parser.add_argument("--rep_filt", action="store",  type=int, dest='rep_filt', help="Remove tags with homomeric runs of nucleotides of length x", default=9 )
parser.add_argument('--minmem', type=int, default=0, dest='minmem', help="Minimum number of reads allowed to comprise a consensus.")
parser.add_argument('--maxmem', type=int, default=100, dest='maxmem', help="Maximum number of reads allowed to comprise a consensus.")
parser.add_argument('--cutoff', type=float, default=0, dest='cutoff', help="Percentage of nucleotides at a given position in a read that must be identical in order for a consensus to be called at that position.")
parser.add_argument('--readlength', type=int, default=81, dest='read_length', help="Length of the input read that is being used.")
parser.add_argument('--readnum', default=None, dest='readnum', help="Read identifier (1 or 2)")
parser.add_argument('-p', action='store_true', dest='pipe', help="Output consensus reads to stdout"  )
parser.add_argument('--ref_filt',  action='store', dest='ref_filt', default=[None],  nargs='*')
o = parser.parse_args()

def consensusMaker (groupedReadsList,  cutoff,  readLength) :
    '''The consensus maker uses a simple "majority rules" algorithm to make a consensus at each base position.  If no nucleotide majority reaches above the minimum theshold (--cutoff), the position is
considered undefined and an 'N' is placed at that position in the read.'''
    nucIdentityList=[0, 0, 0, 0, 0, 0] # In the order of T, C, G, A, N, Total
    nucKeyDict = {0:'T', 1:'C', 2:'G', 3:'A', 4:'N'}
    seqDict = {}
    consensusRead = ''
    
    for i in xrange(readLength) : #Count the types of nucleotides at a position in a read. i is the nucleotide index within a read in groupedReadsList
        for j in xrange(len(groupedReadsList)) : #Do this for every read that comprises a SMI group. j is the read index within groupedReadsList
            
            if groupedReadsList[j][i] == 'T' :
                nucIdentityList[0] += 1
                nucIdentityList[5] += 1
            elif groupedReadsList[j][i] == 'C':
                nucIdentityList[1] += 1
                nucIdentityList[5] += 1
            elif groupedReadsList[j][i] == 'G':
                nucIdentityList[2] += 1
                nucIdentityList[5] += 1 
            elif groupedReadsList[j][i] == 'A':
                nucIdentityList[3] += 1
                nucIdentityList[5] += 1
            elif groupedReadsList[j][i] == 'N':
                nucIdentityList[4] += 1
                nucIdentityList[5] += 1
            
            seqDict[i] = nucIdentityList
        nucIdentityList=[0, 0, 0, 0, 0, 0] #reset for the next nucleotide position
    
    for i in xrange(readLength) :#rebuild consensus read taking into account the cutoff percentage
        
        for j in [0, 1, 2, 3, 4] :
            
            if float(seqDict[i][j])/float(seqDict[i][5]) > cutoff :
                consensusRead += nucKeyDict[j]
                break
            elif float(seqDict[i][j])/float(seqDict[i][5]) < cutoff :
                if j==4:
                    consensusRead += 'N'
                else :
                    pass
            else :
                pass
    return consensusRead

inBam = pysam.Samfile( o.infile, "rb" )
tagFile = open(o.tagfile,  'w')
groupedReadsList = []
barcodedict = {}
readDict = {}
qualScore = 'i'*o.read_length
tagDict=defaultdict( lambda: 0 )
groupdict = defaultdict( lambda: 0 )

if o.ref_filt == [None] :
    ref_filt = [None] #['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', 'X', 'Y']
else :
    ref_filt = o.ref_filt

for chr in ref_filt :
    
    bamEntry = inBam.fetch(until_eof=True)
    
    for line in bamEntry:
        
        #header = line.qname.split('\t')[0]
        tag = line.qname.split('#')[1].split('/')[0]
        
        if ('A'*o.rep_filt in tag) or ('C'*o.rep_filt in tag) or ('G'*o.rep_filt in tag) or ('C'*o.rep_filt in tag) :
            pass 
        else :
            tagDict[ tag ] += 1
        
tagFile.write ( "\n".join( [ "%s\t%d" % ( SMI, tagDict[SMI] ) for SMI in sorted( tagDict.keys(), key=lambda x: tagDict[x], reverse=True ) ] ))
tagDict={} #Clear from memory for now.  In the future I will attempt to use this data structure to do the group naming entirely in memory .

inBam.close()
tagFile.close()
tagFile = open(o.tagfile, 'r')
fastqFile = open(o.fastqfile, 'w')

ctr = 1

for line in tagFile :
    
    line = line.rstrip( '\n' )
    linesplit = line.split ( '\t' )
    barcodedict[ linesplit[ 0 ] ] = ( ctr, int( linesplit[ 1 ] ) )
    ctr += 1
    
tagFile.close()
inBam = pysam.Samfile( o.infile, "rb" )
ctr = 0

for chr in ref_filt :
    
    bamEntry = inBam.fetch(until_eof=True)
    
    for line in bamEntry:
    
        barcode = line.qname.split("#")[1].split('/')[0]
        
        if barcode not in barcodedict : #skip barcodes in bam file that are not found in the tagcounts file
            continue
        
        groupnum, groupmem = barcodedict[ barcode ]
        groupdict[ groupnum ] += 1
        readnum = groupdict[ groupnum ]
    
        if barcode not in readDict and int(readnum) <= o.maxmem and int(groupmem) >= o.minmem : #if a barcode that passes min and max criteria is not already in the read dictionary, add it, along with relevant information, to the read dictionary 
            readDict[barcode] = [groupnum,  readnum,  groupmem,  line.seq]
    
        elif barcode in readDict and int(readnum) <= o.maxmem and int(groupmem) >= o.minmem : #if read has already been added, increment the number of times it's been seen and append the sequence to the list for that barcode
            readDict[barcode][1] = readnum
            readDict[barcode].append(line.seq)
    
            if readDict[barcode][1] == int(groupmem) : #Submit read group to consensusMaker() when all the reads in that group have been encountered.
            
                consensus=consensusMaker(readDict.pop(barcode)[3:],  o.cutoff,  o.read_length)
                
                if o.pipe == True : #print to stdout if the -p option is invoked
                    print('@:%s:%s:%s_%s\n%s\n+\n%s' %(barcode, o.readnum, groupnum, groupmem, consensus, qualScore))
                
                fastqFile.write('@:%s:%s:%s_%s\n%s\n+\n%s\n' %(barcode, o.readnum, groupnum, groupmem, consensus, qualScore))
            else :
                pass
            
inBam.close()
fastqFile.close()
