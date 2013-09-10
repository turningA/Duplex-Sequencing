import sys
import pysam
import re
from collections import defaultdict
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
        "--inRead1", 
        action = "store", 
        dest = "inRead1", 
        help = "Input of read 1 of hairpin sequence"
        )
parser.add_argument(
        "--inRead2", 
        action = "store", 
        dest = "inRead2", 
        help = "Input of read 2 of hairpin sequence"
        )
parser.add_argument(
        "--out_prefix", 
        type = str, 
        default = None, 
        help = "Optional; prefix for the output files"
        )
parser.add_argument(
        "--spacers", 
        type = str, 
        dest = 'spacers', 
        help = "Potential spacers in hairpin, in format \
            'spacer1, spacer2, spacer3, ...'"
            )
parser.add_argument(
        "--slength", 
        type = int, 
        dest = 'slength', 
        help = "Length of the spacer sequence"
        )
parser.add_argument(
        "--blength", 
        type = int, 
        dest = 'blength', 
        help = "Length of the barcode sequence"
        )
parser.add_argument(
        "--plengths", 
        type = str, 
        dest = 'plengths', 
        help = "Length of each primer, coresponds directly to length trimmed \
            off.  Takes the form 'Read1PrimerLength Read2PrimerLength'"
            )
o = parser.parse_args()

#Before running, append "End of File" to the end of each .fq file 
#(case sensitive)
#cat "\nEnd of File">>inRead1.fq
#cat "\nEnd of File">>inRead2.fq

spacers = o.spacers.split(' ')
plengths = o.plengths.split(' ')

#Open the various input/output files
inR1 = open(o.inRead1, 'r')
inR2 = open(o.inRead2, 'r')
outName1 = (
        o.out_prefix + ".r1.fq" if o.out_prefix != None 
        else o.inRead1.replace(".fq", ".r1.fq")
        )
outName2 = (
        o.out_prefix + ".r2.fq" if o.out_prefix != None 
        else o.inRead1.replace(".fq", ".r2.fq")
        )
outR1 = open(outName1, 'w')
outR2 = open(outName2, 'w')

#get the first read pair
lineR = [
        [inR1.readline(),inR1.readline(),inR1.readline(), inR1.readline()],
        [inR2.readline(),inR2.readline(),inR2.readline(), inR2.readline()]
        ]

#examine data to see what trimming should do...
def trimReads (read1, read2, trimPos, primeLength): 
	trimmedRead1 = ["",""]
	trimmedRead2 = ["",""]
	#trim the read
	trimmedRead1[0] = read1[0][int(primeLength[0]):trimPos[0]]
	trimmedRead1[1] = read1[1][int(primeLength[0]):trimPos[0]]
	#trim the mate
	trim2 = (int(trimPos[1]) if trimPos[1] != -1 else 101)
	trimmedRead2[0] = read2[0][int(primeLength[0]) + int(primeLength[1]):trim2]
	trimmedRead2[1] = read2[1][int(primeLength[0]) + int(primeLength[1]):trim2]
	#make sure the two are the same length?
	if len(trimmedRead1[0]) != len(trimmedRead2[0]):
		sys.stderr.write("Warning: Read 1 and Read 2 different lengths.")
	return([trimmedRead1, trimmedRead2])

hdrRenameFxn ='lambda x,y:"%s#%s/%s" % (x.split("#")[0], y, x.split("/")[-1])'
hdrRenameFxn = eval (hdrRenameFxn)

while lineR[0][0] != ""  and lineR[1][0] != "" :
	barcode = []
	trimmedReads = []
	sPos = [-1,-1]
	for spacer in spacers:
		if sPos[0] == -1:
			sPos[0] = lineR[0][1].find(spacer) 
		if sPos[1]==-1:
			sPos[1] = lineR[1][1].find(spacer)
	barcode = lineR[0][1][sPos[0] + o.slength:sPos[0] + o.slength + o.blength]
	trimmedReads = trimReads(
            [lineR[0][1], lineR[0][3]], 
            [lineR[1][1], lineR[1][3]], 
            sPos, 
            plengths
            )
	
	#modify the lines
	lineR[0][0] = hdrRenameFxn(lineR[0][0], barcode)
	lineR[0][1] = trimmedReads[0][0]
	lineR[0][3] = trimmedReads[0][1]
	
	lineR[1][0] = hdrRenameFxn(lineR[1][0], barcode)
	lineR[1][1] = trimmedReads[1][0]
	lineR[1][3] = trimmedReads[1][1]
	
	#write the lines to output
	outR1.write(lineR[0][0])
	outR1.write(lineR[0][1] + "\n")
	outR1.write(lineR[0][2])
	outR1.write(lineR[0][3] + "\n")
	
	outR2.write(lineR[1][0])
	outR2.write(lineR[1][1] + "\n")
	outR2.write(lineR[1][2])
	outR2.write(lineR[1][3] + "\n")
	
	#get next set of lines
	lineR = [
            [inR1.readline(),inR1.readline(),inR1.readline(), inR1.readline()],
            [inR2.readline(),inR2.readline(),inR2.readline(), inR2.readline()]
            ]

inR1.close()
inR2.close()
outR1.close()
outR2.close()
