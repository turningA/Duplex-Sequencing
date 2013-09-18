import sys
from collections import defaultdict

inFile = sys.stdin

varDict = defaultdict(lambda:0)

for line in inFile :
	splitline = line.split('\t')

	if float(splitline[4])/float(splitline[3]) >= .8 :
		if float(splitline[5])/float(splitline[3]) >= .8 :
			varDict[splitline[1]+'T'] += 1
		elif float(splitline[6])/float(splitline[3]) >= .8 :
			varDict[splitline[1]+'C'] += 1
		elif float(splitline[7])/float(splitline[3]) >= .8 :
			varDict[splitline[1]+'G'] += 1
		elif float(splitline[8])/float(splitline[3]) >= .8 :
			varDict[splitline[1]+'A'] += 1
	else :
		continue

for key in varDict.iterkeys() :
	print '%s\t%s' % (key,varDict[key])