Duplex-Sequencing
=================
README

** Duplex Sequencing software package **
** Version 1.21 **
** February 5, 2013 **
** Programs by Scott Kennedy and Mike Schmitt **
** Several steps are based on prior work by Joe Hiatt **

1. DEPENDENCIES

samtools and pysam MUST be installed on your computer for the scripts to work.

http://http://samtools.sourceforge.net/
http://code.google.com/p/pysam/

2. REQUIRED PATHS

Set the DCSPATH and REFPATH variables to point to the locations of your DCS programs, and reference files, respectively.

--Use your favorite text editor (e.g. nano) to open the .bash_profile file in your home directory and add the following lines (this assumes that you are putting the DCS folder in your root directory. Change the path as needed):

DCSPATH= /DCS/programs; export DCSPATH
REFPATH = /DCS/reference; export REFPATH

--Save the file and close the editor.
--Run the following command:

source $HOME/.bash_profile

3. USAGE

--Run the dcs-script.txt file from the directory containing your fastq files. The fastq files must be named seq1.fq and seq2.fq  Recommended command-line usage:

dcs-script.txt 2> dcs-script.se

--The file dcs-script-parallel.txt can also be used. This version processes the read 1 and read 2 files simultaneously. It is faster but requires more memory (minimum 16 GB recommended).

Recommended command-line usage:

dcs-script-parallel.txt 2> dcs-script-parallel.se

4. DATA OUTPUT

--BAM file consisting of DCS reads: seq_both_DCS.bam
--BAM file consisting of SSCS reads: seq_both_sscs.bam

***note that these SSCS reads are all aligned relative to the reference genome, and have thus been reverse-complemented when necessary by the aligner. Thus these SSCS reads do NOT inform whether there is a strand bias due to DNA damage. Doing so requires looking at forward-mapping and reverse-mapping reads separately after the initial alignment. We intend to automate this type of analysis in a future version of our software.

--We have noticed that alignment errors at the ends of reads can result in false mutations. To eliminate these, we also hard-clip the first and last 5 nt of each read after alignment: seq_both_DCS_readgroups_clipped.bam

--text file listing overall mutation frequencies: these are the files having extension .pileup.countmuts

--text file listing position-specific mutation frequencies: files having extension .pileup.mutpos

***the mutpos file is tab-delimited. Output is: reference name, reference base, position number, depth, number of total mutations (excluding indels), number of mutations to T, C, G, A, insertions, deletions

5. QUALITY CONTROL

It is highly recommended to calculate read count statistics from each run for troubleshooting purposes. To do so, run the dcs-stats.txt script from the directory which contains your files:

dcs-stats.txt

Note that this script is designed specifically for mitochondrial DNA data that was generated with the dcs-script-parallel.txt workflow.

You will want to consider the following numbers. If you lose a lot of data at a single step, you can then troubleshoot that step. For example if you lose a lot of reads going from 'mapped' to 'SSCS', the DNA is probably over-duplicated. Consider re-prepping the DNA using a larger amount of input into the PCR.

* total reads
* reads passing filtering criteria (i.e. with a properly located spacer, and with no failed reads within the Duplex Tag sequence)
* total mapped reads
* reads mapped to the mitochondrial genome
* SSCS reads
* DCS reads
* mapped reads to SSCS ratio
* SSCS to DCS ratio

6. DATA ANALYSIS

--The countmuts file can be read directly in a text viewer

--to plot mutation positions in R (with zeros removed):
read.table("/path/seq_both_DCS.bam.pileup.mutpos") -> table
mutcount <- replace(table$V5,table$V5==0, NA)
plot(mutcount/table$V4)
