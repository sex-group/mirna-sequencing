#!/usr/bin/python2
# -*- coding: utf-8 -*-

#########################################################################
#                                                                       #
# Please... See help function or run -h option!                         #
#                                                                       #
#########################################################################



# Importing required modules
import os
import sys
from Bio import SeqIO



# Define functions
# Functions are listed by alphabetic order

# Error
def error(n,arg):
	ERRORLIST={1: "One of the input argument(s) was not recognised.",
		2: "One of the required parameters was not provided.",
		3: "One of the arguments lack its value.",
		4: "Input file does not exist.",
		5: "Output file already exists!",
		6: "Unable to create output files. Check permissions in current directory.",
		}
	print "[Error#"+str(n)+"] "+ERRORLIST[n]
	print "[Error#"+str(n)+"] "+arg
	print "[Error#"+str(n)+"] Please run -h to get some help."
	sys.exit(n)



# File validation
def fileval(checkfile,en,tf):
	# Check if file exists.
	if os.path.isfile(checkfile) == tf:
		error(en,checkfile)



# Get Help
def gethelp(i):
	print """Usage: python2 2017-phd-jan-12-01.py [OPTIONS]
Creates a comprehensive table (TSV file) using data from several other files. Comprehensive table will include Contig name, Gene name (VIT, Intergenic, ...), miRNA database (miRBase, PMRD or both), miRNA name, miRNA Accession, Organism short name (eg vvi), Chromosome where contig mapped, Starting position of the contig, Ending position of the contig, Gene annotation according databases (functional category), Expression (of all the twelve samples), Differentially expression factor (on all combinations), Sequence of the contig and Sequence of the miRNA.

Required arguments:
 -t: TSV file - The file that has all the genes names and some of the sequences (e.g. an output from 2017-phd-jan-09-01-py). Must have header!
 -bm: Blastn output from miRBase (format 6)
 -bp: Blastn output from PMDR (format 6)
 -d: TSV file with information regarding differentially expressed genes
 -fm: Fasta file that originated the miRBase database
 -fp: Fasta file that originated the PMDR database
 -a: Annotation file - The file that has all the genes names and functional categories. Must have header!

NOTE:
 This program was not created in a way that allows you to define the column where each piece of information is. When running with different inputs, it may have to be edited and the output must be confirmed!

Other arguments:
 -o: OUTPUT - The name of the output file that should be produced. Default: standart output.
 -h: Prints this help and exists
 -l: Prints the license statement and exists

Exit status:
 0 if OK
 N If errors were found and printed on screen. The exit status is the error number.

Report bugs to Miguel Ramos <miguelramos22@gmail.com>

Copyright 2017 Miguel Ramos
This program comes with ABSOLUTELY NO WARRANTY; for details run -l option.
This is free software, and you are welcome to redistribute it under certain conditions."""
	sys.exit(0)



# Get License
def getlicense(i):
	print """Copyright 2017 Miguel Ramos <miguelramos22@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>."""
	sys.exit(0)



# Get output
# Default is standard output
global OUTPUT
OUTPUT=False
def getoutput(i):
	if len(sys.argv) <= i+1:
		error(3,'-o')

	global OUTPUT
	OUTPUT=sys.argv[i+1]
	fileval(OUTPUT,5,True)	# Give error if this file already exists

	# Create this file!
	open(OUTPUT, 'w').close()	# Create the file
	fileval(OUTPUT,6,False)		# Confirm if file was created

	print "Output file created: "+OUTPUT
	print "Data will be written to file..."

	# Save header
	savethisline("ContigName\tGeneID\tmiRNA_Database\tmiRNA_Name\tmiRNA_Accession\tmiRNA_Organism\tChromosome\tContigStartChr\tContigEndChr\tGeneName\tGeneName_Older\tGeneFunctionalCategory\tFB_RPKM\tFD_RPKM\tFG_RPKM\tFH_RPKM\tMB_RPKM\tMD_RPKM\tMG_RPKM\tMH_RPKM\tTNB_RPKM\tTND_RPKM\tTNG_RPKM\tTNH_RPKM\tDif_FB_MB\tDif_FB_TNB\tDif_FD_MD\tDif_FD_TND\tDif_FG_MG\tDif_FG_TNG\tDif_FH_MH\tDif_FH_TNH\tDif_MB_TNB\tDif_MD_TND\tDif_MG_TNG\tDif_MH_TNH\tSeq_miRNA\tSeqContig")

	print ""
	return 1				# Advance 1

# Save this line
def savethisline(line):
	# Open globals
	global OUTPUT

	# Save or print
	if OUTPUT is False:
		print line
	else:
		ofile=open(OUTPUT, 'a')
		ofile.write(line+"\n")
		ofile.close()



# Get -t
def gettsvfile(i):
	if len(sys.argv) <= i+1:
		error(3,'-t')

	# Process file 1
	global TSV
	TSV=sys.argv[i+1]
	fileval(TSV,4,False)	# Give error if this file no exists

	# Inform and load information
	print "TSV file located: "+TSV
	print "Loading..."

	global TSV_CONTENT
	TSV_CONTENT=dict()
	with open(TSV) as f:
		# This file contains header
		next(f)
		for line in f:
			fragment=line.strip().split('\t')
			if fragment[1]!="": contig=fragment[1]
			else: contig=fragment[2]
			TSV_CONTENT[contig]=fragment

	print ""
	return 1				# Advance 1 on the arguments

# Get -bm
def getblastmirbase(i):
	if len(sys.argv) <= i+1:
		error(3,'-bm')

	# Process file 1
	global BLASTNMIRBASE
	BLASTNMIRBASE=sys.argv[i+1]
	fileval(BLASTNMIRBASE,4,False)	# Give error if this file no exists

	# Inform and load information
	print "miRBase blastn file located: "+BLASTNMIRBASE
	print "Loading..."

	global BLASTNMIRBASE_CONTENT
	BLASTNMIRBASE_CONTENT=dict()
	with open(BLASTNMIRBASE) as f:
		# This file does not contain header
		for line in f:
			fragment=line.strip().split('\t')
			contig=fragment[0].split('_G_')
			if len(contig)==2: contig='G_'+contig[1]
			else: contig=contig[0]

			# Check if contig was already added
			if contig in BLASTNMIRBASE_CONTENT:
				# It was already added. Replace?
				# Check both evalues and if its lower, print it!
				if float(fragment[10])<float(BLASTNMIRBASE_CONTENT[contig][10]): BLASTNMIRBASE_CONTENT[contig]=fragment
			else: BLASTNMIRBASE_CONTENT[contig]=fragment

	print ""
	return 1					# Advance 1 on the arguments

# Get -bp
def getblastpmrd(i):
	if len(sys.argv) <= i+1:
		error(3,'-bp')

	# Process file 1
	global BLASTNPMRD
	BLASTNPMRD=sys.argv[i+1]
	fileval(BLASTNPMRD,4,False)	# Give error if this file no exists

		# Inform and load information
	print "PMRD blastn file located: "+BLASTNPMRD
	print "Loading..."

	global BLASTNPMRD_CONTENT
	BLASTNPMRD_CONTENT=dict()
	with open(BLASTNPMRD) as f:
		# This file does not contain header
		for line in f:
			fragment=line.strip().split('\t')
			contig=fragment[0].split('_G_')
			if len(contig)==2: contig='G_'+contig[1]
			else: contig=contig[0]

			# Check if contig was already added
			if contig in BLASTNPMRD_CONTENT:
				# It was already added. Replace?
				# Check both evalues and if its lower, print it!
				if float(fragment[10])<float(BLASTNPMRD_CONTENT[contig][10]): BLASTNPMRD_CONTENT[contig]=fragment
			else: BLASTNPMRD_CONTENT[contig]=fragment

	print ""
	return 1					# Advance 1 on the arguments

# Get -d
def getdiffexpress(i):
	if len(sys.argv) <= i+1:
		error(3,'-d')

	# Process file 1
	global DIFEXPRESSED
	DIFEXPRESSED=sys.argv[i+1]
	fileval(DIFEXPRESSED,4,False)	# Give error if this file no exists

	# Inform and load information
	print "Differentially expressed genes file located: "+DIFEXPRESSED
	print "Loading..."

	global DIFEXPRESSED_CONTENT
	DIFEXPRESSED_CONTENT=dict()
	with open(DIFEXPRESSED) as f:
		# This file contains header
		next(f)
		for line in f:
			fragment=line.strip().split('\t')
			DIFEXPRESSED_CONTENT[fragment[1]]=fragment

	print ""
	return 1					# Advance 1 on the arguments

# Get -fm
def getfastamirbase(i):
	if len(sys.argv) <= i+1:
		error(3,'-fm')

	# Process file 1
	global FASTAMIRBASE
	FASTAMIRBASE=sys.argv[i+1]
	fileval(FASTAMIRBASE,4,False)	# Give error if this file no exists

	# Inform and load information
	print "Fasta file from miRBase located: "+FASTAMIRBASE
	print "Loading..."

	global FASTAMIRBASE_CONTENT
	FASTAMIRBASE_CONTENT=SeqIO.to_dict(SeqIO.parse(FASTAMIRBASE, "fasta"))

	print ""
	return 1					# Advance 1 on the arguments

# Get -fp
def getfastapmrd(i):
	if len(sys.argv) <= i+1:
		error(3,'-fp')

	# Process file 1
	global FASTAPMRD
	FASTAPMRD=sys.argv[i+1]
	fileval(FASTAPMRD,4,False)	# Give error if this file no exists

	# Inform and load information
	print "Fasta file from PMRD located: "+FASTAPMRD
	print "Loading..."

	global FASTAPMRD_CONTENT
	FASTAPMRD_CONTENT=SeqIO.index(FASTAPMRD, "fasta")

	print ""
	return 1					# Advance 1 on the arguments

# Get -a
def getannotation(i):
	if len(sys.argv) <= i+1:
		error(3,'-t')

	# Process file 1
	global ANNOTATION
	ANNOTATION=sys.argv[i+1]
	fileval(ANNOTATION,4,False)	# Give error if this file no exists

	# Inform and load information
	print "Annotation file located: "+ANNOTATION
	print "Loading..."

	global ANNOTATION_CONTENT
	ANNOTATION_CONTENT=dict()
	with open(ANNOTATION) as f:
		# This file contains header
		next(f)
		for line in f:
			fragment=line.strip().split('\t')
			ANNOTATION_CONTENT[fragment[0]]=fragment

	print ""
	return 1				# Advance 1 on the arguments



# Start...
print """
Program has started...
"""



# Get input arguments
AGUMENTS={
	'-h': gethelp,
	'-l': getlicense,
	'-t': gettsvfile,
	'-bm': getblastmirbase,
	'-bp': getblastpmrd,
	'-d': getdiffexpress,
	'-fm': getfastamirbase,
	'-fp': getfastapmrd,
	'-a': getannotation,
	'-o': getoutput
	}

n=0
for i in range(1,len(sys.argv)):
	if n>0:
		n=n-1
		continue
	if sys.argv[i] in AGUMENTS.keys():
		n=AGUMENTS[sys.argv[i]](i)
	else:
		error(1,sys.argv[i])



# Check if all arguments were added
# Required variables
REQUIRED=['TSV','BLASTNMIRBASE','BLASTNPMRD','DIFEXPRESSED','FASTAMIRBASE','FASTAPMRD','ANNOTATION']
for var in REQUIRED:
	if var not in globals():
		error(2,var)

# To do...
# Start from miRBase
for hit in BLASTNMIRBASE_CONTENT:
	# See if hit was also found on PMRD
	if hit in BLASTNPMRD_CONTENT:
		if float(BLASTNMIRBASE_CONTENT[hit][10])<float(BLASTNPMRD_CONTENT[hit][10]):
			database="Both:miRBase"
			idatabase="miRBase"
		else:
			database="Both:PMRD"
			idatabase="pmrd"
	else:
		database="miRBase"
		idatabase="miRBase"

	# Get ContigName
	ContigName=hit.split(".")
	if len(ContigName)==1: ContigName=ContigName[0]
	else: ContigName=""

	# Get GeneID
	GeneID=BLASTNMIRBASE_CONTENT[hit][0].split("_G_")[0]
	VIT=GeneID.split(".")[0]

	# Get miRNAName, miRNAAccession, miRNA_Organism and Seq_miRNA
	if idatabase=="miRBase":
		miRNA_Name=BLASTNMIRBASE_CONTENT[hit][1]
		miRNA_Description=FASTAMIRBASE_CONTENT[miRNA_Name].description.split(" ")
		miRNA_Accession=miRNA_Description[1]
		miRNA_Organism=miRNA_Description[2]+" "+miRNA_Description[3]

		Seq_miRNA=str(FASTAMIRBASE_CONTENT[miRNA_Name].seq)
	else:
		miRNA_Name=BLASTNPMRD_CONTENT[hit][1]
		miRNA_Accession=""
		miRNA_Organism=miRNA_Name.split("-")[0]

		Seq_miRNA=str(FASTAPMRD_CONTENT[miRNA_Name].seq)

	# Get Chromossome, RPKMs and sequence
	if ContigName in TSV_CONTENT:
		Chromosome=TSV_CONTENT[ContigName][3]
		ContigStartChr=TSV_CONTENT[ContigName][6]
		ContigEndChr=TSV_CONTENT[ContigName][7]

		FB_RPKM=TSV_CONTENT[ContigName][18]
		FD_RPKM=TSV_CONTENT[ContigName][19]
		FG_RPKM=TSV_CONTENT[ContigName][20]
		FH_RPKM=TSV_CONTENT[ContigName][21]

		MB_RPKM=TSV_CONTENT[ContigName][22]
		MD_RPKM=TSV_CONTENT[ContigName][23]
		MG_RPKM=TSV_CONTENT[ContigName][24]
		MH_RPKM=TSV_CONTENT[ContigName][25]

		TNB_RPKM=TSV_CONTENT[ContigName][26]
		TND_RPKM=TSV_CONTENT[ContigName][27]
		TNG_RPKM=TSV_CONTENT[ContigName][28]
		TNH_RPKM=TSV_CONTENT[ContigName][29]

		SeqContig=TSV_CONTENT[ContigName][30]

	else:
		Chromosome=TSV_CONTENT[VIT][3]
		ContigStartChr=TSV_CONTENT[VIT][6]
		ContigEndChr=TSV_CONTENT[VIT][7]

		FB_RPKM=TSV_CONTENT[VIT][18]
		FD_RPKM=TSV_CONTENT[VIT][19]
		FG_RPKM=TSV_CONTENT[VIT][20]
		FH_RPKM=TSV_CONTENT[VIT][21]

		MB_RPKM=TSV_CONTENT[VIT][22]
		MD_RPKM=TSV_CONTENT[VIT][23]
		MG_RPKM=TSV_CONTENT[VIT][24]
		MH_RPKM=TSV_CONTENT[VIT][25]

		TNB_RPKM=TSV_CONTENT[VIT][26]
		TND_RPKM=TSV_CONTENT[VIT][27]
		TNG_RPKM=TSV_CONTENT[VIT][28]
		TNH_RPKM=TSV_CONTENT[VIT][29]

		SeqContig=""

	# Get GeneName, GeneName_Older and GeneFunctionalCategory
	GeneName=""
	GeneName_Older=""
	GeneFunctionalCategory=""

	if VIT in ANNOTATION_CONTENT:
		if len(ANNOTATION_CONTENT[VIT])>=4: GeneName=ANNOTATION_CONTENT[VIT][3]
		if len(ANNOTATION_CONTENT[VIT])>=5: GeneName_Older=ANNOTATION_CONTENT[VIT][4]
		if len(ANNOTATION_CONTENT[VIT])>=6: GeneFunctionalCategory=ANNOTATION_CONTENT[VIT][5]

	# Get Dif Expression index!
	if ContigName in DIFEXPRESSED_CONTENT:
		Dif_FB_MB=DIFEXPRESSED_CONTENT[ContigName][19]
		Dif_FB_TNB=DIFEXPRESSED_CONTENT[ContigName][20]
		Dif_FD_MD=DIFEXPRESSED_CONTENT[ContigName][21]
		Dif_FD_TND=DIFEXPRESSED_CONTENT[ContigName][22]
		Dif_FG_MG=DIFEXPRESSED_CONTENT[ContigName][23]
		Dif_FG_TNG=DIFEXPRESSED_CONTENT[ContigName][24]
		Dif_FH_MH=DIFEXPRESSED_CONTENT[ContigName][25]
		Dif_FH_TNH=DIFEXPRESSED_CONTENT[ContigName][26]
		Dif_MB_TNB=DIFEXPRESSED_CONTENT[ContigName][27]
		Dif_MD_TND=DIFEXPRESSED_CONTENT[ContigName][28]
		Dif_MG_TNG=DIFEXPRESSED_CONTENT[ContigName][29]
		Dif_MH_TNH=DIFEXPRESSED_CONTENT[ContigName][30]
	else:
		Dif_FB_MB=""
		Dif_FB_TNB=""
		Dif_FD_MD=""
		Dif_FD_TND=""
		Dif_FG_MG=""
		Dif_FG_TNG=""
		Dif_FH_MH=""
		Dif_FH_TNH=""
		Dif_MB_TNB=""
		Dif_MD_TND=""
		Dif_MG_TNG=""
		Dif_MH_TNH=""

	# If it is also present on PMRD, delete it from there! Even if PMRD was not used (if it was best, we would have used it!)
	if hit in BLASTNPMRD_CONTENT: del BLASTNPMRD_CONTENT[hit]


	# Save this line
	savethisline(ContigName+"\t"+GeneID+"\t"+database+"\t"+miRNA_Name+"\t"+miRNA_Accession+"\t"+miRNA_Organism+"\t"+Chromosome+"\t"+ContigStartChr+"\t"+ContigEndChr+"\t"+GeneName+"\t"+GeneName_Older+"\t"+GeneFunctionalCategory+"\t"+FB_RPKM+"\t"+FD_RPKM+"\t"+FG_RPKM+"\t"+FH_RPKM+"\t"+MB_RPKM+"\t"+MD_RPKM+"\t"+MG_RPKM+"\t"+MH_RPKM+"\t"+TNB_RPKM+"\t"+TND_RPKM+"\t"+TNG_RPKM+"\t"+TNH_RPKM+"\t"+Dif_FB_MB+"\t"+Dif_FB_TNB+"\t"+Dif_FD_MD+"\t"+Dif_FD_TND+"\t"+Dif_FG_MG+"\t"+Dif_FG_TNG+"\t"+Dif_FH_MH+"\t"+Dif_FH_TNH+"\t"+Dif_MB_TNB+"\t"+Dif_MD_TND+"\t"+Dif_MG_TNG+"\t"+Dif_MH_TNH+"\t"+Seq_miRNA+"\t"+SeqContig)





# Go to PMRD
# All mirbase occurences have been processed. These are exlussive
for hit in BLASTNPMRD_CONTENT:
	# See if hit was also found on PMRD
	database="PMRD"

	# Get ContigName
	ContigName=hit.split(".")
	if len(ContigName)==1: ContigName=ContigName[0]
	else: ContigName=""

	# Get GeneID
	GeneID=BLASTNPMRD_CONTENT[hit][0].split("_G_")[0]
	VIT=GeneID.split(".")[0]

	# Get miRNAName, miRNAAccession, miRNA_Organism and Seq_miRNA
	miRNA_Name=BLASTNPMRD_CONTENT[hit][1]
	miRNA_Accession=""
	miRNA_Organism=miRNA_Name.split("-")[0]

	Seq_miRNA=str(FASTAPMRD_CONTENT[miRNA_Name].seq)

	# Get Chromossome, RPKMs and sequence
	if ContigName in TSV_CONTENT:
		Chromosome=TSV_CONTENT[ContigName][3]
		ContigStartChr=TSV_CONTENT[ContigName][6]
		ContigEndChr=TSV_CONTENT[ContigName][7]

		FB_RPKM=TSV_CONTENT[ContigName][18]
		FD_RPKM=TSV_CONTENT[ContigName][19]
		FG_RPKM=TSV_CONTENT[ContigName][20]
		FH_RPKM=TSV_CONTENT[ContigName][21]

		MB_RPKM=TSV_CONTENT[ContigName][22]
		MD_RPKM=TSV_CONTENT[ContigName][23]
		MG_RPKM=TSV_CONTENT[ContigName][24]
		MH_RPKM=TSV_CONTENT[ContigName][25]

		TNB_RPKM=TSV_CONTENT[ContigName][26]
		TND_RPKM=TSV_CONTENT[ContigName][27]
		TNG_RPKM=TSV_CONTENT[ContigName][28]
		TNH_RPKM=TSV_CONTENT[ContigName][29]

		SeqContig=TSV_CONTENT[ContigName][30]

	else:
		Chromosome=TSV_CONTENT[VIT][3]
		ContigStartChr=TSV_CONTENT[VIT][6]
		ContigEndChr=TSV_CONTENT[VIT][7]

		FB_RPKM=TSV_CONTENT[VIT][18]
		FD_RPKM=TSV_CONTENT[VIT][19]
		FG_RPKM=TSV_CONTENT[VIT][20]
		FH_RPKM=TSV_CONTENT[VIT][21]

		MB_RPKM=TSV_CONTENT[VIT][22]
		MD_RPKM=TSV_CONTENT[VIT][23]
		MG_RPKM=TSV_CONTENT[VIT][24]
		MH_RPKM=TSV_CONTENT[VIT][25]

		TNB_RPKM=TSV_CONTENT[VIT][26]
		TND_RPKM=TSV_CONTENT[VIT][27]
		TNG_RPKM=TSV_CONTENT[VIT][28]
		TNH_RPKM=TSV_CONTENT[VIT][29]

		SeqContig=""

	# Get GeneName, GeneName_Older and GeneFunctionalCategory
	GeneName=""
	GeneName_Older=""
	GeneFunctionalCategory=""

	if VIT in ANNOTATION_CONTENT:
		if len(ANNOTATION_CONTENT[VIT])>=4: GeneName=ANNOTATION_CONTENT[VIT][3]
		if len(ANNOTATION_CONTENT[VIT])>=5: GeneName_Older=ANNOTATION_CONTENT[VIT][4]
		if len(ANNOTATION_CONTENT[VIT])>=6: GeneFunctionalCategory=ANNOTATION_CONTENT[VIT][5]

	# Get Dif Expression index!
	if ContigName in DIFEXPRESSED_CONTENT:
		Dif_FB_MB=DIFEXPRESSED_CONTENT[ContigName][19]
		Dif_FB_TNB=DIFEXPRESSED_CONTENT[ContigName][20]
		Dif_FD_MD=DIFEXPRESSED_CONTENT[ContigName][21]
		Dif_FD_TND=DIFEXPRESSED_CONTENT[ContigName][22]
		Dif_FG_MG=DIFEXPRESSED_CONTENT[ContigName][23]
		Dif_FG_TNG=DIFEXPRESSED_CONTENT[ContigName][24]
		Dif_FH_MH=DIFEXPRESSED_CONTENT[ContigName][25]
		Dif_FH_TNH=DIFEXPRESSED_CONTENT[ContigName][26]
		Dif_MB_TNB=DIFEXPRESSED_CONTENT[ContigName][27]
		Dif_MD_TND=DIFEXPRESSED_CONTENT[ContigName][28]
		Dif_MG_TNG=DIFEXPRESSED_CONTENT[ContigName][29]
		Dif_MH_TNH=DIFEXPRESSED_CONTENT[ContigName][30]
	else:
		Dif_FB_MB=""
		Dif_FB_TNB=""
		Dif_FD_MD=""
		Dif_FD_TND=""
		Dif_FG_MG=""
		Dif_FG_TNG=""
		Dif_FH_MH=""
		Dif_FH_TNH=""
		Dif_MB_TNB=""
		Dif_MD_TND=""
		Dif_MG_TNG=""
		Dif_MH_TNH=""

	# Save this line
	savethisline(ContigName+"\t"+GeneID+"\t"+database+"\t"+miRNA_Name+"\t"+miRNA_Accession+"\t"+miRNA_Organism+"\t"+Chromosome+"\t"+ContigStartChr+"\t"+ContigEndChr+"\t"+GeneName+"\t"+GeneName_Older+"\t"+GeneFunctionalCategory+"\t"+FB_RPKM+"\t"+FD_RPKM+"\t"+FG_RPKM+"\t"+FH_RPKM+"\t"+MB_RPKM+"\t"+MD_RPKM+"\t"+MG_RPKM+"\t"+MH_RPKM+"\t"+TNB_RPKM+"\t"+TND_RPKM+"\t"+TNG_RPKM+"\t"+TNH_RPKM+"\t"+Dif_FB_MB+"\t"+Dif_FB_TNB+"\t"+Dif_FD_MD+"\t"+Dif_FD_TND+"\t"+Dif_FG_MG+"\t"+Dif_FG_TNG+"\t"+Dif_FH_MH+"\t"+Dif_FH_TNH+"\t"+Dif_MB_TNB+"\t"+Dif_MD_TND+"\t"+Dif_MG_TNG+"\t"+Dif_MH_TNH+"\t"+Seq_miRNA+"\t"+SeqContig)

# The end!
print """
... Program has finished!
"""
