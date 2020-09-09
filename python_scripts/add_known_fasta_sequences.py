#!/usr/bin/python2
# -*- coding: utf-8 -*-

#########################################################################
#                                                                       #
# Gets fasta sequences                                                  #
#                                                                       #
# Using the output of 2017-phd-jan-09-01.py and a fasta file with the   #
# reference sequences, it provides a fasta file with all the sequences. #
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
	print """Usage: python2 2017-phd-jan-10-01.py [OPTIONS]
When providing a TSV with sequences, and a fasta file, this script will generate a fasta file with the fasta sequences from the TSV when possible, and on the other cases, it will provide the reference sequences.

Required arguments:
 -t: TSV file - The file that has all the genes names and some of the sequences (e.g. an output from 2017-phd-jan-09-01-py). Must have header!
 -f: REFERENCE FASTA FILE - A fasta file that should have all sequences that the -t does not.


Other arguments:
 -o: OUTPUT - The name of the output file that should be produced. Default: standart output.
 -cp: The number of the column of the -t file that has the contig name. First column is 0. Default: 1.
 -gp: The number of the column of the -t file that has the gene name. First column is 0. Default: 2.
 -sp: The number of the column of the -t file that has the sequence. First column is 0. Default: 30.
 -a: When present, this option will require absolute match between the gene name and the sequence name. If not set, all isoforms will be added to the file!
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

	return 1				# Advance 1



# Get cp - column where gene name is on -t file
global CP
CP=1
def getcp(i):
	if len(sys.argv) <= i+1:
		error(3,'-tp')

	global CP
	CP=sys.argv[i+1]
	return 1				# Advance 1

# Get gp - column where gene name is on -t file
global GP
GP=2
def getgp(i):
	if len(sys.argv) <= i+1:
		error(3,'-tp')

	global GP
	GP=sys.argv[i+1]
	return 1				# Advance 1

# Set -a (absolute)
global ABSOLUTE
ABSOLUTE=False
def seta(i):
	global ABSOLUTE
	ABSOLUTE=True			# Set on
	return 1				# Advance 1

# Get sp - column where sequence should be on -t file
global SP
SP=30
def getsp(i):
	if len(sys.argv) <= i+1:
		error(3,'-rp')

	global SP
	SP=sys.argv[i+1]
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
def gettsv(i):
	if len(sys.argv) <= i+1:
		error(3,'-t')

	# Process file 1
	global TSV
	TSV=sys.argv[i+1]
	fileval(TSV,4,False)	# Give error if this file no exists
	return 1				# Advance 1 on the arguments




# Get -f
def getfasta(i):
	if len(sys.argv) <= i+1:
		error(3,'-f')

	# Process file 1
	global FASTA
	FASTA=sys.argv[i+1]
	fileval(FASTA,4,False)	# Give error if this file no exists
	return 1					# Advance 1 on the arguments




# Get input arguments
AGUMENTS={
	'-h': gethelp,
	'-t': gettsv,
	'-f': getfasta,
	'-l': getlicense,
	'-o': getoutput,
	'-cp': getcp,
	'-gp': getgp,
	'-sp': getsp,
	'-a': seta
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
REQUIRED=['TSV','FASTA']
for var in REQUIRED:
	if var not in globals():
		error(2,var)



# No errors found so far...
# Start...
print """
Program has started...
"""

# Open tsv file and start extracting sequences...
print "... loading sequences from TSV file..."
tfile=open(TSV,'r')
next(tfile)
l=0
lwd=0
lwf=0
for line in tfile:
	# Check if it have sequence
	fragment=line.strip().split('\t')
	if len(fragment)>SP:
		if fragment[GP]!="": savethisline(">"+fragment[GP]+"_"+fragment[CP]+"\n"+fragment[SP])
		else: savethisline(">"+fragment[CP]+"\n"+fragment[SP])
		lwd+=1

	else:
		# Get this sequence from the fasta file...
		for record in SeqIO.parse(FASTA, "fasta"):
			if record.id.startswith(fragment[GP]+"."):
				savethisline(">"+str(record.id)+"\n"+str(record.seq))
				lwf+=1

	# Count for logging
	l+=1

tfile.close();

# Checkpoint. How many differente sequences were retrieved
print "... Number of identified sequences: "+str(l)+" ..."
print "... Number of retrieved fasta sequences (direct): "+str(lwd)+" ..."
print "... Number of retrieved fasta sequences (from reference): "+str(lwf)+" ..."
print "... Number of retrieved fasta sequences (total): "+str(int(lwd+lwd))+" ..."

# The end!
print """
... Program has finished!
"""
