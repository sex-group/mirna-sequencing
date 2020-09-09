#!/usr/bin/python2
# -*- coding: utf-8 -*-

#########################################################################
#                                                                       #
# Add genes that are not present on the main file                       #
#                                                                       #
# Giving two files, one of them with missing genes, this scrip should   #
# add the remaining genes to the original file                          #
#                                                                       #
#########################################################################



# Importing required modules
import os
import sys



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
	print """Usage: python2 2017-phd-jan-09-01.py [OPTIONS]
From two TSV with several genes, in which the second file may have more genes than the first, this script adds the genes that are on the sencon file to the first one.

Required arguments:
 -t: INPUT TEMPLATE FILE - The file that may have genes missing.
 -r: INPUT REFERENCE FILE - The file wich have the genes to me added. This file should have no header!


Other arguments:
 -o: OUTPUT - The name of the output file that should be produced (default: standart output).
 -tp: On the -t TSV file, the number of the column where the gene name is. First column is 0 (default: 2).
 -rp: On the -r TSV file, the number of the column where the gene name is. First column is 0 (default: 0).
 -c: Instruction to complete the line to be added inside double quotation marks (default: "NoExpress\t\t[!rp!0]\t[!rp!2]\t\t[!rp!3]\t[!rp!4]\t\t\t\t\t\t\t\t\t\t\t\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t")
 	[!rp!X] - Xth column of the -rp file.
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



# Get column from t file
# Default is 2
global TP
TP=2
def gettp(i):
	if len(sys.argv) <= i+1:
		error(3,'-tp')

	global TP
	TP=sys.argv[i+1]
	return 1				# Advance 1

# Get column from r file
# Default is 0
global RP
RP=0
def getrp(i):
	if len(sys.argv) <= i+1:
		error(3,'-rp')

	global RP
	RP=sys.argv[i+1]
	return 1				# Advance 1

# Get info about line completion when new sequences are found
# Default is detailed on help!
global C
C="NoExpress\t\t[!rp!0]\t[!rp!2]\t\t\t[!rp!3]\t[!rp!4]\t\t\t\t\t\t\t\t\t\t\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t"
def getc(i):
	if len(sys.argv) <= i+1:
		error(3,'-c')

	global C
	C=sys.argv[i+1]
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



# Get template
def gettemplate(i):
	if len(sys.argv) <= i+1:
		error(3,'-t')

	# Process file 1
	global TEMPLATE
	TEMPLATE=sys.argv[i+1]
	fileval(TEMPLATE,4,False)	# Give error if this file no exists
	return 1					# Advance 1 on the arguments




# Get template
def getreference(i):
	if len(sys.argv) <= i+1:
		error(3,'-r')

	# Process file 1
	global REFERENCE
	REFERENCE=sys.argv[i+1]
	fileval(REFERENCE,4,False)	# Give error if this file no exists
	return 1					# Advance 1 on the arguments




# Get input arguments
AGUMENTS={
	'-h': gethelp,
	'-t': gettemplate,
	'-r': getreference,
	'-l': getlicense,
	'-o': getoutput,
	'-tp': gettp,
	'-rp': getrp,
	'-c': getc
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
REQUIRED=['TEMPLATE','REFERENCE']
for var in REQUIRED:
	if var not in globals():
		error(2,var)



# No errors found so far...
# Start...
print """
Program has started...
"""



# Open template file and import all sequences
print "... loading sequences from template file..."
tfile=open(TEMPLATE,'r')
sequences=set()
l=0
for line in tfile:
	# Add to the set
	fragment=line.strip().split('\t')
	sequences.add(fragment[TP])

	# Count for logging
	l+=1

	# Save to output
	savethisline(line.strip())

tfile.close();

# Checkpoint. How many differente sequences were retrieved
print "... Number of identified lines: "+str(l)+" ..."
print "... Number of different sequences: "+str(len(sequences))+" (includes header if present) ..."



# Loop on the reference file
print "... loading sequences from reference file..."
rfile=open(REFERENCE,'r')
it=0
ia=0
for line in rfile:
	# Check if sequence is present on sequences
	fragment=line.strip().split('\t')
	if fragment[RP] not in sequences:
		# Add this sequence
		# Prepare line
		line=C

		for i in range(len(fragment)):
			line=line.replace("[!rp!"+str(i)+"]",fragment[i])

		# Add it!
		savethisline(line)

		# Count it for logging
		ia+=1

	# Count for logging
	it+=1

rfile.close();

# Checkpoint... How many different sequences were found? And new?
print "... Number of identified lines: "+str(it)+" ..."
print "... Number of sequences added: "+str(ia)+" ..."
print "... Number of lines that the final file should have: "+str(int(l+ia))+" ..."

# The end!
print """
... Program has finished!
"""
