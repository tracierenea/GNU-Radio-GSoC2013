#!/usr/bin/python
import string
from numpy import *

def readAlistFile(filename, verbose=0):
	# This function reads in an alist file and returns the
	# corresponding parity check matrix H. The format of alist files
	# is desribed at:
	# http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html

	myfile = open(filename,'r')
	data = myfile.readlines()

	size      = string.split(data[0])
	numRows   = int(size[0])
	numCols   = int(size[1])
	weights   = string.split(data[1])
	rowWeight = weights[0]
	colWeight = weights[1]

	if verbose: 
		print 'Matrix size: rows:', numRows, 'columns:', numCols
		print 'Row weight:', rowWeight, '/ Column Weight:', colWeight

	H = zeros((numRows,numCols))
	for lineNumber in arange(4,4+numRows):
		indices = string.split(data[lineNumber])
		for index in indices:
			H[lineNumber-4,int(index)-1] = 1

	# the subsequent lines in the file list the indices for where
	# the 1s are in the columns, but this is redundant information

	return H

def writeAlistFile(filename,H,verbose=0):
	# This function reads in a parity check matrix H and writes the
	# corresponding alist file. The format of alist files is
	# desribed at: 
	# http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html

	myfile = open(filename,'w')

	numRows = H.shape[0]
	numCols = H.shape[1]

	if verbose: 
		print 'Writing to file:', filename
		print 'Matrix size: rows:', numRows, 'columns:', numCols

	tempstring = `numRows` + ' ' + `numCols` + '\n'
	myfile.write(tempstring)

	tempstring1 = ''
	tempstring2 = ''
	maxRowWeight = 0
	for rowNum in arange(numRows):
		nonzeros = array(H[rowNum,:].nonzero())
		rowWeight = nonzeros.shape[1]
		if rowWeight > maxRowWeight:
			maxRowWeight = rowWeight
		tempstring1 = tempstring1 + `rowWeight` + ' '
		for tempArray in nonzeros:
			for index in tempArray:
				tempstring2 = tempstring2 + `index+1` + ' '
			tempstring2 = tempstring2 + '\n'
	tempstring1 = tempstring1 + '\n'

	tempstring3 = ''
	tempstring4 = ''
	maxColWeight = 0
	for colNum in arange(numCols):
		nonzeros = array(H[:,colNum].nonzero())
		colWeight = nonzeros.shape[1]
		if colWeight > maxColWeight:
			maxColWeight = colWeight
		tempstring3 = tempstring3 + `colWeight` + ' '
		for tempArray in nonzeros:
			for index in tempArray:
				tempstring4 = tempstring4 + `index+1` + ' '
			tempstring4 = tempstring4 + '\n'
	tempstring3 = tempstring3 + '\n'

	tempstring = `maxRowWeight` + ' ' + `maxColWeight` + '\n'
	myfile.write(tempstring)
	myfile.write(tempstring1) # all the row weights
	myfile.write(tempstring3) # all the column weights
	myfile.write(tempstring2) # the nonzero indices for each row
	myfile.write(tempstring4) # the nonzero indices for each column
	myfile.close()

	if verbose: print 'File closed.'
