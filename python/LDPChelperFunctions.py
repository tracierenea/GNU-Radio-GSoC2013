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
