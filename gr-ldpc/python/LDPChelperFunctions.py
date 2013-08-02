#!/usr/bin/python
import string
import numpy as np

def writeAlistFile(filename,H,verbose=0):
	# This function reads in a parity check matrix H and writes the
	# corresponding alist file. The format of alist files is
	# desribed at: 
	# http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html

	myfile = open(filename,'w')

	numRows = H.shape[0]
	numCols = H.shape[1]

	tempstring = `numCols` + ' ' + `numRows` + '\n'
	myfile.write(tempstring)

	tempstring1 = ''
	tempstring2 = ''
	maxRowWeight = 0
	for rowNum in np.arange(numRows):
		nonzeros = np.array(H[rowNum,:].nonzero())
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
	for colNum in np.arange(numCols):
		nonzeros = np.array(H[:,colNum].nonzero())
		colWeight = nonzeros.shape[1]
		if colWeight > maxColWeight:
			maxColWeight = colWeight
		tempstring3 = tempstring3 + `colWeight` + ' '
		for tempArray in nonzeros:
			for index in tempArray:
				tempstring4 = tempstring4 + `index+1` + ' '
			tempstring4 = tempstring4 + '\n'
	tempstring3 = tempstring3 + '\n'

	tempstring = `maxColWeight` + ' ' + `maxRowWeight` + '\n'
	myfile.write(tempstring)
	myfile.write(tempstring3) # all the column weights
	myfile.write(tempstring1) # all the row weights
	myfile.write(tempstring4) # the nonzero indices for each column
	myfile.write(tempstring2) # the nonzero indices for each row
	myfile.close()
