#!/usr/bin/python

from numpy import *
from numpy.random import shuffle, randint
from numpy.linalg import inv, det

verbose = 1

# FIXME - I need to go in an specify argument size constraints
# to ensure that the matrix multiplication works correctly. 

def matrixMultiplierEncoder(H,A):

	# Using notation from Ziemer & Tranter's Principles of 
	# Communications, 5th edition
	size   = H.shape
	r      = size[0]
	k      = size[1]/2
	Hp     = H[0:r,0:k]
	Ir     = identity(r,int)
	Ik     = identity(k,int)
	G      = concatenate((Ik,Hp))
	Hcheck = concatenate((Hp,Ir),axis=1)

	if verbose: printEncodeMatrices(H,Hp,Ir,Ik,G,Hcheck)
	T = dot(G,A) % 2  # codeword vector T. use mod 2 operations
	return T

def printEncodeMatrices(H,Hp,Ir,Ik,G,Hcheck):
	print 'H (parity-check matrix): \n', H
	print '\nHp:\n', Hp
	print '\nIr:\n', Ir
	print '\nIk:\n', Ik
	print '\nG (generator matrix):\n', G
	print '\nHcheck (should be the same):\n', Hcheck , '\n'

def calcSyndrome(H,codeword):
	syndrome = dot(H,codeword) % 2	# use modulo 2 operations
	return syndrome

def haveMatch(syndrome):
	check = 1
	if any(syndrome):
		# the matrix is not all zeros, it's not a codeword
		check = 0
	return check

def singleParityErrorFix(H,testWord):
	# this is from section 10.3.6 of Principles of Communications,
	# 5th ed (Wiley)

	# create a copy since argument testWord array is mutable
	codeword = testWord.copy() 
	syndrome = calcSyndrome(H,codeword)
	if haveMatch(syndrome): 
		print '  - Received word doesn\'t seem to be erroneous.',
		print 'No changes being made.'
		return codeword

	[numRows, numColumns] = H.shape

	#check syndrome against the columns in H.
	for columnNum in arange(numColumns):
		if array_equiv(H[:,columnNum], syndrome.transpose()):
			# now flip the associated symbol bit.
			# i'm sure there is a better way with ...
			if   codeword[columnNum] == 0: codeword[columnNum] = 1
			elif codeword[columnNum] == 1: codeword[columnNum] = 0
			break

	return codeword

def bitFlipDecoder(maxIterations, H, codeword):
	receivedCodeword = codeword
	sizeOfReceivedCodeword = receivedCodeword.shape
	n = sizeOfReceivedCodeword[0] # number of symbols/bit in codeword

	syndrome = calcSyndrome(H,receivedCodeword)
	if haveMatch(syndrome):
		if verbose: print 'Valid codeword. No bit flips required.'
		return receivedCodeword
	else:
		if verbose: print 'Evaluating codeword:', receivedCodeword
		testCodeword = receivedCodeword
		iteration = 1;
	
	while iteration <= maxIterations:
		if verbose: print 'Iteration:', iteration

		# For each of the n bits in the codeword, determine how many
		# of the unsatisfied parity checks involve that bit. To do 
		# this: 
		# First find the nonzero entries in the syndrome. The entry 
		# numbers correspond to rows of interest in H.
		rowsToLookAtInH = array(nonzero(syndrome)) 

		# Second, for each bit, determine how many of unsatisfied 
		# parity checks involve this bit and store this count in an 
		# array. 
		counts = zeros_like(receivedCodeword)

		for row in rowsToLookAtInH.transpose():
			for bitNumber in arange(n):
				if H[row,bitNumber] > 0:
					counts[bitNumber] = counts[bitNumber] + 1

		# Next, determine which bit(s) is associated with the most 
		# unsatisfied parity checks, and flip it/them
		if verbose: print 'counts:', counts
		bitsToFlip = where(counts==counts.max())
		for bitNumber in bitsToFlip[0]:
			if verbose: print 'We need to flip bit:', bitNumber
			testCodeword[bitNumber] = bitwise_xor\
			                          (testCodeword[bitNumber],1)

		if verbose: print 'New codeword candidate: ',testCodeword
		syndrome = calcSyndrome(H,testCodeword)
		if haveMatch(syndrome):
			if verbose: 
				print 'Codeword declared to be:', testCodeword
			return testCodeword
		else:
			iteration = iteration + 1

	else:
		if verbose: 
			print 'Max iteration count of', maxIterations,
			print 'has been reached without finding valid codeword.',
			print 'Returning received codeword.'
		return receivedCodeword


def regularLDPCcodeConstructor(n,p,q):
	# n = codeword length
	# p = column weight
	# q = row weight

	# Following Gallager's approach where we create p submatrices. 
	# Reference: Turbo Coding for Satellite and Wireless 
	# Communications, sec 9.3

	# for this algorithm, n/p must be an integer, because we need the
	# number of rows in eacn submatrix to be a whole number
	ratioTest = (n*1.0)/q
	if ratioTest%1 != 0:
		print 'The ratio of inputs n/q must be a whole number.'
		return

	# FIX: There should probably be other guidelines for n/p/q,
	# but I have not found any specifics in the literature....

	# First one first: 
	m = (n*p)/q  # number of rows in H matrix
	submatrix1 = zeros((m/p,n))  

	for row in arange(m/p):

		# these lines have the ones going from bottom to top, left
		# to right across the matrix. Switching to top to bottom, 
		# left to right, because I think it will be faster to then
		# convert it to systematic form
		# range1 = n-(row+1)*q
		# range2 = n-row*q

		range1 = row*q
		range2 = (row+1)*q 
		submatrix1[row,range1:range2] = 1

	H = submatrix1

	# create the other submatrices and vertically stack them on
	submatrixNum = 2
	newColumnOrder = arange(n)

	while submatrixNum <= p:
		submatrix = zeros((m/p,n))
		shuffle(newColumnOrder)

		for columnNum in arange(n):
			submatrix[:,columnNum] = \
			                  submatrix1[:,newColumnOrder[columnNum]]

		H = vstack((H,submatrix))

		submatrixNum = submatrixNum + 1 

	# double check the row weight and column weights
	size = H.shape
	rows = size[0]
	cols = size[1]

	# check the row weights
	for rowNum in arange(rows):
		nonzeros = array(H[rowNum,:].nonzero())
		if nonzeros.shape[1] != q:
			print 'Row', rowNum, 'has incorrect weight!'
			return

	# check the column weights
	for columnNum in arange(cols):
		nonzeros = array(H[:,columnNum].nonzero())
		if nonzeros.shape[1] != p:
			print 'Row', columnNum, 'has incorrect weight!'
			return

	# At this point, we can use this matrix for parity checks, 
	# but it is not in systematic form so it can't be used
	# to find the parity matrix P, and therefore can't be used
	# for encoding via the generation matrix. 

	# Putting H into its systematic form requires row operations. 
	# I cannot figure out how to write this into an algorithm for
	# matrices of unknown size, and also using mod 2 operations.
	# <in work....>

	# There are some published steps that require taking the 
	# determinant but taking the determinant of a huge matrix
	# won't work.....

	return H

def greedyUpperTriangulation(H):
	# This function performs row/column permutations to bring
	# H into approximate upper triangular form via greedy 
	# upper triangulation method outlined in Modern Coding 
	# Theory Appendix 1, Section A.2
	H_t  = H.copy()
	size = H_t.shape
	n = size[1]
	k = n - size[0]
	g = t = 0
	while t != (n-k-g):
		H_residual = H_t[t:n-k-g,t:n]
		size       = H_residual.shape
		numRows    = size[0]
		numCols    = size[1]

		minResidualDegrees = zeros((1,numCols))

		for colNum in arange(numCols):
			nonZeroElements = array(H_residual[:,colNum].nonzero())
			minResidualDegrees[0,colNum] = nonZeroElements.shape[1]

		# find the minimum nonzero residual degree
		nonZeroElementIndices = minResidualDegrees.nonzero()
		nonZeroElements=minResidualDegrees[nonZeroElementIndices[0],\
		                                   nonZeroElementIndices[1]]
		minimumResidualDegree = nonZeroElements.min()

		# get indices of all of the columns in H_t that have degree
		# equal to the min positive residual degree, then pick at
		# random column c
		indices = (minResidualDegrees == minimumResidualDegree)\
		                                 .nonzero()[1]
		indices = indices + t
		if indices.shape[0] == 1:
			columnC = indices[0]
		else:
			randomIndex = randint(0,indices.shape[0],(1,1))[0][0]
			columnC = indices[randomIndex]

			#Htemp = H_t.copy()
		Htemp = H_t.copy()

		if minimumResidualDegree == 1:
			# This is the 'extend' case
			rowThatContainsNonZero = H_residual[:,columnC-t]\
			                                   .nonzero()[0][0]
			
			# swap column c with column t (book says t+1 but we 
			# index from 0, not 1)
			Htemp[:,columnC] = H_t[:,t]
			Htemp[:,t] = H_t[:,columnC]
			H_t = Htemp.copy()
			Htemp = H_t.copy()
			# swap row r with row t (book says t+1 but we index from 
			# 0, not 1)
			Htemp[rowThatContainsNonZero + t,:] = H_t[t,:]
			Htemp[t,:] = H_t[rowThatContainsNonZero + t,:]
			H_t = Htemp.copy()
			Htemp = H_t.copy()
		else:
			# This is the 'choose' case
			rowsThatContainNonZeros = H_residual[:,columnC-t]\
			                                    .nonzero()[0]
			
			# swap column c with column t (book says t+1 but we 
			# index from 0, not 1)
			Htemp[:,columnC] = H_t[:,t]
			Htemp[:,t] = H_t[:,columnC]
			H_t = Htemp.copy()
			Htemp = H_t.copy()

			# swap row r1 with row t
			r1 = rowsThatContainNonZeros[0]
			Htemp[r1+t,:] = H_t[t,:]
			Htemp[t,:] = H_t[r1+t,:]
			numRowsLeft = rowsThatContainNonZeros.shape[0]-1
			H_t = Htemp.copy()
			Htemp = H_t.copy()

			# move the other rows that contain nonZero entries to the
			# bottom of the matrix. We can't just swap them, 
			# otherwise we will be pulling up rows that we pushed 
			# down before. So, use a rotation method.
			for index in arange (1,numRowsLeft+1):
				rowInH_residual = rowsThatContainNonZeros[index]
				rowInH_t = rowInH_residual + t - index +1
				m = n-k
				# move the row with the nonzero element to the 
				# bottom; don't update H_t
				Htemp[m-1,:] = H_t[rowInH_t,:] 
				# now rotate the bottom rows up
				sub_index = 1
				while sub_index < (m - rowInH_t):
					Htemp[m-sub_index-1,:] = H_t[m-sub_index,:]
					sub_index = sub_index+1
				H_t = Htemp.copy()
				Htemp = H_t.copy()

			# save temp H as new H_t
			H_t = Htemp.copy()
			Htemp = H_t.copy()
			g = g + (minimumResidualDegree - 1)

		t = t + 1

	return [H_t, g, t]

def invMod2(squareMatrix):
	inverse = inv(squareMatrix)
	temp    = det(squareMatrix)*inverse
	invMod2array = temp % 2
	t = squareMatrix.shape[0]

	if ((dot(squareMatrix,invMod2array) % 2) - eye(t,t)).any():
		if verbose:print 'Error in mod 2 inverse calculation! Det =',
		if verbose: print det(squareMatrix) % 2
		# FIXME is this the most appropriate error to raise?
		raise linalg.linalg.LinAlgError
	else:
		return invMod2array

def swapColumns(a,b,arrayIn):
	arrayOut = arrayIn.copy()
	arrayOut[:,a] = arrayIn[:,b]
	arrayOut[:,b] = arrayIn[:,a]
	return arrayOut

def moveRowToBottom(i,arrayIn):
	arrayOut = arrayIn.copy()
	numRows = arrayOut.shape[0]
	# push the specified row to the bottom
	arrayOut[numRows-1] = arrayIn[i,:]
	# now rotate the bottom rows up
	index = 2
	while (numRows-index) >= i:
		arrayOut[numRows-index,:] = arrayIn[numRows-index+1]
		index = index + 1
	return arrayOut

def getSystematicGmatrix(H):
	tempArray = H.copy()
	numRows = tempArray.shape[0]
	numColumns = tempArray.shape[1] 
	limit = numRows
	rank = 0
	i = 0

	# create an array to save the column permutations
	columnOrder = arange(numColumns).reshape(1,numColumns)

	# create an array to save the row permutations. we just need
	# this to know which dependent rows to delete
	rowOrder = arange(numRows).reshape(numRows,1)

	while i < limit: 
		# Flag indicating that the row contains a non-zero entry
		found  = False
		for j in arange(i, numColumns):
			if tempArray[i, j] == 1: 
				# Encountered a non-zero entry at (i, j)
				found =  True 
				# Increment rank by 1
				rank = rank + 1    
				# make the entry at (i,i) be 1   
				tempArray = swapColumns(j,i,tempArray)  
				# keep track of the column swapping
				columnOrder = swapColumns(j,i,columnOrder)
				break
		if found == True:
			for k in arange(0,numRows): 
				if k == i: continue
				# Checking for 1's 
				if tempArray[k, i] == 1:
					# add row i to row k
					tempArray[k,:] = tempArray[k,:] + tempArray[i,:]
					# Addition is mod2
					tempArray = tempArray.copy() % 2 
					# All the entries above & below (i, i) are now 0 
			i = i + 1
		if found == False:
			# push the row of 0s to the bottom, and move the bottom
			# rows up (sort of a rotation thing)
			tempArray = moveRowToBottom(i,tempArray)
			# decrease limit since we just found a row of 0s
			limit -= 1 
			# keep track of row swapping
			rowOrder = moveRowToBottom(i,rowOrder)

	G =  tempArray[0:i,:]

	# we don't need the dependent rows
	finalRowOrder = rowOrder[0:i]

	if verbose: 
		print 'G.shape:', G.shape
		print 'rank:', rank
		print 'The redundant junk left on the bottom was size:'
		print tempArray[i:numRows,:].shape
		print 'New column order for H:', columnOrder
		print 'New row order for H:\n', rowOrder
		print 'but we can delete these rows from H:\n',
		print rowOrder[i:numRows]
		print 'so final row order is:\n', finalRowOrder
		print finalRowOrder.shape

	# let's reorder H, per the permutations taken above
	# first, put rows in order, omitting the dependent rows
	newNumberOfRowsForH = finalRowOrder.shape[0]
	newH = zeros((newNumberOfRowsForH, numColumns))
	for index in arange(newNumberOfRowsForH):
		newH[index,:] = H[finalRowOrder[index],:]

	# next, put the columns in order
	tempHarray = newH.copy()
	for index in arange(numColumns):
		newH[:,index] = tempHarray[:,columnOrder[0,index]]

	if verbose:
		print 'original H.shape:', H.shape
		print 'newH.shape:', newH.shape
		
	return [G, newH]

def printArrayToFile(arrayName,filename):
	numRows = arrayName.shape[0]
	numColumns = arrayName.shape[1]
	myfile = open(filename,'w')
	for row in arange(numRows):
		for col in arange(numColumns):
			element = arrayName[row,col].astype(int)
			tempstring = `element` + ' '
			myfile.write(tempstring)
		myfile.write('\n')
	myfile.close()