#!/usr/bin/python

from numpy import *
from numpy.random import shuffle, randint
from numpy.linalg import inv, det

verbose = 1

# FIX - I need to go in an specify argument size constraints
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
		# of the unsatisfied parity checks involve that bit. To do this: 

		# First find the nonzero entries in the syndrome. The entry numbers 
		# correspond to rows of interest in H.
		rowsToLookAtInH = array(nonzero(syndrome)) 

		# Second, for each bit, determine how many of unsatisfied parity 
		# checks involve this bit and store this count in an array. 
		counts = zeros_like(receivedCodeword)

		for row in rowsToLookAtInH.transpose():
			for bitNumber in arange(n):
				if H[row,bitNumber] > 0:
					counts[bitNumber] = counts[bitNumber] + 1

		# Next, determine which bit(s) is associated with the most 
		# unsatisfied parity checks, and flip it/them
		print 'counts:', counts
		bitsToFlip = where(counts==counts.max())
		for bitNumber in bitsToFlip[0]:
			if verbose: print 'We need to flip bit:', bitNumber
			testCodeword[bitNumber] = bitwise_xor(testCodeword[bitNumber],1)

		if verbose: print 'New codeword candidate: ',testCodeword
		syndrome = calcSyndrome(H,testCodeword)
		if haveMatch(syndrome):
			if verbose: print 'Codeword declared to be:', testCodeword
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
	# Reference: Turbo Coding for Satellite and Wireless Communications, sec 9.3

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
		range1 = n-(row+1)*q
		range2 = n-row*q
		submatrix1[row,range1:range2] = 1
		# originally I had the 1s going top to bottom, but then you 
		# get stuck with a singular C2 matrix. Not sure yet why all the 
		# books have it that way...

	H = submatrix1

	# create the other submatrices and vertically stack them on
	submatrixNum = 2
	newColumnOrder = arange(n)

	while submatrixNum <= p:
		submatrix = zeros((m/p,n))
		shuffle(newColumnOrder)

		for columnNum in arange(n):
			submatrix[:,columnNum] = submatrix1[:,newColumnOrder[columnNum]]

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

		minimumResidualDegree = minResidualDegrees.min()
		if minimumResidualDegree == 0:
			if any(minResidualDegrees):
				minimumResidualDegree = 1

		# get indices of all of the columns in H_t that have degree
		# equal to the min positive residual degree, then pick at
		# random column c
		indices = (minResidualDegrees == minimumResidualDegree).nonzero()[1]
		if indices.shape[0] == 1:
			columnC = indices[0]
		else:
			randomIndex = randint(0,indices.shape[0],(1,1))[0][0]
			columnC = indices[randomIndex]

		if minimumResidualDegree == 1:
			rowThatContainsNonZero = H_residual[:,columnC].nonzero()[0][0]
			Htemp = H_t.copy()
			
			# swap column c with column t (book says t+1 but we index from 0, not 1)
			Htemp[:,columnC+t] = H_t[:,t]
			Htemp[:,t] = H_t[:,columnC+t]
			H_t = Htemp.copy()
			Htemp = H_t.copy()
			# swap row r with row t (book says t+1 but we index from 0, not 1)
			Htemp[rowThatContainsNonZero + t,:] = H_t[t,:]
			Htemp[t,:] = H_t[rowThatContainsNonZero + t,:]
			H_t = Htemp.copy()
			Htemp = H_t.copy()
		else:
			rowsThatContainNonZeros = H_residual[:,columnC].nonzero()[0]
			Htemp = H_t.copy()
			# swap column c with column t (book says t+1 but we index from 0, not 1)
			Htemp[:,columnC+t] = H_t[:,t]
			Htemp[:,t] = H_t[:,columnC+t]
			H_t = Htemp.copy()
			Htemp = H_t.copy()

			# swap row r1 with row t
			r1 = rowsThatContainNonZeros[0]
			Htemp[r1,:] = H_t[t,:]
			Htemp[t,:] = H_t[r1,:]
			numRowsLeft = rowsThatContainNonZeros.shape[0]-1

			# move the other rows that contain nonZero entries to the
			# bottom of the matrix
			for index in arange (1,numRowsLeft+1):
				rowInH_residual = rowsThatContainNonZeros[index]
				
				if (k-(rowInH_residual+t)) <= numRowsLeft:
					# this row is already at the bottom
					continue
				else: 

					tempIndex = index
					possibleBottomIndex = k-tempIndex
					while Htemp[possibleBottomIndex,t] == 1:
						tempIndex = tempIndex+1
						possibleBottomIndex = k-tempIndex
					bottomRow = possibleBottomIndex			
					Htemp[bottomRow,:] = H_t[rowInH_residual+t,:]
					Htemp[rowInH_residual+t,:] = H_t[bottomRow,:]
				H_t = Htemp.copy()
				Htemp = H_t.copy()

			# save temp H as new H_t
			H_t = Htemp.copy()
			Htemp = H_t.copy()
			g = g + (minimumResidualDegree - 1)

		t = t + 1

	return H_t
	# FIX will also need to return gap (g) and parameter t