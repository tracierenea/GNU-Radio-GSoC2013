#!/usr/bin/python

from numpy import *

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