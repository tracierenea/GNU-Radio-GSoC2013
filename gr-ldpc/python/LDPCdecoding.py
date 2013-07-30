#!/usr/bin/python

from numpy import *
from numpy.random import shuffle, randint
from numpy.linalg import inv, det
from LDPCpreprocessing import calcSyndrome, haveMatch

# 0 gives no debug output, 1 gives a little, 2 gives a lot
verbose = 0

def bitFlipDecoder(maxIterations, H, codeword):
	receivedCodeword = codeword.copy()
	sizeOfReceivedCodeword = receivedCodeword.shape
	n = sizeOfReceivedCodeword[0] # number of symbols/bit in codeword

	syndrome = calcSyndrome(H,receivedCodeword)
	if haveMatch(syndrome):
		if verbose: 
			print 'bitFlipDecoder: Received valid codeword. No bit',
			print ' flips required.'
		return receivedCodeword
	else:
		if verbose > 1: 
			print 'Evaluating codeword:\n', receivedCodeword
		testCodeword = receivedCodeword
		iteration = 1;
	
	while iteration <= maxIterations:
		if verbose >1: print 'bitFlipDecoder: iteration:', iteration

		# For each of the n bits in the codeword, determine how many
		# of the unsatisfied parity checks involve that bit. To do 
		# this: 
		# First find the nonzero entries in the syndrome. The entry 
		# numbers correspond to rows of interest in H.
		rowsToLookAtInH = array(nonzero(syndrome))[0] 

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
		if verbose > 1: print 'counts:\n', counts
		bitsToFlip = where(counts==counts.max())
		for bitNumber in bitsToFlip[0]:
			if verbose > 1: print 'We need to flip bit:', bitNumber
			testCodeword[bitNumber] = \
			             bitwise_xor(int(testCodeword[bitNumber]),1)

		if verbose > 1: 
			print 'New codeword candidate:\n',testCodeword
		syndrome = calcSyndrome(H,testCodeword)
		if haveMatch(syndrome):
			if verbose: print 'Valid codeword found.'
			if verbose>1: 
				print 'Codeword declared to be:\n', testCodeword
			return testCodeword
		else:
			iteration = iteration + 1

	else:
		if verbose: 
			print 'Max iteration count of', maxIterations,
			print 'has been reached without finding valid codeword.',
			print 'Returning received codeword.'
		return receivedCodeword
