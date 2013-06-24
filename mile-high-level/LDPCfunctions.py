#!/usr/bin/python

from numpy import *

verbose = 1

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

def bilFlipDecoder(maxIterations, H, codeword):
	# deleting stuff that didn't work. in progress...
	return codeword