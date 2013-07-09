#!/usr/bin/python

from numpy import *
from LDPCfunctions import *

#####################################################################

####### test the encoder
print ' --- Encoder test:'
# In the future we will build a H matrix in another module, but for 
# now, start with some assumed small parity check matrix H
H = array([[1, 1, 0, 1, 0, 0],
	       [0, 1, 1, 0, 1, 0],
	       [1, 0, 1, 0, 0, 1]])
A = array([[1], [1], [0]]) # column vector of k information symbols
codeword = matrixMultiplierEncoder(H,A)
print codeword


####### test the single parity error decoder
print '\n --- Single parity error decode test:'
print 'H matrix: \n', H
print '\nCASE 1:\n'

transmittedWord = array([[1], [1], [0], [0], [1], [1]])
print 'transmitted word:\n', transmittedWord

# verify no mods to a correctly received codeword  ######
receivedWord = array([[1], [1], [0], [0], [1], [1]])
print 'received word  (no error case):\n', receivedWord

testCodeword = singleParityErrorFix(H,receivedWord)

if not array_equiv(testCodeword,receivedWord):
	print 'new candidate codeword:\n', testCodeword
	syndrome = calcSyndrome(H,testCodeword)
	if haveMatch(syndrome):
		print '  - Received word had a single parity error. Corrected.'
	else:
		print '  - Single parity correction not succesful.'

# verify case of one symbol error  ######
receivedWord = array([[1], [1], [1], [0], [1], [1]])
print 'received word  (1 symbol error):\n', receivedWord

testCodeword = singleParityErrorFix(H,receivedWord)

if not array_equiv(testCodeword,receivedWord):
	print 'new candidate codeword:\n', testCodeword
	syndrome = calcSyndrome(H,testCodeword)
	if haveMatch(syndrome):
		print '  - Received word had a single parity error. Corrected.'
	else:
		print '  - Single parity correction not succesful.'


# verify that this doesn't work with 2 symbol errors
print '\nCASE 2:\n'
transmittedWord = array([[1], [1], [0], [0], [1], [1]])
print 'transmitted word:\n', transmittedWord

receivedWord = array([[1], [1], [1], [0], [0], [1]])
print 'received word  (2 symbol errors):\n', receivedWord

testCodeword = singleParityErrorFix(H,receivedWord)
print 'new candidate codeword:\n', testCodeword

if array_equiv(testCodeword,receivedWord):
	syndrome = calcSyndrome(H,testCodeword)
	if haveMatch(syndrome):
		print '  - Received word had a single parity error. Corrected.'
		print '  - Corrected codeword is:\n', testCodeword
	else:
		print '  - Single parity correction not succesful.'
print 'The corrected codeword does not match transmitted word.'

# the codeword found was 111000, not correct

####### test the bit flip decoder (aka hard decoding) proposed by Gallager

# From Fundamentals of Error Correcting Codes, Example 15.6.1
# This is a parity check matrix for a (16,3,4) LDPC code, so 
# note that there are 4 1's in every row and 3 1's in every column.
# We will get this from a generator function later, but for now, just
# use this example from the text

H = array([[1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	       [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
	       [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
	       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
	       [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0],
	       [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0],
	       [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
	       [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1],
	       [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
	       [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1],
	       [0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0],
	       [0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]])

maxIterations = 20
print '\n\nBit flip algorithm test:'

receivedCodeword = array([0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0])
print '\nExample part a:'
bitFlipDecoder(maxIterations,H,receivedCodeword)
# the solution should be: 1100101001100000. Two errors corrected.

receivedCodeword = array([1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
print '\nExample part b:'
bitFlipDecoder(maxIterations,H,receivedCodeword)
# the solution should be: 1100101001100000. Two errors corrected

receivedCodeword = array([0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0])
print '\nExample part c:'
bitFlipDecoder(maxIterations,H,receivedCodeword)
# this received word will cause the function to get caught in a loop of 
# flipping the same bits back and forth. No solution will be reached, even
# if maxIterations is raised.

####### test the regular LDPC code constructor
print 'Test 1 - (20,3,4) code:'
H = regularLDPCcodeConstructor(20,3,4)
print H

print '\nTest 2 - (16,4,3) code:'
H = regularLDPCcodeConstructor(16,3,4)
print H

# This one should return an error
print '\nTest 3 - (16,3,3) code:'
H = regularLDPCcodeConstructor(16,3,3)
print H

print '\nTest 4 - (1000,10,20) code:'
H = regularLDPCcodeConstructor(1000,10,20)
print H
print H.shape
