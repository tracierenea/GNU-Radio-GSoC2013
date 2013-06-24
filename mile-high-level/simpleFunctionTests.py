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



















