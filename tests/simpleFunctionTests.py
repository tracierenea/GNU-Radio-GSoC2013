#!/usr/bin/python

from numpy import *
from LDPCdecoding import *
from LDPCpreprocessing import *
from LDPChelperFunctions import readAlistFile

####### test the bit flip decoder (aka hard decoding) 

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
transmittedCodeword = array([1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0])
print '\n\n --- Bit flip algorithm test:'

receivedCodeword= array([0,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0])
print '\n  Example part a, two bits fliped:\ncodeword received:'
print receivedCodeword
print 'codeword solution:\n', 
print bitFlipDecoder(maxIterations,H,receivedCodeword)
print 'transmitted codeword:\n',transmittedCodeword


receivedCodeword = array([1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0])
print '\n  Example part b, two bits flipped:\ncodeword received:'
print receivedCodeword
print 'codeword solution:\n', 
print bitFlipDecoder(maxIterations,H,receivedCodeword)
print 'transmitted codeword:\n',transmittedCodeword

receivedCodeword = array([0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,0])
print '\n Example part c, 3 bits flipped:\ncodeword received:'
print receivedCodeword
print 'codeword solution:\n', 
print bitFlipDecoder(maxIterations,H,receivedCodeword)
print 'transmitted codeword:\n',transmittedCodeword
# this received word will cause the function to get caught in a loop 
# of flipping the same bits back and forth. No solution will be 
# reached, even if maxIterations is raised.

####### test the regular LDPC code constructor
print '\n\n --- Code construtor test:'
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

####### test function readAlistFile (read an alist file)
print '\nTest 5 - readAlistFile test'
H = readAlistFile("12.4.3.111.alist",1)

# this true H, and the alist file is from:
# http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html
trueH = array([[0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0],
	           [0,0,0,1,0,0,1,0,1,0,0,0,1,0,0,0],
	           [0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0],
	           [0,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0],
	           [0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,1],
	           [1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0],
	           [0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,0],
	           [0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,1],
	           [1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1],
	           [0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0],
	           [0,1,0,0,0,0,0,0,0,0,1,0,1,0,1,0],
	           [1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0]])

if allclose(H,trueH):
	print 'Successful test: H matrix created matches true H.'
else:
	print 'Test not successful:',
	print 'H matrix created does not match true H.'
