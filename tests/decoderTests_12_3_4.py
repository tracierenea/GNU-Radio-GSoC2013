#!/usr/bin/python
from LDPCfunctions import *
from numpy import *
from numpy.random import random_integers

# this (12,3,4) H is from Example A.12 in Modern Coding Theory.
# It's full rank
H = array([[0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0],  # 1
	       [1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0],  # 2
	       [0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1],  # 3
	       [0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0],  # 4
	       [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1],  # 5
	       [1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1]]) # 6

# this returns a regular H but it is not in systematic form
[G, newH] = getSystematicGmatrix(H)

# get the systematic form of H
k = G.shape[0]
n = G.shape[1]
m = n - k
P = G[0:k,k:n]
I_m = identity(m,int)
newH = hstack((P.transpose(),I_m))
s = random_integers(0,1,k).reshape(k,1)
x = dot(G.transpose(),s) % 2		# codeword

myfile = open('DecodingTestOutput_12_3_4.txt','w')
maxIterationsLimit = 20
numBits = x.shape[0]
tempstring = 'Testing codeword length: ' + `numBits`
myfile.write(tempstring)
myfile.write('\n\nFlipping bits\tDecoded Correctly?\n')

# Flipping just one bit first
totalCount = 0
correctlyDecoded = 0
notCorrectlyDecoded = 0
transmittedCodeword = x.copy().astype(int)
for index in arange(numBits):
	totalCount = totalCount + 1
	tempstring = `index`
	myfile.write(tempstring)
	receivedCodeword = transmittedCodeword.copy()
	receivedCodeword[index] =\
	             logical_xor(transmittedCodeword[index],1)
	correctedCodeword = \
	    bitFlipDecoder(maxIterationsLimit, newH, receivedCodeword)
	if allclose(correctedCodeword,transmittedCodeword):
		myfile.write('\t\tyes\t')
		correctlyDecoded = correctlyDecoded + 1
	else:
		myfile.write('\t\tNO\t')
		notCorrectlyDecoded = notCorrectlyDecoded + 1
	tempstring = `transmittedCodeword.transpose()`
	myfile.write(tempstring)
	myfile.write('<-- transmitted\n\t\t\t')
	tempstring = `receivedCodeword.transpose()`
	myfile.write(tempstring)
	myfile.write('<-- received\n\t\t\t')
	tempstring = `correctedCodeword.transpose()`
	myfile.write(tempstring)
	myfile.write('<-- corrected codeword\n\n')
tempstring = 'Tests ran with 1 bit flipped: ' + `totalCount` + '\nCorrectly decoded: ' + `correctlyDecoded` + ', Not correctly decoded: ' + `notCorrectlyDecoded` + '\n\n'
myfile.write(tempstring)

# flipping 2 bits
myfile.write('\n\nFlipping bits\tDecoded Correctly?\n')
totalCount = 0
correctlyDecoded = 0
notCorrectlyDecoded = 0
for index1 in arange(numBits):
	for index2 in arange(index1+1,numBits):
		if index1 == index2: continue
		totalCount = totalCount + 1
		tempstring = `index1` + '  ' + `index2`
		myfile.write(tempstring)
		receivedCodeword = transmittedCodeword.copy()
		receivedCodeword[index1] =\
	             logical_xor(transmittedCodeword[index1],1)
		receivedCodeword[index2] =\
	             logical_xor(transmittedCodeword[index2],1)          
		correctedCodeword = \
	    bitFlipDecoder(maxIterationsLimit, newH, receivedCodeword)
		if allclose(correctedCodeword,transmittedCodeword):
			myfile.write('\t\tyes\t')
			correctlyDecoded = correctlyDecoded + 1
		else:
			myfile.write('\t\tNO\t')
			notCorrectlyDecoded = notCorrectlyDecoded + 1
		tempstring = `transmittedCodeword.transpose()`
		myfile.write(tempstring)
		myfile.write('<-- transmitted\n\t\t\t')
		tempstring = `receivedCodeword.transpose()`
		myfile.write(tempstring)
		myfile.write('<-- received\n\t\t\t')
		tempstring = `correctedCodeword.transpose()`
		myfile.write(tempstring)
		myfile.write('<-- corrected codeword\n\n')
tempstring = 'Tests ran with 2 bits flipped: ' + `totalCount` + '\nCorrectly decoded: ' + `correctlyDecoded` + ', Not correctly decoded: ' + `notCorrectlyDecoded`
myfile.write(tempstring)

myfile.close()