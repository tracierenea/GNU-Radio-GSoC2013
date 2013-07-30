#!/usr/bin/python
from numpy import *
from LDPCencoding import RichardsonUrbankeEncoder
from LDPCdecoding import *
from LDPCpreprocessing import *
from numpy.random import random_integers
import sys

n = int(sys.argv[1])
p = int(sys.argv[2])
q = int(sys.argv[3])
print 'Running test for (n,p,q) = (',n,',',p,',',q,')'

H = regularLDPCcodeConstructor(n,p,q)
newH = getFullRankHmatrix(H)
[invT,invPhi,E,A,B,D,bestH,n,k,g]=getParametersForEncoding(newH,100)

print '\n\nTest of encoding k =',k,
print 'randomly generated information bits'
print 'H.shape:', bestH.shape
print 'gap:',g
print 'Codeword size:', bestH.shape[1]
rate = divide((k*1.0),bestH.shape[1])
print 'Rate:','%2.3F' % rate
maxIterationsLimit = 100

passCount = failCount = 0
for testNum in arange(10):
	print '\nTest number:', testNum

	s = random_integers(0,1,k).reshape(k,1)
	x = RichardsonUrbankeEncoder(invT,invPhi,E,A,B,D,bestH,n,k,g,s)
	numBits = x.shape[0]
	transmittedCodeword = x.copy()
	
	# verify:
	testArray = dot(bestH,x) % 2
	
	# so many rounding errors! I'm seeing 2.0000000 in testArray 
	# after doing % 2
	for index in arange(testArray.shape[0]):
		if (abs(2-testArray[index,0])) < 0.01:
			# this is a 2, and 2%2 = 0
			testArray[index,0] = 0
		elif (abs(1-testArray[index,0])) < 0.01:
			# this is a 1
			testArray[index,0] = 1
		elif (abs(0-testArray[index,0])) < 0.01:
			# this is a 0
			testArray[index,0] = 0
		else: 
			print 'What value is this???', testArray[index,0]

	if testArray.any():
		failCount += 1
	else: 
		# if codeword passed parity check, perform bit flip decoding
		# tests
		passCount += 1
		
		totalCount = 0
		correctlyDecoded = 0
		notCorrectlyDecoded = 0
		for index in arange(numBits):
			totalCount += 1
			receivedCodeword = transmittedCodeword.copy()
			receivedCodeword[index] =\
			                logical_xor(transmittedCodeword[index],1)
			correctedCodeword =\
		bitFlipDecoder(maxIterationsLimit, bestH,receivedCodeword)
			if allclose(correctedCodeword,transmittedCodeword):
				correctlyDecoded += 1
			else:
				notCorrectlyDecoded += 1
		successRate = divide(correctlyDecoded, totalCount) * 100
		print 'Tested flipping every bit, total:',totalCount,
		print '\nCorrectly decoded:',correctlyDecoded,
		print '/ Not correctly decoded:',notCorrectlyDecoded,
		print '/ Success rate:','%2.3F'% successRate

		totalCount = 0
		correctlyDecoded = 0
		notCorrectlyDecoded = 0
		# I think 10000 is a sufficent number of tests...
		while totalCount < 10000: 
			index1 = randint(0,numBits-1)
			index2 = randint(0,numBits-1)
			while index2 == index1:
				index2 = randint(0,numBits-1)
			totalCount += 1
			receivedCodeword = transmittedCodeword.copy()
			receivedCodeword[index1] =\
				           logical_xor(transmittedCodeword[index1],1)
			receivedCodeword[index2] =\
 				          logical_xor(transmittedCodeword[index2],1)
			correctedCodeword = \
			bitFlipDecoder(maxIterationsLimit,bestH,receivedCodeword)
			if allclose(correctedCodeword,transmittedCodeword):
				correctlyDecoded += 1
			else:
				notCorrectlyDecoded += 1
			
			successRate=divide(correctlyDecoded*1.0, totalCount)*100
		print 'Tests ran with 2 bits flipped:',totalCount,
		print '\nCorrectly decoded:',correctlyDecoded,
		print '/ Not correctly decoded:',notCorrectlyDecoded,
		print '/ Success rate:','%2.3F'% successRate

print passCount, 'encoding parity check tests passed and', 
print failCount, 'tests failed.'
