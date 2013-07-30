#!/usr/bin/python
from LDPCpreprocessing import *
from LDPCencoding import RichardsonUrbankeEncoder
from numpy import *
from numpy.random import random_integers

############ these are the preprocessing steps  #####################

H = regularLDPCcodeConstructor(1200,3,6)
newH = getFullRankHmatrix(H)
[invT,invPhi,E,A,B,D,bestH,n,k,g] =getParametersForEncoding(newH,100)

############ this is all real-time encoding #########################

print '\n\nTest of encoding k=',k,'random information bits'
print 'H.shape:', bestH.shape
print 'gap:',g
print 'Codeword size:', bestH.shape[1]
rate = divide((k*1.0),bestH.shape[1])
print 'Rate:','%2.3F' % rate

passCount = failCount = 0
for testNum in arange(500):

	s = random_integers(0,1,k).reshape(k,1)
	x = RichardsonUrbankeEncoder(invT,invPhi,E,A,B,D,bestH,n,k,g,s)
	
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
		passCount += 1

print passCount, 'tests passed and', failCount, 'tests failed.'
