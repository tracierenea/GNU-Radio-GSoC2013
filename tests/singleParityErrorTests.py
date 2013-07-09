#!/usr/bin/python

from numpy import *
from LDPCfunctions import *

H = array([[1, 1, 0, 1, 0, 0],
	       [0, 1, 1, 0, 1, 0],
	       [1, 0, 1, 0, 0, 1]])

transmittedWord = array([[1], [1], [0], [0], [1], [1]])
print 'transmitted word:\n', transmittedWord

[numRows, numColumns] = transmittedWord.shape
for rowNum in arange(numRows):
	print '\ntest #', rowNum

	# create one parity symbol error
	receivedWord = transmittedWord.copy()
	if   receivedWord[rowNum] == 0: receivedWord[rowNum] = 1
	elif receivedWord[rowNum] == 1: receivedWord[rowNum] = 0
	print 'Received word with one parity error:\n', receivedWord

	# try to get a correct codeword
	testCodeword = singleParityErrorFix(H,receivedWord)
	print 'new candidate codeword:\n', testCodeword

	# if we got something new, test it
	if not array_equiv(testCodeword,receivedWord):
		syndrome = calcSyndrome(H,testCodeword)
		if haveMatch(syndrome):
			print 'Received word had a single parity error. Corrected.'
			if array_equiv(testCodeword,transmittedWord):
				print 'Corrected codeword matches transmitted word.'
			else:
				print '\'Corrected codeword does not match transmitted word!'
		else:
			print '  - Single parity correction not succesful.'
	# if we didn't get something new, something's not working
	else:
		print 'Didn\'t get a new word from singleParityErrorFix function'


	# so this little test shows that for this little symbol word, 
	# an error on any one of the symbols can be corrected 
