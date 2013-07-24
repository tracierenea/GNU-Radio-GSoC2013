#!/usr/bin/python

from numpy import *
from numpy.random import shuffle, randint
from numpy.linalg import inv, det
from itertools import permutations

# 0 gives no debug output, 1 gives a little, 2 gives a lot
verbose = 0

def calcSyndrome(H,codeword):
	syndrome = dot(H,codeword) % 2	# use modulo 2 operations
	return syndrome

def haveMatch(syndrome):
	check = 1
	if any(syndrome):
		# the matrix is not all zeros, it's not a codeword
		check = 0
	return check

def regularLDPCcodeConstructor(n,p,q):
	# n = codeword length
	# p = column weight
	# q = row weight

	# Following Gallager's approach where we create p submatrices. 
	# Reference: Turbo Coding for Satellite and Wireless 
	# Communications, sec 9.3

	# for this algorithm, n/p must be an integer, because we need the
	# number of rows in eacn submatrix to be a whole number
	ratioTest = (n*1.0)/q
	if ratioTest%1 != 0:
		print 'The ratio of inputs n/q must be a whole number.'
		return

	# FIX: There should probably be other guidelines for n/p/q,
	# but I have not found any specifics in the literature....

	# First one first: 
	m = (n*p)/q  # number of rows in H matrix
	submatrix1 = zeros((m/p,n))  

	for row in arange(m/p):

		# these lines have the ones going from bottom to top, left
		# to right across the matrix. Switching to top to bottom, 
		# left to right, because I think it will be faster to then
		# convert it to systematic form
		# range1 = n-(row+1)*q
		# range2 = n-row*q

		range1 = row*q
		range2 = (row+1)*q 
		submatrix1[row,range1:range2] = 1

	H = submatrix1

	# create the other submatrices and vertically stack them on
	submatrixNum = 2
	newColumnOrder = arange(n)

	while submatrixNum <= p:
		submatrix = zeros((m/p,n))
		shuffle(newColumnOrder)

		for columnNum in arange(n):
			submatrix[:,columnNum] = \
			                  submatrix1[:,newColumnOrder[columnNum]]

		H = vstack((H,submatrix))

		submatrixNum = submatrixNum + 1 

	# double check the row weight and column weights
	size = H.shape
	rows = size[0]
	cols = size[1]

	# check the row weights
	for rowNum in arange(rows):
		nonzeros = array(H[rowNum,:].nonzero())
		if nonzeros.shape[1] != q:
			print 'Row', rowNum, 'has incorrect weight!'
			return

	# check the column weights
	for columnNum in arange(cols):
		nonzeros = array(H[:,columnNum].nonzero())
		if nonzeros.shape[1] != p:
			print 'Row', columnNum, 'has incorrect weight!'
			return

	# At this point, we can use this matrix for parity checks, 
	# but it is not in systematic form so it can't be used
	# to find the parity matrix P, and therefore can't be used
	# for encoding via the generation matrix. 

	# Putting H into its systematic form requires row operations. 
	# I cannot figure out how to write this into an algorithm for
	# matrices of unknown size, and also using mod 2 operations.
	# <in work....>

	# There are some published steps that require taking the 
	# determinant but taking the determinant of a huge matrix
	# won't work.....

	return H

def greedyUpperTriangulation(H):
	# This function performs row/column permutations to bring
	# H into approximate upper triangular form via greedy 
	# upper triangulation method outlined in Modern Coding 
	# Theory Appendix 1, Section A.2
	H_t  = H.copy()

	# Per email from Dr. Urbanke, author of this textbook, this
	# algorithm requires H to be full rank
	if linalg.matrix_rank(H_t) != H_t.shape[0]:
			print 'Rank of H:',linalg.matrix_rank(tempArray)
			print 'H has', H_t.shape[0], 'rows'
			print 'Error: H must be full rank.'
			return

	size = H_t.shape
	n = size[1]
	k = n - size[0]
	g = t = 0

	while t != (n-k-g):
		H_residual = H_t[t:n-k-g,t:n]
		size       = H_residual.shape
		numRows    = size[0]
		numCols    = size[1]

		minResidualDegrees = zeros((1,numCols))

		for colNum in arange(numCols):
			nonZeroElements = array(H_residual[:,colNum].nonzero())
			minResidualDegrees[0,colNum] = nonZeroElements.shape[1]

		# find the minimum nonzero residual degree
		nonZeroElementIndices = minResidualDegrees.nonzero()
		nonZeroElements=minResidualDegrees[nonZeroElementIndices[0],\
		                                   nonZeroElementIndices[1]]
		minimumResidualDegree = nonZeroElements.min()

		# get indices of all of the columns in H_t that have degree
		# equal to the min positive residual degree, then pick at
		# random column c
		indices = (minResidualDegrees == minimumResidualDegree)\
		                                 .nonzero()[1]
		indices = indices + t
		if indices.shape[0] == 1:
			columnC = indices[0]
		else:
			randomIndex = randint(0,indices.shape[0],(1,1))[0][0]
			columnC = indices[randomIndex]

			#Htemp = H_t.copy()
		Htemp = H_t.copy()

		if minimumResidualDegree == 1:
			# This is the 'extend' case
			rowThatContainsNonZero = H_residual[:,columnC-t]\
			                                   .nonzero()[0][0]
			
			# swap column c with column t (book says t+1 but we 
			# index from 0, not 1)
			Htemp[:,columnC] = H_t[:,t]
			Htemp[:,t] = H_t[:,columnC]
			H_t = Htemp.copy()
			Htemp = H_t.copy()
			# swap row r with row t (book says t+1 but we index from 
			# 0, not 1)
			Htemp[rowThatContainsNonZero + t,:] = H_t[t,:]
			Htemp[t,:] = H_t[rowThatContainsNonZero + t,:]
			H_t = Htemp.copy()
			Htemp = H_t.copy()
		else:
			# This is the 'choose' case
			rowsThatContainNonZeros = H_residual[:,columnC-t]\
			                                    .nonzero()[0]
			
			# swap column c with column t (book says t+1 but we 
			# index from 0, not 1)
			Htemp[:,columnC] = H_t[:,t]
			Htemp[:,t] = H_t[:,columnC]
			H_t = Htemp.copy()
			Htemp = H_t.copy()

			# swap row r1 with row t
			r1 = rowsThatContainNonZeros[0]
			Htemp[r1+t,:] = H_t[t,:]
			Htemp[t,:] = H_t[r1+t,:]
			numRowsLeft = rowsThatContainNonZeros.shape[0]-1
			H_t = Htemp.copy()
			Htemp = H_t.copy()

			# move the other rows that contain nonZero entries to the
			# bottom of the matrix. We can't just swap them, 
			# otherwise we will be pulling up rows that we pushed 
			# down before. So, use a rotation method.
			for index in arange (1,numRowsLeft+1):
				rowInH_residual = rowsThatContainNonZeros[index]
				rowInH_t = rowInH_residual + t - index +1
				m = n-k
				# move the row with the nonzero element to the 
				# bottom; don't update H_t
				Htemp[m-1,:] = H_t[rowInH_t,:] 
				# now rotate the bottom rows up
				sub_index = 1
				while sub_index < (m - rowInH_t):
					Htemp[m-sub_index-1,:] = H_t[m-sub_index,:]
					sub_index = sub_index+1
				H_t = Htemp.copy()
				Htemp = H_t.copy()

			# save temp H as new H_t
			H_t = Htemp.copy()
			Htemp = H_t.copy()
			g = g + (minimumResidualDegree - 1)

		t = t + 1

	if g == 0:
		if verbose: print 'Error: gap is 0.'
		return

	# we need to ensure phi is nonsingular
	T = H_t[0:t, 0:t]
	E = H_t[t:t+g,0:t]
	A = H_t[0:t,t:t+g]
	C = H_t[t:t+g,t:t+g]
	D = H_t[t:t+g,t+g:n] 

	invTmod2array = invMod2(T)
	temp1  = dot(E,invTmod2array) % 2
	temp2  = dot(temp1,A) % 2
	phi    = (C - temp2) % 2
	if phi.any():
		try:
			# try to take the inverse of phi
			invPhi = invMod2(phi)
		except linalg.linalg.LinAlgError:
			# phi is singular
			if verbose > 1: print 'Initial phi is singular'
		else:
			# phi is nonsingular, so we need to use this version of H
			if verbose > 1: print 'Initial phi is nonsingular'
			return [H_t, g, t]
	else:
		if verbose: print 'Initial phi is empty or all zeros:\n', phi

	# if the C and D submatrices are all zeros, there is no point in
	# shuffling them around in an attempt to find a good phi
	if not (C.any() or D.any()):
		if verbose: 
			print 'C and D are all zeros. There is no hope in',
			print ' finding a nonsingular phi matrix. '
		return

	# we can't look at every row/column pertubation possibility
	# because there would be (n-t)! column shuffles and g! row
	# shuffles. g has gotten up to 12 in tests, so 12! would still
	# take quite some time. Instead, we will just pick an arbitrary 
	# number of max iterations to perform, then break
	maxIterations = 300
	iterationCount = 0
	columnsToShuffle = arange(t,n) 
	rowsToShuffle = arange(t,t+g)

	while iterationCount < maxIterations:
		if verbose > 1: print 'iterationCount:', iterationCount
		tempH = H_t.copy()

		shuffle(columnsToShuffle)
		shuffle(rowsToShuffle)
		index = 0
		for newDestinationColumnNumber in arange(t,n):
			oldColumnNumber = columnsToShuffle[index]
			tempH[:,newDestinationColumnNumber] = \
			                                   H_t[:,oldColumnNumber]
			index +=1

		tempH2 = tempH.copy()
		index = 0
		for newDesinationRowNumber in arange(t,t+g):
			oldRowNumber = rowsToShuffle[index]
			tempH[newDesinationRowNumber,:] = tempH2[oldRowNumber,:]
			index +=1
		
		# now test this new H matrix
		H_t = tempH.copy()
		T = H_t[0:t, 0:t]
		E = H_t[t:t+g,0:t]
		A = H_t[0:t,t:t+g]
		C = H_t[t:t+g,t:t+g]
		invTmod2array = invMod2(T)
		temp1  = dot(E,invTmod2array) % 2
		temp2  = dot(temp1,A) % 2
		phi    = (C - temp2) % 2
		if phi.any():
			try:
				# try to take the inverse of phi
				invPhi = invMod2(phi)
			except linalg.linalg.LinAlgError:
				# phi is singular
				if verbose > 1: print 'Phi is still singular'
			else:
				# phi is nonsingular, so we're done
				if verbose: 
					print 'Found a nonsingular phi'
					print 'iterationCount = ', iterationCount
				return [H_t, g, t]
		else:
			if verbose > 1: print 'phi is all zeros'

		iterationCount +=1

	# if we've reached this point, then we haven't found a
	# version of H that has a nonsingular phi
	if verbose: print '--- Error: nonsingular phi matrix not found.'


def invMod2(squareMatrix):
	A = squareMatrix.copy()

	# special case for one element array [1]
	if A.size == 1 and A[0] == 1:
		return array([1])

	Ainverse = inv(A)
	B = dot(det(A),Ainverse)
	C = B % 2
	t = A.shape[0]

	if ((dot(A,C) % 2) - eye(t,t)).any():
		if verbose:print 'Error in mod 2 inverse calculation!'
		# FIXME is this the most appropriate error to raise?
		raise linalg.linalg.LinAlgError
	else:
		return C

def swapColumns(a,b,arrayIn):
	arrayOut = arrayIn.copy()
	arrayOut[:,a] = arrayIn[:,b]
	arrayOut[:,b] = arrayIn[:,a]
	return arrayOut

def moveRowToBottom(i,arrayIn):
	arrayOut = arrayIn.copy()
	numRows = arrayOut.shape[0]
	# push the specified row to the bottom
	arrayOut[numRows-1] = arrayIn[i,:]
	# now rotate the bottom rows up
	index = 2
	while (numRows-index) >= i:
		arrayOut[numRows-index,:] = arrayIn[numRows-index+1]
		index = index + 1
	return arrayOut

def getSystematicGmatrix(H):
	# This function finds the systematic form of the generator
	# matrix G. The form is G = [I P] where I is an identity matrix
	# and P is the parity submatrix. If the H matrix provided 
	# is not full rank, then dependent rows will be deleted.
	tempArray = H.copy()
	numRows = tempArray.shape[0]
	numColumns = tempArray.shape[1] 
	limit = numRows
	rank = 0
	i = 0

	while i < limit: 
		# Flag indicating that the row contains a non-zero entry
		found  = False
		for j in arange(i, numColumns):
			if tempArray[i, j] == 1: 
				# Encountered a non-zero entry at (i, j)
				found =  True 
				# Increment rank by 1
				rank = rank + 1    
				# make the entry at (i,i) be 1   
				tempArray = swapColumns(j,i,tempArray)  
				break
		if found == True:
			for k in arange(0,numRows): 
				if k == i: continue
				# Checking for 1's 
				if tempArray[k, i] == 1:
					# add row i to row k
					tempArray[k,:] = tempArray[k,:] + tempArray[i,:]
					# Addition is mod2
					tempArray = tempArray.copy() % 2 
					# All the entries above & below (i, i) are now 0 
			i = i + 1
		if found == False:
			# push the row of 0s to the bottom, and move the bottom
			# rows up (sort of a rotation thing)
			tempArray = moveRowToBottom(i,tempArray)
			# decrease limit since we just found a row of 0s
			limit -= 1 

	# the rows below i are the dependent rows, which we discard
	G = tempArray[0:i,:]

	return G

def getFullRankHmatrix(H):
	# This function accepts a parity check matrix H and, if it is not
	# already full rank, will determine which rows are dependent and 
	# remove them. This updated matrix will be returned.
	tempArray = H.copy()
	if linalg.matrix_rank(tempArray) == tempArray.shape[0]:
		if verbose > 1: print 'Returning H; it is already full rank.'
		return tempArray
	
	numRows = tempArray.shape[0]
	numColumns = tempArray.shape[1] 
	limit = numRows
	rank = 0
	i = 0

	# create an array to save the column permutations
	columnOrder = arange(numColumns).reshape(1,numColumns)

	# create an array to save the row permutations. we just need
	# this to know which dependent rows to delete
	rowOrder = arange(numRows).reshape(numRows,1)

	while i < limit: 
		if verbose > 1: print 'In getFullRankHmatrix; i:', i
		# Flag indicating that the row contains a non-zero entry
		found  = False
		for j in arange(i, numColumns):
			if tempArray[i, j] == 1: 
				# Encountered a non-zero entry at (i, j)
				found =  True 
				# Increment rank by 1
				rank = rank + 1    
				# make the entry at (i,i) be 1   
				tempArray = swapColumns(j,i,tempArray)  
				# keep track of the column swapping
				columnOrder = swapColumns(j,i,columnOrder)
				break
		if found == True:
			for k in arange(0,numRows): 
				if k == i: continue
				# Checking for 1's 
				if tempArray[k, i] == 1:
					# add row i to row k
					tempArray[k,:] = tempArray[k,:] + tempArray[i,:]
					# Addition is mod2
					tempArray = tempArray.copy() % 2 
					# All the entries above & below (i, i) are now 0 
			i = i + 1
		if found == False:
			# push the row of 0s to the bottom, and move the bottom
			# rows up (sort of a rotation thing)
			tempArray = moveRowToBottom(i,tempArray)
			# decrease limit since we just found a row of 0s
			limit -= 1 
			# keep track of row swapping
			rowOrder = moveRowToBottom(i,rowOrder)

	# we don't need the dependent rows
	finalRowOrder = rowOrder[0:i]

	# let's reorder H, per the permutations taken above
	# first, put rows in order, omitting the dependent rows
	newNumberOfRowsForH = finalRowOrder.shape[0]
	newH = zeros((newNumberOfRowsForH, numColumns))
	for index in arange(newNumberOfRowsForH):
		newH[index,:] = H[finalRowOrder[index],:]

	# next, put the columns in order
	tempHarray = newH.copy()
	for index in arange(numColumns):
		newH[:,index] = tempHarray[:,columnOrder[0,index]]

	if verbose:
		print 'original H.shape:', H.shape
		print 'newH.shape:', newH.shape

	return newH

def printArrayToFile(arrayName,filename):
	# this function is handy because this nunpy command is not
	# working for me: set_printoptions(threshold='nan')
	numRows = arrayName.shape[0]
	numColumns = arrayName.shape[1]
	myfile = open(filename,'w')
	for row in arange(numRows):
		for col in arange(numColumns):
			element = arrayName[row,col].astype(int)
			tempstring = `element` + ' '
			myfile.write(tempstring)
		myfile.write('\n')
	myfile.close()
