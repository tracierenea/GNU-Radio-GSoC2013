#!/usr/bin/python
from LDPCfunctions import *
from numpy import *
from numpy.random import random_integers

# Example of how to do encoding per Appendix A in Modern Coding
# Theory (Richardson/Urbanke). 

# The H matrices created by regularLDPCcodeconstructor are
# not full rank. Using this simple textbook example just to
# show that it works if the matrix is full rank.

H = array([[0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0],  # 1
	       [1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0],  # 2
	       [0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1],  # 3
	       [0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0],  # 4
	       [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1],  # 5
	       [1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1]]) # 6
print 'Origiinal H.shape:', H.shape
printArrayToFile(H,'H_12_3_6_original.txt')
[G, newH] = getSystematicGmatrix(H)
printArrayToFile(G,'G_from_12_3_6_H_matrix.txt')
printArrayToFile(newH,'newH_12_3_6.txt')
print 'linalg.matrix_rank(newH):', linalg.matrix_rank(newH)

############ this is all preprocessing #############################

numIterations = 10
print 'Running permutation algorithm', numIterations,
print 'times to ensure the lowest \ngap is found and the',
print 'resulting phi is nonsingular.'

# set this arbitrarily high to force whole loop on first run
g = 10**10
flagFoundBetterH = 0

for index in arange(numIterations):
 	print '============== Index:', index
	[betterH, gap, t]  = greedyUpperTriangulation(newH)
	print 'gap:', gap, ' g:', g

	if gap < g:
		T = betterH[0:t, 0:t]
		E = betterH[t:t+gap,0:t]
		A = betterH[0:t,t:t+gap]
		C = betterH[t:t+gap,t:t+gap]
		invTmod2array = invMod2(T)
		temp1  = dot(E,invTmod2array) % 2
		temp2  = dot(temp1,A) % 2
		phi    = (C - temp2) % 2
		print 'phi:\n', phi
		# if phi is not an empty matrix or a matrix of 0s, press on
		if phi.any():
			try:
				# try to take the inverse of phi
				invPhi = invMod2(phi)
			except linalg.linalg.LinAlgError:
				# phi is singular
				print 'phi is singular'
			else:
				# phi is nonsingular, so this is our new candidate
				print 'Found betterH woohoo!'
				flagFoundBetterH = 1

				finalH = betterH
				finalGap = gap
				g = gap
				final_t = t
 	else:
 		print 'Gap:', gap, 'is not larger than curent g:', g

if flagFoundBetterH: 
	print 'New H matrix has gap g:', finalGap, 'and t =', final_t
else: 
	print '\n\nDid not find betterH! We have to quit here...\n\n'

############ this is all real-time encoding #########################
if flagFoundBetterH:
	n = finalH.shape[1]
	k = n - finalH.shape[0]  
	s = random_integers(0,1,k).reshape(k,1)
	print 's, filled with k=', k, 'random information bits:\n',s
	print 's.shape:', s.shape

	T = finalH[0:final_t, 0:final_t]
	E = finalH[final_t:final_t+finalGap,0:final_t] #
	A = finalH[0:final_t,final_t:final_t+finalGap] #
	B = finalH[0:final_t,final_t+finalGap:n] #
	D = finalH[final_t:final_t+finalGap,final_t+finalGap:n] #
	invTmod2array = invMod2(T) #

	# compute p1 (this method has lowest complexity)
	a = dot(B,s) % 2
	b = dot(invTmod2array,a) % 2
	c = dot(E,b) % 2
	d = dot(D,s) % 2
	e = d + c % 2
	p2 = dot(invPhi,e) % 2

	# compute p2 (this method has lowest complexity)
	a = dot(A, p2) % 2
	b = dot(B, s) % 2
	c = a + b % 2
	p1 = dot(invTmod2array,c) % 2

	# concatenate to get codeword x
	x = vstack((p1, p2, s))

	# verify:
	testArray = dot(finalH,x) % 2
	if testArray.any():
		print '\n\nCodeword did not pass parity check!!!\n'
	else: 
		print '\n\nCodeword passed partiy check.\n' 