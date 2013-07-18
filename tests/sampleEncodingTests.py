#!/usr/bin/python
from LDPCfunctions import *
from numpy import *
from numpy.random import random_integers

# Example of how to do encoding per Appendix A in Modern Coding
# Theory (Richardson/Urbanke). 

# The H matrices created by regularLDPCcodeconstructor are
# not full rank. Using getSystematicGmatrix to get a version that is
# full rank, then sending it on to the encoding algorithm 

H = regularLDPCcodeConstructor(1920,4,6)
# H = regularLDPCcodeConstructor(20,3,4)
print 'Origiinal H.shape:', H.shape
printArrayToFile(H,'H_1920_4_6_original.txt')
[G, newH] = getSystematicGmatrix(H)
printArrayToFile(G,'G_from_1920_4_6_H_matrix.txt')
print 'linalg.matrix_rank(newH):', linalg.matrix_rank(newH)

# Hp = newH[0:newH.shape[0],0:newH.shape[0]]
# print 'Hp.shape:', Hp.shape
# print 'Inverse of Hp?\n', invMod2(Hp)

############ this is all preprocessing #############################

numIterations = 20
print 'Running permutation algorithm', numIterations,
print 'times to ensure the lowest \ngap is found and the',
print 'resulting phi is nonsingular.'

# set this arbitrarily high to force whole loop on first run
g = 10**10
flagFoundBetterH = 0

for index in arange(numIterations):
 	print '============== Index:', index
	[betterH, gap, t]  = greedyUpperTriangulation(newH)

	if gap < g:
		n = betterH.shape[1]
		T = betterH[0:t, 0:t]
		E = betterH[t:t+gap,0:t]
		A = betterH[0:t,t:t+gap]
		C = betterH[t:t+gap,t:t+gap]
		invTmod2array = invMod2(T)
		temp1  = dot(E,invTmod2array) % 2
		temp2  = dot(temp1,A) % 2
		phi    = (C - temp2) % 2
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
	print 's, filled with k=', k, 'random information bits:\n'
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