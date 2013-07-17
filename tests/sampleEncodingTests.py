#!/usr/bin/python
from LDPCfunctions import *
from numpy import *
from numpy.random import random_integers

# Example of how to do encoding per Appendix A in Modern Coding
# Theory (Richardson/Urbanke). 

# The H matrices created by regularLDPCcodeconstructor are
# not full rank. Using getSystematicGmatrix to get a version that is
# full rank, then sending it on to the encoding algorithm 

# H = regularLDPCcodeConstructor(1920,4,6)
H = regularLDPCcodeConstructor(20,3,4)
print 'Origiinal H.shape:', H.shape
printArrayToFile(H,'H_20_3_4_original.txt')
[G, newH] = getSystematicGmatrix(H)
printArrayToFile(G,'G_from_20_3_4_H_matrix.txt')
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
		B = betterH[0:t,t+gap:n]
		D = betterH[t:t+gap,t+gap:n]
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
 	else:
 		print 'Gap:', gap, 'is not larger than curent g:', g

if flagFoundBetterH: 
	print 'New H matrix has gap g:', gap, 'and t =', t
	print 'T:\n', T, '\n', T.shape
	print 'E:\n', E, '\n', E.shape
	print 'A:\n', A, '\n', A.shape
	print 'C:\n', C, '\n', C.shape
	print 'B:\n', B, '\n', B.shape
	print 'D:\n', D, '\n', D.shape
else: 
	print '\n\nDid not find betterH! We have to quit here...\n\n'

############ this is all real-time encoding #########################
if flagFoundBetterH:
	n = newH.shape[1]
	k = n - newH.shape[0] 
	s = random_integers(0,1,k).reshape(k,1)
	print 's, filled with k=', k, 'information bits:\n', s
	print 's.shape:', s.shape

	# compute p1 (this method has lowest complexity)
	print 'invPhi.shape:', invPhi.shape
	a = dot(B,s) % 2
	print 'a.shape:', a.shape
	b = dot(invTmod2array,a) % 2
	print 'b.shape:', b.shape
	c = dot(E,b) % 2
	print 'c.shape:', c.shape
	d = dot(D,s) % 2
	print 'd.shape:', d.shape
	e = d + c % 2
	print 'e.shape:', e.shape
	p2 = dot(invPhi,e) % 2
	print 'p2.shape:', p2.shape
	print 'p2:\n', p2

	# compute p2 (this method has lowest complexity)
	a = dot(A, p2) % 2
	b = dot(B, s) % 2
	c = a + b % 2
	p1 = dot(invTmod2array,c) % 2
	print 'p1:\n', p1

	# concatenate to get codeword x
	x = vstack((p1, p2, s))
	print 'Codeword x:\n', x

	# verify:
	print 'This should be zeros:\n', dot(betterH,x) % 2
	