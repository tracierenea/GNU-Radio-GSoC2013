#!/usr/bin/python
from numpy import *

def matrixMultiplierEncoder(G,s):
	# This function accepts a generator matrix in systematic form
	# [I P] and the message s, and outputs the codeword. Using 
	# notation from Turbo Coding for Satellite and Wireless
	# Communication by Soleymani, Gao, and Vilaipornsawai
	shapeG = G.shape
	k = shapeG[0]
	n = shapeG[1]
	print G[0:k,0:k]
	if any(G[0:k,0:k]-identity(k)):
		if verbose: 'Error: G not in systematic form.'
		return
	if shapeG[0] != s.shape[0]:
		print 'Error: provided message is not the right size.'
		return
		
	t = dot(G.transpose(),s) % 2

	return t

def RichardsonUrbankeEncoder(invT,invPhi,E,A,B,D,bestH,n,k,g,s):
	# This method is from Richardson and Urbanke's Modern Coding
	# Theory, Appendix A

	# compute p1 and p2 (this method has lowest complexity)
	a = dot(B,s) % 2
	b = dot(invT,a) % 2
	c = dot(E,b) % 2
	d = dot(D,s) % 2
	e = d + c % 2
	p2 = dot(invPhi,e) % 2
	a = dot(A, p2) % 2
	a.resize((n-k-g),1)
	b = dot(B, s) % 2
	c = a + b % 2
	p1 = dot(invT,c) % 2

	# concatenate to get codeword x
	x = vstack((p1, p2, s))

	return x
