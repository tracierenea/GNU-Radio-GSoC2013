#!/usr/bin/python

from numpy import *

verbose = 0

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

	print G
	P = G[:,k:n]
	print 'P:\n',P
	print G.transpose().shape
	print s.shape
	t = dot(G.transpose(),s) % 2
	print 't:\n',t

	return t