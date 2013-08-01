#!/usr/bin/env python
# 
# Copyright 2013 Tracie Perez.
# 
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published 
# by the Free Software Foundation; either version 3, or (at your 
# option) any later version.
# 
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

class LDPC_parity_check_matrix:
	""" Class for a LDPC parity check matrix """
	def __init__(self, alist_filename = None, n_p_q = None ):
		if (alist_filename != None):
			self.H = read_alist_file(alist_filename)
		elif (n_p_q != None):
			self.H = regular_LDPC_code_contructor(n_p_q)
		else:
			print 'Error: provide either an alist filename or', 
			print 'parameters for constructing regular LDPC parity',
			print 'check matrix.'

		self.rank = linalg.matrix_rank(H)
		self.numRows = H.shape[0]
		self.n = H.shape[1]
		self.k = self.n -self.numRows


	def read_alist_file(filename):
		# This function reads in an alist file and returns the
		# corresponding parity check matrix H. The format of alist
		# files is desribed at:
		# http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html

		myfile = open(filename,'r')
		data = myfile.readlines()

		size      = string.split(data[0])
		numRows   = int(size[0])
		numCols   = int(size[1])
		weights   = string.split(data[1])
		rowWeight = weights[0]
		colWeight = weights[1]

		H = zeros((numRows,numCols))
		for lineNumber in arange(4,4+numRows):
			indices = string.split(data[lineNumber])
			for index in indices:
				H[lineNumber-4,int(index)-1] = 1

		# the subsequent lines in the file list the indices for where
		# the 1s are in the columns, but this is redundant 
		# information

		return H

	def regular_LDPC_code_contructor(n_p_q):

		# Following Gallager's approach where we create p submatrices. 
		# Reference: Turbo Coding for Satellite and Wireless 
		# Communications, sec 9.3

		# Note, the matrices computed from this algorithm will never
		# have full rank. Reference Gallager's Dissertation. They 
		# will have rank = (number of rows - p + 1). To convert it
		# to full rank, use the function getFullRankHmatrix

		n = n_p_q[0] 	# codeword length
		p = n_p_q[1] 	# column weight
		q = n_p_q[2]	# row weight
		# FIX: There should probably be other guidelines for n/p/q,
		# but I have not found any specifics in the literature....

		# For this algorithm, n/p must be an integer, because the
		# number of rows in eacn submatrix must be a whole number.
		ratioTest = (n*1.0)/q
		if ratioTest%1 != 0:
			print 'Error in regular_LDPC_code_contructor: The'
			print 'ratio of inputs n/q must be a whole number.'
			return

		# First submatrix first: 
		m = (n*p)/q  # number of rows in H matrix
		submatrix1 = zeros((m/p,n))  
		for row in arange(m/p):
			range1 = row*q
			range2 = (row+1)*q 
			submatrix1[row,range1:range2] = 1
		H = submatrix1

		# Create the other submatrices and vertically stack them on.
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

		# Double check the row weight and column weights.
		size = H.shape
		rows = size[0]
		cols = size[1]

		# Check the row weights.
		for rowNum in arange(rows):
			nonzeros = array(H[rowNum,:].nonzero())
			if nonzeros.shape[1] != q:
				print 'Row', rowNum, 'has incorrect weight!'
				return

		# Check the column weights
		for columnNum in arange(cols):
			nonzeros = array(H[:,columnNum].nonzero())
			if nonzeros.shape[1] != p:
				print 'Row', columnNum, 'has incorrect weight!'
				return

		return H
		