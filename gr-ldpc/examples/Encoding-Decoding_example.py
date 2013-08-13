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

from LDPC_H_matrix import LDPC_parity_check_matrix
from LDPCpreprocessing import *
import numpy as np
from numpy.random import random_integers
from gnuradio import gr, blocks
from bit_flip_decoder_ss import *
from Richardson_Urbanke_encoder_ss import *

class my_top_block (gr.top_block):
	def __init__(self):
		gr.top_block.__init__(self)

		# FIXME When I try larger matrices, I get allocation errors
		# from the console.

		# Load the parity check matrix.
		# FIXME create an example script to show how to generate
		# the parity check matrix.
		parity_check_matrix=LDPC_parity_check_matrix(alist_filename=\
		          "../../alist_files/H_100_3_5_encoding-ready.alist")

		# The gap g for each matrix is listed in the README file in
		# the alist_files directory.
		g = 2
		t = parity_check_matrix.numRows - g
		H = parity_check_matrix.H
		[invTmod2array,invPhi,E,A,B,D,n,k] = \
		             extractUpperTriangulationMatrixParameters(H,t,g)

		# Simulate actual data by generating an array of random 
		# numbers of length k.
		k = parity_check_matrix.k
		dataword = random_integers(0,1,(1,k))[0].reshape(1,k)
		# Need to make it a tuple
		dataword_tuple = ()
		for index in range(k):
			dataword_tuple=dataword_tuple+(int(dataword[0][index]),)
		
		# Create signal source.
		src = blocks.vector_source_i(dataword_tuple)
		
		# This is a preprocessing step required for the encoder.
		# You don't need to do this repeatedly during real-time
		# encoding. 
		[invT,invPhi,E,A,B,D,n,k] = \
					 extractUpperTriangulationMatrixParameters(H,t,g)

		str2vec = blocks.stream_to_vector(4, k)
		vec2str = blocks.vector_to_stream(4, n)

		# Setup encoder and decoder.
		encoder = Richardson_Urbanke_encoder_ss(invT,invPhi,E,A,B,D,\
		                                        n,k,g)
		decoder = bit_flip_decoder_ss(parity_check_matrix)
		
		# FIXME need to add a channel block to introduce some errors
		# and test the decoder

		# Connect the blocks and run.
		dst = blocks.vector_sink_i()
		self.connect(src, str2vec, encoder, decoder, vec2str, dst)
		self.run ()
		result_data = dst.data()
		
		# FIXME I need to turn this into a block: pull the dataword 
		# from the codeword.
		dataword_result = np.zeros((1,k), int)
		for index in range(k):
			dataword_result[0,index] = result_data[n-k+index]
		
		# Check resulting dataword against original dataword
		test = dataword - dataword_result
		if (test.any()):
			print 'Test failed.'
		else:
			print 'Test passed.'
		
if __name__ == '__main__':
	try:
		my_top_block().run()
	except [[KeyboardInterrupt]]:
		pass
