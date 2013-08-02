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

from gnuradio import gr, gr_unittest, blocks
from generator_matrix_encoder_ss import generator_matrix_encoder_ss
import numpy as np

class qa_generator_matrix_encoder_ss (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
    	# This test case is from Fundamentals of Error-Correcting
    	# Codes by Huffman and Pless. First part of Example 15.7.1.

    	# Generator matrix from textbook example
    	G = np.array([[1,0,0,0,0,1,1,1],
	                  [0,1,0,0,1,0,1,1],
	                  [0,0,1,0,1,1,0,1],
	                  [0,0,0,1,1,1,1,0]])
    	# this is the dataword to be encoded
    	dataword = (1,0,1,1)
    	# this is what the encoder should produce
    	codeword = (1,0,1,1,0,1,0,0)
    	k = len(dataword)
    	n = len(codeword)
    	src = blocks.vector_source_i(dataword)
        str2vec = blocks.stream_to_vector(4, k)
        vec2str = blocks.vector_to_stream(4, n)
        encoder = generator_matrix_encoder_ss(G)
        dst = blocks.vector_sink_i()
        self.tb.connect(src, str2vec, encoder, vec2str, dst)
        self.tb.run ()
        result_data = dst.data()
        # check data
        self.assertTupleEqual(codeword,result_data)


if __name__ == '__main__':
    gr_unittest.run(qa_generator_matrix_encoder_ss, "qa_generator_matrix_encoder_ss.xml")
