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
from bit_flip_decoder_ss import bit_flip_decoder_ss
from LDPC_H_matrix import LDPC_parity_check_matrix

class qa_bit_flip_decoder_ss (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
    	# This test case is from Fundamentals of Error-Correcting
    	# Codes by Huffman and Pless. Example 15.6.1.

        parity_check_matrix = LDPC_parity_check_matrix(
            alist_filename="qa_bit_flip_decoder_ss_test_001_t.alist")
    	# transmitted codeword is the truth data
        transmitted_codeword = (1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0)
        # received codeword contains two bits flipped (2 errors)
        received_codeword    = (1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0)
        src = blocks.vector_source_i(received_codeword)
        str2vec = blocks.stream_to_vector(4, parity_check_matrix.n)
        vec2str = blocks.vector_to_stream(4, parity_check_matrix.n)
        decoder = bit_flip_decoder_ss(parity_check_matrix)
        dst = blocks.vector_sink_i()
        self.tb.connect(src, str2vec, decoder, vec2str, dst)
        self.tb.run ()
        result_data = dst.data()
        # check data
        self.assertTupleEqual(transmitted_codeword,result_data)

    def test_002_t (self):
        # This test case is from Fundamentals of Error-Correcting
        # Codes by Huffman and Pless. Example 15.6.1.

        parity_check_matrix = LDPC_parity_check_matrix(
            alist_filename="qa_bit_flip_decoder_ss_test_001_t.alist")
        # transmitted codeword is the truth data
        transmitted_codeword = (1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0)
        # received codeword contains two bits flipped (2 errors)
        received_codeword    = (0,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0)
        src = blocks.vector_source_i(received_codeword)
        str2vec = blocks.stream_to_vector(4, parity_check_matrix.n)
        vec2str = blocks.vector_to_stream(4, parity_check_matrix.n)
        decoder = bit_flip_decoder_ss(parity_check_matrix)
        dst = blocks.vector_sink_i()
        self.tb.connect(src, str2vec, decoder, vec2str, dst)
        self.tb.run ()
        result_data = dst.data()
        # check data
        self.assertTupleEqual(transmitted_codeword,result_data)

    def test_003_t (self):
        # This test case is from Fundamentals of Error-Correcting
        # Codes by Huffman and Pless. Example 15.6.1.

        parity_check_matrix = LDPC_parity_check_matrix(
            alist_filename="qa_bit_flip_decoder_ss_test_001_t.alist")
        # transmitted codeword is the truth data
        transmitted_codeword = (1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0)
        # received codeword contains one bit flipped (1 error)
        received_codeword    = (1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,1)
        src = blocks.vector_source_i(received_codeword)
        str2vec = blocks.stream_to_vector(4, parity_check_matrix.n)
        vec2str = blocks.vector_to_stream(4, parity_check_matrix.n)
        decoder = bit_flip_decoder_ss(parity_check_matrix)
        dst = blocks.vector_sink_i()
        self.tb.connect(src, str2vec, decoder, vec2str, dst)
        self.tb.run ()
        result_data = dst.data()
        # check data
        self.assertTupleEqual(transmitted_codeword,result_data)    

if __name__ == '__main__':
    gr_unittest.run(qa_bit_flip_decoder_ss, "qa_bit_flip_decoder_ss.xml")
