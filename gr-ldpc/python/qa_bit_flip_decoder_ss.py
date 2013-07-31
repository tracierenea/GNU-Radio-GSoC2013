#!/usr/bin/env python
# 
# Copyright 2013 <+YOU OR YOUR COMPANY+>.
# 
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
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

from gnuradio import gr, gr_unittest
from bit_flip_decoder_ss import bit_flip_decoder_ss

class qa_bit_flip_decoder_ss (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
    	# transmitted codeword is the truth data
        transmitted_codeword = (1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0)
        # received codeword contains two bits flipped (2 errors)
        received_codeword    = (1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0)
        # source of shorts
        src = gr.vector_source_s(received_codeword)
        # solution is the codeword found by decoder
        solution = qa_bit_flip_decoder_ss()
        # short sink that writes to a vector
        dst = gr.vector_sink_s()

        self.tb.connect(Src,solution)
        self.tb.connect(solution,dst)
        self.tb.run ()
        result_data = dst.data()

        # check data
        self.assertTupleEqual(transmitted_codeword,result_data)



if __name__ == '__main__':
    gr_unittest.run(qa_bit_flip_decoder_ss, "qa_bit_flip_decoder_ss.xml")
