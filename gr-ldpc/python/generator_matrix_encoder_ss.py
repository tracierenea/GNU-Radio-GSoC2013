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

import numpy as np
from gnuradio import gr

class generator_matrix_encoder_ss(gr.sync_block):
    """
    This is a standard encoding method that is discussed in many
    error-correcting textbooks. One is: Turbo Coding for Satellite
    and Wireless Communications by Soleymani, Gao, and 
    Vilaipornsawai. 

    The generator matrix G must be in systematic
    form. Use function getSystematicGmatrix if necessary.
    """
    def __init__(self,G):
        self.G = G
        self.k = G.shape[0]
        self.n = G.shape[1]
        if np.any(G[0:self.k,0:self.k]-np.identity(self.k)):
            print 'Error: G not in systematic form.'
        gr.sync_block.__init__(self,
            name="generator_matrix_encoder_ss",
            in_sig=[(np.int32,self.k)],
            out_sig=[(np.int32,self.n)])

    def work(self, input_items, output_items):
        # s is the word of k information bits to be encoded
        s = np.zeros((self.k,1), dtype=int)
        for i in np.arange(len(input_items[0][0])):
            s[i,0] = input_items[0][0][i]

        if self.k != s.shape[0]:
            print 'Error: provided message is not the right size.'
            return
        
        codeword = np.dot(self.G.transpose(),s) % 2

        # FIXME use the slice operator when assigning to output items
        # out[:] = in0

        for i in np.arange(len(output_items[0][0])):
            output_items[0][0][i] = codeword[i,0]        

        return len(output_items[0])
