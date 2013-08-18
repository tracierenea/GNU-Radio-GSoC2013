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
from LDPC_H_matrix import LDPC_parity_check_matrix


class bit_flip_decoder_ss(gr.sync_block):
    """
    Bit flip decoding block
    """
    def __init__(self,LDPC_parity_check_matrix, max_iterations=100):
        self.parity_check_matrix = LDPC_parity_check_matrix
        self.max_iterations = max_iterations
        gr.sync_block.__init__(self,
            name="bit_flip_decoder_ss",
            in_sig=[(np.int32,self.parity_check_matrix.n)],
            out_sig=[(np.int32, self.parity_check_matrix.n)])


    def work(self, input_items, output_items):
        n = self.parity_check_matrix.n
        H = self.parity_check_matrix.H
        codeword = np.zeros((n,1), dtype=int)
        found_valid_codeword = 0
        for i in np.arange(len(input_items[0][0])):
            codeword[i,0] = input_items[0][0][i]
        syndrome =self.calculate_syndrome(H, codeword)
        iteration = 1

        # If the received codeword is valid, no further processing
        # is required. 
        if self.have_match(syndrome):
            iteration = self.max_iterations
            found_valid_codeword = 1

        while iteration < self.max_iterations:
            # For each of the n bits in the codeword, determine how 
            # many of the unsatisfied parity checks involve that bit.
            # To do this: First find the nonzero entries in the 
            # syndrome. The entry numbers correspond to rows of 
            # interest in H.
            rowsToLookAtInH = np.array(np.nonzero(syndrome))[0] 

            # Second, for each bit, determine how many of the 
            # unsatisfied parity checks involve this bit and store 
            # this count in an array. 
            counts = np.zeros_like(codeword)
            for row in rowsToLookAtInH.transpose():
                for bitNumber in np.arange(n):
                    if H[row,bitNumber] > 0:
                        counts[bitNumber] = counts[bitNumber] + 1
            
            # Next, determine which bit(s) is associated with the 
            # most unsatisfied parity checks, and flip it/them
            bitsToFlip = np.where(counts==counts.max())[0]

            for bitNumber in bitsToFlip:
                codeword[bitNumber] = \
                        np.bitwise_xor(int(codeword[bitNumber]),1)

            # note which bits were flipped for check on next iteratio
            lastBitsFlipped = bitsToFlip

            syndrome = self.calculate_syndrome(H,codeword)
            if self.have_match(syndrome):
                # If we've found a valid codeword, stop the loop
                iteration = self.max_iterations
                found_valid_codeword = 1
            else:
                iteration = iteration + 1

        if not found_valid_codeword:
            print 'Error in bit_flip_decoder_ss: valid codeword not'
            print 'found.'

        for i in np.arange(len(output_items[0][0])):
            output_items[0][0][i] = codeword[i,0]
 
        # FIXME use the slice operator when assigning to output items
        # out[:] = in0  
        return len(output_items[0])

    def calculate_syndrome(self,H, codeword):
        """
        Given a codeword and parity check matrix, calculates the
        syndrome.
        """
        syndrome = np.dot(H, codeword) % 2  # Use modulo 2 operations
        return syndrome

    def have_match(self,syndrome):
        """
        Given a syndrome, determines if the associated codeword is
        valid.
        """
        check = 1
        if np.any(syndrome):
            # The matrix is not all zeros, so it's not a codeword
            check = 0
        return check
