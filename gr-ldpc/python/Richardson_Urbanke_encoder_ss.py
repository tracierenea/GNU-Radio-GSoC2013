#!/usr/bin/env python
# 
# Copyright 2013 Tracie Perez
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

class Richardson_Urbanke_encoder_ss(gr.sync_block):
    """
    This encoding method is from Richardson and Urbanke's Modern 
    Coding Theory, Appendix A. Arguments to the constructor 
    should be taken care of during a preprocessing step using the 
    function getParametersForEncoding so that the speed of real-time 
    encoding is minimized. 
    """
    def __init__(self,invT,invPhi,E,A,B,D,n,k,g):
        self.inverseT   = invT
        self.inversePhi = invPhi
        self.E = E
        self.A = A
        self.B = B
        self.D = D
        self.n = n
        self.k = k
        self.gap = g
        gr.sync_block.__init__(self,
            name="Richardson_Urbanke_encoder_ss",
            in_sig=[(np.int32,self.k)],
            out_sig=[(np.int32,self.n)])

    def work(self, input_items, output_items):
        # s is the word of k information bits to be encoded
        s = np.zeros((self.k,1), dtype=int)
        for i in np.arange(len(input_items[0][0])):
            s[i,0] = input_items[0][0][i]

        # compute p1 and p2 (this method has lowest complexity)
        a = np.dot(self.B, s) % 2
        b = np.dot(self.inverseT, a) % 2
        c = np.dot(self.E, b) % 2
        d = np.dot(self.D, s) % 2
        e = d + c % 2
        p2 = np.dot(self.inversePhi, e) % 2
        a = np.dot(self.A, p2) % 2
        a.resize((self.n-self.k-self.gap), 1)
        b = np.dot(self.B, s) % 2
        c = a + b % 2
        p1 = np.dot(self.inverseT, c) % 2

        # concatenate to get codeword
        codeword = np.vstack((p1, p2, s))

        # need to add this to account for rounding errors seeing in
        # mod 2 operations. These margins are pretty generous.
        tempTest = zeros_like(codeword)
        for colNum in arange(codeword.shape[1]):
            for rowNum in arange(codeword.shape[0]):
                value = codeword[rowNum,colNum]
                if (abs(1-value)) < 0.2:
                    tempTest[rowNum,colNum] = 1
                elif (abs(2-value)) < 0.2:
                    tempTest[rowNum,colNum] = 0
                elif (abs(0-value)) < 0.2:
                    tempTest[rowNum,colNum] = 0
        codeword = tempTest.copy()

        # FIXME use the slice operator when assigning to output items
        # out[:] = in0

        for i in np.arange(len(output_items[0][0])):
            output_items[0][0][i] = codeword[i,0]
        return len(output_items[0])
