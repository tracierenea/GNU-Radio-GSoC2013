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

# This is an example of how to generate a parity check matrix prior 
# to real-time encoding using the Richardson_Urbanke_encoder_ss block

# First, generate a regular LDPC parity check matrix. Specify
# the properties desired. For example:

n = 1500	# number of columns, corresponds to codeword length
p = 3		# column weight
q = 5 		# row weight

parity_check_matrix = LDPC_parity_check_matrix(n_p_q = [n,p,q])

# In order to use the Richardson_Urbanke_encoder_ss, the matrix must
# be full rank. The matrices generated by the 
# regular_LDPC_code_contructor function will never be full rank. So,
# use the getFullRankHmatrix function. 

newH = getFullRankHmatrix(parity_check_matrix.H)

# At this point, the matrix is no longer regular.

# Next, some preprocessing steps need to be performed as described
# Richardson and Urbanke in Modern Coding Theory, Appendix A. This
# can take a while...

[invT,invPhi,E,A,B,D,bestH,n,k,g] =getParametersForEncoding(newH,100)

# Print out some of the resulting properties. You will need the gap 
# value later when setting up the flow graph.
print 'Completed preprocessing steps. Parity check matrix',
print 'properties:'
print 'Size of H:', bestH.shape
print 'Rank of H:', linalg.matrix_rank(bestH)
print 'n:', n, '\tk:', k, '\tgap:', g
print 'Rate (k/n):', (k*1.0)/n

# Finally, save the matrix to an alist file for future use: 
alist_filename = "H_"+`n`+"_"+`p`+"_"+`q`+"_encoding-ready.alist"
parity_check_matrix.writeAlistFile(alist_filename,bestH)
print 'Matrix saved to alist file:', alist_filename

# Now refer to Encoding-Decoding_example.py for how to put together
# the flow graph.