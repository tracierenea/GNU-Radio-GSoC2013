/* -*- c++ -*- */

#define LDPC_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "ldpc_swig_doc.i"

%{
#include "ldpc/decoder_bit_flip_ss.h"
%}


%include "ldpc/decoder_bit_flip_ss.h"
GR_SWIG_BLOCK_MAGIC2(ldpc, decoder_bit_flip_ss);
