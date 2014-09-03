wignerVille.jl
==============

An implementation of the Wigner-Ville transform in julia
ported from c-code. 

This module is a straight-forward implementation of:

1. The basic Wigner-Ville (WV) transform <br>
2. The pseudo-WV (PWV) <br>
3. And the smoothed-PWD (SPWD)

in the programming language julia.

Initial stage...

The final version calculates the SPWV in parallel using 
open.blas and fftw bindings.


