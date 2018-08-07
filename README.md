# Lumos
Dedicated repository for the Lumos scheduling technique implementation and evaluation scripts used in this paper:

**LUMOS: A Fast and Efficient Optical Circuit Switch Scheduling Technique**

Ariel Livshits and Shay Vargaftik

*IEEE Communication Letters, 2018.*

[Link to IEEE Digital Library Location](https://ieeexplore.ieee.org/document/8423619/)

Please cite the appropriate paper if you use this code.

The repository contains implementations of four OCS scheduling algorithms:

  -- Birkhoff–von Neumann Decomposition (BvN)*
  
  -- Solstice 
  
  -- Eclipse 
  
  -- Lumos
  
  **C Implemntation:**
  
  Files:
   SMG.c/.h - skewed and sparse matrix generator (synthetic input)
    
   Algorithms.c/.h - implementation of BvN, Solstice, Eclipse and Lumos
    
   HopcroftKarp.c/.h - implementation of the Hopcroft-Karp bipartite matching algorithm [Wikipedia Article](https://en.wikipedia.org/wiki/Hopcroft%E2%80%93Karp_algorithm)
    
   hungarian.c/.h - implementation of the linear minimum assignment (Hungarian algorithm) [Wikipedia Article](https://en.wikipedia.org/wiki/Hopcroft%E2%80%93Karp_algorithm) credit to: Cyrill Stachniss, 2004 [GitHub Repository](https://github.com/losvald/libhungarian)
    
   lms_structs.c/.h and Resources.h - implementation of the 'lms_mat_t' sparse matrix representation. Denote the number of non-zero entries in the sparse matrix as 'nz', then every operation of 'lms_mat_t' has a theoretical running time of O(nz).
    
   main.c - the 'main' function of the program 'decomp'. 'main' recieves parameters from 'analysis.py', then calls upon the sparse matrix generator to get the input, then calls upon the scheduling algorithms to decompse the matrix. Finally, 'main' measures CPU run-time (in nano-seconds), as well as decomposition length and demand completion time. The output is a file containing mean and standard deviation of the three measurements, mentioned above, over the number of trails given.
    
**analysis.py (python 2.7):** 

Sets the parameters for 'decomp'. For example:

  	params['trials'] = 100
  	params['delta'] = 25
  	params['matrix_dim'] = [96, 160, 384]
  	params['small_flow_num'] = 48
  	params['large_flow_num'] = 16
  	params['small_coeff'] = [1,16]
  	params['large_coeff'] = [16,100]
    
'trails' - number of matrices to decompose

'delta' - the switch reconfiguration time penalty

'matrix_dim' - list of switch radices (data will be collected for 'trails' per switch radix)

'small_flow_num' - number of small flows per input port of switch

'large_flow_num' - number of large flows per input port of switch

'small_coeff' - range from which the small flow demands are chosen (uniformally at random)

'large_coeff' - range from which the large flow demands are chosen (uniformally at random)
    
The output of the script is a '.db' file generated by the python 'shelve' package which contains all data output of 'decomp'
Required python packages: numpy and shelve

makefile - makefile for the target program: 'decomp'

(*) - BvN was not used in the evaluation of this paper, however it is a known baseline for the demand covering problem.
