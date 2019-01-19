lib_aln_inexact_matching 
========================

lib_aln_inexact_matching library allows you to use the algorithm implemented in the aligner tool `BWA-backtrack`_ to solve the `approximate string matching problem`_.

Description
-----------

`Burrows-Wheeler Aligner`_ (BWA) is a software created by bioinformatic scientists `Heng Li`_ and `Richard Durbin`_ to performs `local alignment`_ of short reads on a reference genome. 
It's widely used in bioinformatics and consists of the following three independent algorithms: BWA-backtrack, BWA-SW and BWA-MEM. 
I refer to the `BWA documentation`_ for more information on which algorithm to choose based on the reads dataset that you want to align.

As the name of the software suggests, all three algorithms use the `Burrows-Wheeler transform`_ (BWT) as the basic data structure. 
The BWT, in addition to other auxilliary data structures, is a space-efficient way to index a genome(more generally any string). 
This new data structure is called `FM-Index`_ or BWT-index. 

To be precise, BWA uses a variant of the FM-Index called `FMD-Index`_, used almost exclusively in bioinformatics. 
FMD-Index of a DNA sequence is the FM-index of the concatenation of the sequence and its reverse complement.
For practical purposes, this allows forward and backward strand aligned simultaneously, which is faster than aligning the two strands separately. 

BWA-backtrack is the original alignment algorithm from the BWA program’s release in 2009. The core of the aligner is based on a very efficient 
algorithm that solves the problem of inexact matching.  
Considering that very often we need to search a string *S* (for example a k-mer) within a reference whit at most *m* mismatches, I created this library 
to make this algorithm available, eliminating everything concerning the alignment present in the original version.

Moreover, thanks to the FMD-Index, it's possible to efficiently simultaneously not only search *S* with at most *m* mismatches, but also its reverse complement.

Requirements
------------

- CMake
- zlib

Usage
-----

Include `lib_aln_inexact_matching.h` in your source code and compile it with `lib_aln_inexact_matching` library.

Please see `documentation`_ for details concerning the public interface of the library.

Tests
-----

Inside the **src** directory, is present the file example.c which simulates the use of the library.
If you want to run test, please go to directory **lib_aln_inexact_matching** and then use `make` to complile.
Next, run the executable `example`.

**NOTE** 

- The reference genome and the index that will be created are located within the directory **data**.
- For simplicity, the pattern to be searched and the number of mismatches are defined directly inside example.c.
 	 
Citation
--------

* Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*, **25**, 1754-1760. [PMID:[`19451168`_]]. 

License
-------

Copyright (c) 2019 Mattia Marcolin.

'lib_aln_inexact_matching' is released under `GPLv3`_.

.. _Burrows-Wheeler Aligner: https://github.com/lh3/bwa
.. _local alignment: https://en.wikipedia.org/wiki/Sequence_alignment
.. _approximate string matching problem: https://en.wikipedia.org/wiki/Approximate_string_matching
.. _BWA documentation: http://bio-bwa.sourceforge.net/bwa.shtml
.. _Burrows-Wheeler transform: http://www.hpl.hp.com/techreports/Compaq-DEC/SRC-RR-124.pdf
.. _FM-Index: https://ieeexplore.ieee.org/abstract/document/892127
.. _FMD-Index: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3389770/
.. _BWA-backtrack: https://www.ncbi.nlm.nih.gov/pubmed/19451168
.. _GPLv3: https://en.wikipedia.org/wiki/GNU_General_Public_License
.. _Heng Li: http://lh3lh3.users.sourceforge.net/
.. _Richard Durbin: https://www.sanger.ac.uk/people/directory/durbin-richard
.. _documentation: https://github.com/mattiamarcolin/lib_aln_inexact_matching/blob/master/DOCUMENTATION.rst
.. _19451168: https://www.ncbi.nlm.nih.gov/pubmed/19451168


