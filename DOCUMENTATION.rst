======================================= 
lib_aln_inexact_matching: documentation
=======================================

.. contents ::

Definitions
===========

type search
-----------

The library allows you to operate the following two search modes.

- **ARBITRARY_HIT**: The search ends when a legal hit is found.     
- **ALL_HITS**: The research space is exhaustively analyzed.

type output
-----------
Defines whether among the hits returned by the algorithm there may be some that must be considered the reverse complement.

- **ALLOW_REV_COMP**: In the output returned by the algorithm, there may be hits that must be considered reverse complement.
- **NO_REV_COMP**: In the output returned by the algorithm, there aren't hits that must be considered the reverse complement.

**NOTE**: To be careful, because of the use of the FMD-index, the algorithm is not faster if the reverse complement is not allowed.

type bwt index algorithm
------------------------

Define the algorithm for constructing BWT index:

- **BWTALGO_AUTO**: The best algorithm is performed between *BWTALGO_BWTSW* and *BWTALGO_IS* based on the reference genome.
- **BWTALGO_BWTSW**: Algorithm implemented in BWT-SW. This method works with the whole human genome, but it does not work with database smaller than 10MB and it is usually slower than IS.
- **BWTALGO_IS**: IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast, but does not work with database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for IS algorithm are reimplemented by Yuta Mori.

**NOTE**: This information are present inside `BWA documentation`_.

search_result
-------------
Data structure containing all the information regarding each hit found::

    typedef struct 
    {
        char* hit;
	uint64_t* positions_to_ref;
	uint32_t num_occur;
	uint32_t* different_positions;
	uint8_t n_mismatches;
	bool is_rev_comp;
     } search_result;

:hit: Pointer to the first character of the hit.     
:positions_to_ref: Pointer to the first initial position of hit inside reference genome.
:num_occur: Number of hit's occurrences.
:different_positions: Pointer to the first position where the hit and the input patter differ.
:n_mismatches: Number of mismatches between the hit and the input pattern.
:is_rev_comp: This field assume value false if hit must be considered forward, true if hit must be considered reverse complement. 
 
Functions
=========

lib_aln_index
----------------

Index database sequences in the FASTA format::

    void lib_aln_index(char* path_genome, char* prefix, int algo_type);

:path_genome: Path where is locate database sequences in the FASTA format.
:prefix: Prefix of the output database.
:algo_type: Index construction algorithm. Permissible value are *BWTALGO_AUTO*, *BWTALGO_BWTSW* and *BWTALGO_IS*. 

lib_aln_idx_load
----------------

Loads in RAM the FMD-index of the reference genome::

    bwaidx_t* lib_aln_idx_load(const char *path_genome)

:path_genome: Path where is locate database sequences in the FASTA format.  

lib_aln_idx_destroy
-------------------

Deallocates memory uses by the index::

 void lib_aln_idx_destroy(bwaidx_t *idx);

:idx: FMD-index related to the reference genome.

lib_aln_bound_backtracking
--------------------------

Method that solve approximate pattern matching problem::

 search_result** lib_aln_bound_backtracking(const bwaidx_t *idx, const char* pattern_input, 
				            const uint8_t max_mismatches, uint32_t* num_hits,
					    const uint8_t type_search, const uint8_t type_output);

:idx: FMD-index related to the reference genome.
:pattern_input: Patter to be searched in the index.
:max_mismatches: Max number of mismatches between hit found by algorithm and pattern_input.
:num_hits: Number of admissible hits found.
:type_search: Defines the type of search. Permissible value are *ARBITRARY_HIT* and *ALL_HITS*.
:type_output: Defines the type of output. Permissible value are *ALLOW_REV_COMP* and *NO_REV_COMP*.

lib_aln_sr_destroy
------------------

Deallocates memory use to store the results of the search::

   void lib_aln_sr_destroy(search_result** result, uint32_t num_hits);

:result: Output of method lib_aln_bound_backtracking.
:num_hits: Number of hits inside result.

.. _BWA documentation: http://bio-bwa.sourceforge.net/bwa.shtml


