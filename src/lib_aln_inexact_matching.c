/* The MIT License

 Copyright (c) 2019 Mattia Marcolin.

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "lib_aln_inexact_matching.h"
#include "bounded_backtracking_seach.h"
#include "fmdindex_load.h"

//Necessary only for GPLv3 version
int fmd_idx_build(const char *fa, const char *prefix, int algo_type, int block_size);

//Derived from bwa_idx_load_from_disk in bwa.c
bwaidx_t* lib_aln_idx_load(const char *path_genome)
{
	bwaidx_t *idx;
	char *prefix;
	prefix = bwa_idx_infer_prefix(path_genome);
	if (prefix == 0)
		return 0;

	idx = calloc(1, sizeof(bwaidx_t));
	idx->bwt = bwa_idx_load_bwt(path_genome);
	if (idx->bwt == 0)
		return 0;
	int i, c;

	idx->bns = bns_restore(prefix);
	if (idx->bns == 0)
		return 0;
	for (i = c = 0; i < idx->bns->n_seqs; ++i)
		if (idx->bns->anns[i].is_alt)
			++c;

	idx->pac = calloc(idx->bns->l_pac / 4 + 1, 1);
	err_fread_noeof(idx->pac, 1, idx->bns->l_pac / 4 + 1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
	err_fclose(idx->bns->fp_pac);
	idx->bns->fp_pac = 0;

	free(prefix);
	return idx;
}

//Derived from bwa_idx_destroy in bwa.c
void lib_aln_idx_destroy(bwaidx_t *idx)
{
	if (idx == 0)
		return;

	free(idx->bwt);
	free(idx->bns->anns);
	free(idx->bns);
	free(idx);
}

search_result** lib_aln_bound_backtracking(const bwaidx_t *idx, const char* pattern_input, uint8_t max_mismatches, uint32_t* numHit,
											const uint8_t type_search, const uint8_t type_output)
{
	//Check that the input parameters are valid
	controlParam(idx, pattern_input, type_search, type_output);

	//Size of pattern input
	size_t pattern_len = strlen(pattern_input);

	if (pattern_len < max_mismatches)
		max_mismatches = pattern_len;

	/*
	 * The pattern to be searched is the r.c. of the input string and each base is converted to integer.
	 * A=0, C=1, G=2, T=3
	 */
	uint8_t* pattern_to_search = create_pattern_to_search(pattern_input, pattern_len);

	//Store all information regard the search in seq
	input_query* seq = init_bwt_seq(pattern_to_search, pattern_len, max_mismatches, type_search, type_output);

	//Core method that solves the approximate pattern matching problem
	search_result** returnM = get_approximate_match(idx, seq, numHit);

	free(pattern_to_search);
	free(seq);

	return returnM;
}

void lib_aln_sr_destroy(search_result** result, uint32_t numHit)
{
	for (int i = 0; i < numHit; i++)
		free(result[i]);

	free(result);
}

//To improve and only for  GPLv3 version..
void lib_aln_index(const char* path_genome, const char* prefix, int algo_type)
{
	if (path_genome == 0)
	{
		fprintf(stderr, "Miss pattern to search.\n");
		exit(EXIT_FAILURE);
	}
	else if (path_genome == 0)
	{
		fprintf(stderr, "Miss prefix.\n");
		exit(EXIT_FAILURE);
	}
	else if (algo_type < 0 || algo_type > 4)
	{
		fprintf(stderr, "Miss type of algorithm to apply.\n");
		exit(EXIT_FAILURE);
	}

	fmd_idx_build(path_genome, prefix, algo_type, 10000000); //10000000 is block_size
}

