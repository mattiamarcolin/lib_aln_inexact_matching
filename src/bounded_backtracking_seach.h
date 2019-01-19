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

#ifndef _BACKTRACKING_SEARCH_H
#define _BACKTRACKING_SEARCH_H

#include "occ.h"
#include "sa.h"

#ifndef BWAIDX_T
#define BWAIDX_T

typedef struct
{
	bwt_t *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

#endif

//TYPE_SEARCH
#ifndef TYPE_SEARCH_
#define TYPE_SEARCH_

#define ARBITRARY_HIT 0x0
#define ALL_HITS 0x1
//#define UNIQUE_BEST_HIT[TO DO]
//#define ALL_BEST_HIT[TO DO]

#endif

//TYPE_OUTPUT
#ifndef TYPE_OUTPUT_
#define TYPE_OUTPUT_

#define ALLOW_REV_COMP 0x0
#define NO_REV_COMP 0x1

#endif

typedef struct
{
	ubyte_t *seq;
	uint32_t len :20, type_search :1, type_output :1;
	uint8_t max_diff;
} input_query;

typedef struct
{
	uint64_t w;
	int bid;
} bwt_width_t;

typedef struct
{
	uint32_t info; //i
	uint8_t n_mm :8;
	int last_diff_pos;
	uint64_t k, l; // (k,l) is the SA region of [i,n-1]
} entry_t;

typedef struct
{
	uint32_t n_entries_substack, m_entries;
	entry_t *start_substack;
} substack_t;

typedef struct
{
	int n_sub_stacks, best, n_entries;
	substack_t *stacks;
} stack_t;

#ifndef SEARCH_RESULT
#define SEARCH_RESULT

typedef struct
{
	char* hit;
	uint64_t* positions_to_ref;
	uint32_t num_occur;
	uint32_t* different_positions;
	uint8_t n_mismatches;
	bool is_rev_comp;
} search_result;

#endif

#ifdef __cplusplus
extern "C"
{
#endif

	void controlParam(const bwaidx_t *idx, const char* pattern_input, const uint8_t type_search, const uint8_t type_output);

	input_query* init_bwt_seq(uint8_t* pattern_to_search, const size_t pattern_len, const uint8_t nmismatch, const uint8_t type_search,
								const uint8_t type_output);

	uint8_t* create_pattern_to_search(const char* pattern_input, const size_t pattern_len);

	search_result** get_approximate_match(const bwaidx_t* idx, input_query* info_seq, uint32_t* numHit);
#ifdef __cplusplus
}
#endif

#endif
