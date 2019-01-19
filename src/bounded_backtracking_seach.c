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
#include <inttypes.h>
#include <stdbool.h>
#include "bounded_backtracking_seach.h"
#include "sa.h"
#include "occ.h"
#include "bntseq.h"

//Only for debug: 0 no output,3 output
int verbose_bound_backtracking_search = 0;

static int cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);

static stack_t * init_stack(int nmismatch);
static void reset_stack(stack_t *stack);
static void destroy_stack(stack_t *stack);

static inline void pop(stack_t *stack, entry_t *e);

static inline void push(stack_t *stack, int i, uint64_t k, uint64_t l, int n_mm, int is_diff);
static inline void shadow(int x, int len, uint64_t max, int last_diff_pos, bwt_width_t *w);

static int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, uint64_t *k0, uint64_t *l0);
static uint64_t* get_pos_from_sa_interval(const bwaidx_t* idx, const uint64_t k, const uint64_t l, input_query* info_seq, uint64_t* numer_forward,
											uint64_t* numer_rc);

static void get_string_to_pos(const bwaidx_t *idx, input_query* info_seq, search_result* sr);
static void add_entry_to_result(const bwaidx_t* idx, search_result** rs, uint8_t n_mm, uint64_t n_f, uint64_t n_revC, uint64_t* pos,
								input_query* info_seq, int* aln_n);

static int comp(const void * elem1, const void * elem2);

//Only for debug
static void print_pattern_to_search(const ubyte_t * seq, int len);
static void print_character(int i);

//Modified version of bwt_match_gap in bwtgap.c file
search_result** get_approximate_match(const bwaidx_t* idx, input_query* info_seq, uint32_t* num_hits)
{
	//Warning: Hit means an admissible match (with respect to the input parameters) between the pattern and the reference

	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr,"Start get_approximate_match\n");

	/*
	 * Number of permissible match found.
	 * Warning: if the current hit found and its r.v. are present inside the reference, it's increased by two
	 */
	int n_hit_found = 0;

	//Initially, the structure containing the info relative to the hits found has size 4.
	int m_hit_found = 4;
	search_result** returnM = (search_result**) malloc(m_hit_found * sizeof(search_result*));
	for (int i = 0; i < m_hit_found; i++)
		*(returnM + i) = (search_result *) malloc(sizeof(search_result));

	/*
	 * Heap-like data structure to keep partial hits. It is prioritized on the number of mismatch
	 * inside the partial hits(also called entry): less mismatch have a entry first is extract
	 */
	stack_t *stack = init_stack(info_seq->max_diff);

	//see cal_width implementation
	bwt_width_t* width = (bwt_width_t*) malloc((info_seq->len + 1) * sizeof(bwt_width_t));
	cal_width(idx->bwt, info_seq->len, info_seq->seq, width);

	//Get bwt 
	bwt_t* bwt = idx->bwt;

	//Max number(s) of mismatch(es) between a hit and the reference
	const int max_diff = info_seq->max_diff;

	//Print pattern to search
	if (verbose_bound_backtracking_search > 2)
	{
		fprintf(stderr,"\nThe search pattern is the r.c. of input pattern: ");
		print_pattern_to_search(info_seq->seq, info_seq->len);
	}
	//Reset stack
	reset_stack(stack);

	push(stack, info_seq->len, 0, bwt->seq_len, 0, 0);

	//With priority on the number of mismatches, partial hits calculated previously are extracted
	while (stack->n_entries)
	{
		//Get the best entry among the previously calculated partial hits
		entry_t e;
		pop(stack, &e);

		//Defines whether the current entry 'e' is a hit
		bool hit_found = false;

		/*
		 * The number of mismatch of the current entry regard the substring [i,n-1] of the
		 * input pattern.
		 */
		int i = e.info;

		//(k,l) is the SA region of [i,n-1]
		uint64_t k = e.k;
		uint64_t l = e.l;

		/*
		 * For the current entry 'e' previously calculated, n_mm_yet_allow is the max number of
		 * mismatch that 'e' can still contain
		 */
		int n_mm_yet_allow = max_diff - e.n_mm;

		if (verbose_bound_backtracking_search > 2)
		{
			fprintf(stderr,"\n*New entry is extracted, this has: ");
			fprintf(stderr,"\ni: %d\n", i);
			fprintf(stderr,"k: %ld\n", k);
			fprintf(stderr,"l: %ld\n", l);
			fprintf(stderr,"n_mm_yet_allow: %d\n\n", n_mm_yet_allow);
		}

		//This control should be useless..
		if (n_mm_yet_allow < 0)
			continue;

		if (i > 0 && n_mm_yet_allow < width[i - 1].bid)
		{
			if (verbose_bound_backtracking_search > 2)
			{
				fprintf(stderr,"°°Used width!°°\n");
				fprintf(stderr,"n_mm_yet_allow: %d\n", n_mm_yet_allow);
				fprintf(stderr,"width[i - 1].bid: %d\n", width[i - 1].bid);
			}
			continue;
		}

		// Check whether a hit is found
		if (i == 0)
		{
			//This means that the length of the current entry is exactly equal to the length of the input patter
			hit_found = true;
		}
		else if (n_mm_yet_allow == 0)
		{
			/*
			 * This means that the current partial hit have the max possible number of mismatch,
			 * so the only possible hit is composed by 'e' and the miss base must be equal to the
			 * base of the input patter
			 */
			if (verbose_bound_backtracking_search > 2)
				fprintf(stderr,"bwt_match_exact_alt\n");

			if (bwt_match_exact_alt(bwt, i, info_seq->seq, &k, &l))
			{
				if (verbose_bound_backtracking_search > 2)
				{
					fprintf(stderr,"hit found\n");
					fprintf(stderr,"k: %ld\n", k);
					fprintf(stderr,"l: %ld\n", l);
				}
				hit_found = true;
			}
			else
			{
				if (verbose_bound_backtracking_search > 2)
					fprintf(stderr,"Hit not found from bwt_match_exact_alt\n");

				continue; // no hit, skip
			}
		}
		if (hit_found)
		{
			/*
			 * Inside the suffix interval [k,l], some position is relative
			 * of a hit that is forward respect the pattern input, other
			 * is r.c. respect the pattern input
			 */
			uint64_t numer_forward = 0;
			uint64_t numer_revC = 0;

			//This method is used to reduce the search space
			shadow(l - k + 1, info_seq->len, bwt->seq_len, e.last_diff_pos, width);

			//We might have to add both the direct and the reverse complement of the current entry
			if (!(n_hit_found + 2 <= m_hit_found))
			{
				m_hit_found <<= 1;
				(returnM) = (search_result**) realloc(returnM, m_hit_found * sizeof(search_result*));
				for (int i = n_hit_found; i < m_hit_found; i++)
					*(returnM + i) = (search_result *) malloc(sizeof(search_result));

				if (verbose_bound_backtracking_search > 2)
					fprintf(stderr,"Realloc successfully completed\n");
			}
			if (verbose_bound_backtracking_search > 2)
			{
				printf("k: %"PRIu64 "\n", k);
				printf("l: %"PRIu64 "\n", l);
			}

			//For each index j that belongs to the suffix interval [k,l], get_pos_from_sa_interval
			//return all SA(j) valid
			uint64_t* position_to_ref = get_pos_from_sa_interval(idx, k, l, info_seq, &numer_forward, &numer_revC);

			//If all positions corresponding to suffix interval [k, l] are not valid continue
			if ((numer_forward + numer_revC) == 0)
				continue;

			//Add hit(s) in returnM and update n_hit_found
			add_entry_to_result(idx, returnM, e.n_mm, numer_forward, numer_revC, position_to_ref, info_seq, &n_hit_found);

			free(position_to_ref);

			if (info_seq->type_search == ARBITRARY_HIT)
				break;

			continue;
		}

		//Increase(not decrease..) the size of the current partial hit
		--i;

		/*
		 * Give k - 1 and l of the current entry 'e', compute
		 * O(x,k-1) and O(x,l), where x can be A,C,G,T.
		 *
		 * For more detail please see section (2.3) and (2.4) of original paper:
		 * https://academic.oup.com/bioinformatics/article/25/14/1754/225615
		 *
		 */
		uint64_t cnt_k[4], cnt_l[4];
		bwt_2occ4(bwt, k - 1, l, cnt_k, cnt_l);

		//Try to extend the current partial hit whit every possible base
		for (int j = 1; j <= 4; ++j)
		{
			int c = (info_seq->seq[i] + j) & 3;

			if (verbose_bound_backtracking_search > 2)
			{
				fprintf(stderr,">Try to extend the current partial hit whit: ");
				print_character(c);
			}
			int is_mm = (j != 4 || info_seq->seq[i] > 3);
			k = bwt->L2[c] + cnt_k[c] + 1;
			l = bwt->L2[c] + cnt_l[c];

			if (k <= l)
			{
				if (verbose_bound_backtracking_search > 2)
				{
					if (is_mm)
						fprintf(stderr, "Found(whit mismatch)-->push inside stack this partial hit whit this values: \n");
					else
						fprintf(stderr, "Found(no mismatch)-->push inside stack this partial hit whit this values: \n");

					fprintf(stderr,"i %d\n", i);
					fprintf(stderr,"k %ld:\n", k);
					fprintf(stderr,"l %ld:\n", l);
					fprintf(stderr,"num. mismatch inside this partial hit is: %d\n", e.n_mm + is_mm);
				}
				//The partial hit is add in the reference
				push(stack, i, k, l, e.n_mm + is_mm, is_mm);
			}
			else
			{
				if (verbose_bound_backtracking_search > 2)
					fprintf(stderr, "Not found-->discard\n");
			}
		}
	}

	*(num_hits) = n_hit_found;

	if (verbose_bound_backtracking_search > 2)
	{
		fprintf(stderr, "\nNumber of hit found: %d\n", n_hit_found);
		fprintf(stderr, "End get_approximate_match\n");
	}

	//free memory
	destroy_stack(stack);
	free(width);

	return returnM;
}

input_query* init_bwt_seq(uint8_t* pattern_to_search, const size_t pattern_len, const uint8_t nmismatch, const uint8_t type_search,
							const uint8_t type_output)
{
	input_query *seq = (input_query*) calloc(1, sizeof(input_query));
	seq->seq = pattern_to_search;
	seq->len = pattern_len;
	seq->max_diff = nmismatch;
	seq->type_search = type_search;
	seq->type_output = type_output;
	return seq;
}

uint8_t* create_pattern_to_search(const char* pattern_input, const size_t pattern_len)
{
	//Code pattern input
	uint8_t* pattern_to_search = (uint8_t*) malloc(pattern_len * sizeof(int8_t));

	uint8_t base_to_int[] =

	{0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0xFF, 0x01, 0xFF, 0xFF,
		0xFF, 0x02, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x04, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0xFF, 0x01, 0xFF, 0xFF, 0xFF, 0x02, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x04, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF};

	for (int j = 0; j < pattern_len; ++j)
	{
		pattern_to_search[j] = base_to_int[(size_t) pattern_input[pattern_len - j - 1]]; //reverse
		pattern_to_search[j] = pattern_to_search[j] > 3 ? 4 : 3 - pattern_to_search[j]; //complement
	}
	return pattern_to_search;
}

/*
 *Warning: this method works on the reverse (not complement) of the patter T given by the user.
 *
 *width[i].bid represents the lower bound of the number of differences in P_r[0,i], where
 *P_r is the reverse of the input pattern. This is useful for making the search space smaller.
 */
//Derived from bwt_cal_width in bwtaln.c file
static int cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr, "\n->Inside bwt_cal_width\n");

	//[k,l] is the S.A interval values
	uint64_t k, l;

	//ok e ol represent represent the values of the rank function
	uint64_t ok, ol;

	int i, bid;
	bid = 0;

	k = 0;
	l = bwt->seq_len;

	for (i = 0; i < len; ++i)
	{
		//Before i have complement
		ubyte_t c = 3 - str[i];

		if (verbose_bound_backtracking_search > 2)
			fprintf(stderr, "Current character considered: %c\n", "ACGTN"[c]);

		if (c < 4)
		{
			//This method calculates the rank function(usually indicated with O) of c
			bwt_2occ(bwt, k - 1, l, c, &ok, &ol);

			//L2[c] is the Count function (usually indicated with C)
			k = bwt->L2[c] + ok + 1;
			l = bwt->L2[c] + ol;

			if (verbose_bound_backtracking_search > 2)
			{
				fprintf(stderr, "k: %ld\n", k);
				fprintf(stderr, "l: %ld\n", l);
			}
		}
		if (k > l || c > 3)
		{
			// then restart
			if (verbose_bound_backtracking_search > 2)
				fprintf(stderr, "Restart\n");

			k = 0;
			l = bwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;

		if (verbose_bound_backtracking_search > 2)
		{
			fprintf(stderr, "--->width[%d].w = %ld\n", i, l - k + 1);
			fprintf(stderr, "--->width[%d].bid = %d\n", i, bid);
		}
	}

	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

//It's the same of gap_reset_stack in bwtgap.c file
static void reset_stack(stack_t *stack)
{
	int i;
	for (i = 0; i != stack->n_sub_stacks; ++i)
		stack->stacks[i].n_entries_substack = 0;
	stack->best = stack->n_sub_stacks;
	stack->n_entries = 0;
}

//Derived from init_stack2 in bwtgap.c file
static stack_t * init_stack(int nmismatch)
{
	stack_t *stack = (stack_t*) calloc(1, sizeof(stack_t));

	//Each partial hit, based on the number of mismatch it contains, will be "clustered together",
	//in substack(so we can have at most nmismatch + 1 substack)
	stack->n_sub_stacks = nmismatch + 1;
	stack->stacks = (substack_t*) calloc(stack->n_sub_stacks, sizeof(stack_t));
	return stack;
}

//Derived from gap_push in bwtgap.c file
static inline void push(stack_t *stack, int i, uint64_t k, uint64_t l, int n_mm, int is_diff)
{
	//Get pointer to substack relative to partial hit whith n_mm mismatch
	substack_t *q = stack->stacks + n_mm;

	if (q->n_entries_substack == q->m_entries)
	{
		q->m_entries = q->m_entries ? q->m_entries << 1 : 4;
		q->start_substack = (entry_t*) realloc(q->start_substack, sizeof(entry_t) * q->m_entries);
	}

	//Add entry at the end of the substack
	entry_t *p = q->start_substack + q->n_entries_substack;
	p->info = (uint32_t) i;
	p->k = k;
	p->l = l;
	p->n_mm = n_mm;
	p->last_diff_pos = is_diff ? i : 0;

	//Increase the total number of entries whit n_mm mismatches
	++(q->n_entries_substack);

	//Increase the total number of entries in general
	++(stack->n_entries);

	if (stack->best > n_mm)
	{
		if (verbose_bound_backtracking_search > 2)
			fprintf(stderr, "Update stack-> best\n");
		stack->best = n_mm;
	}
}

//Derived from gap_pop in bwtgap.c file
static inline void pop(stack_t *stack, entry_t *e)
{
	/*
	 * Substack containing partial hits whit number of mismatch is minimum,
	 * all partial hit that NOT belong to this memory region have more mismatch.
	 */
	substack_t *q = stack->stacks + stack->best;

	//Of this substack, we take the last one memorized
	*e = q->start_substack[q->n_entries_substack - 1];

	//Decrease the total number of entries whit n_mm mismatches
	--(q->n_entries_substack);
	//Decrease the total number of entries in general
	--(stack->n_entries);

	/*
	 * If the number of partial hit whit n_mm mismatches now is zero but
	 * the number of entry is not zero, update best for next pop operation
	 */
	if (q->n_entries_substack == 0 && stack->n_entries)
	{
		// reset best
		int i;
		for (i = stack->best + 1; i < stack->n_sub_stacks; ++i)
			if (stack->stacks[i].n_entries_substack != 0)
				break;
		stack->best = i;
	}
	else if (stack->n_entries == 0)
	{
		//Necessary for the first iteration
		stack->best = stack->n_sub_stacks;
	}
}

//It's the same of original function gap_shadow present in bwtgap.c
static inline void shadow(int x, int len, uint64_t max, int last_diff_pos, bwt_width_t *w)
{
	if (verbose_bound_backtracking_search > 2)
	{
		fprintf(stderr, "->Inside shadow\n");
		fprintf(stderr, "last_diff_pos: %d\n", last_diff_pos);
	}
	int i, j;
	for (i = j = 0; i < last_diff_pos; ++i)
	{
		if (verbose_bound_backtracking_search > 2)
		{
			fprintf(stderr, "Current values:\n");
			fprintf(stderr, "w[%d].w: %ld\n", i, w[i].w);
			fprintf(stderr, "w[%d].bid: %d\n", i, w[i].bid);
			fprintf(stderr, "x = l-k+1: %d\n", x);
		}
		if (w[i].w > x)
		{
			w[i].w -= x;
			if (verbose_bound_backtracking_search > 2)
				fprintf(stderr, "Update: w[%d] = %ld: ", i, w[i].w);
		}
		else if (w[i].w == x)
		{
			w[i].bid = 1;
			w[i].w = max - (++j);
			if (verbose_bound_backtracking_search > 2)
			{
				fprintf(stderr, "\nUpdate bid (1)\n");
				fprintf(stderr, "Update w[%d].w: %ld\n", i, w[i].w);
			}
		} // else should not happen
	}

	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr, "->Exit gap_shadow\n");
}

//It's the same of original, present in bwt.c
static int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, uint64_t *k0, uint64_t *l0)
{
	int i;
	uint64_t k, l, ok, ol;
	k = *k0;
	l = *l0;
	for (i = len - 1; i >= 0; --i)
	{
		ubyte_t c = str[i];
		if (c > 3)
			return 0; // there is an N here. no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l)
			return 0; // no match
	}
	*k0 = k;
	*l0 = l;
	return l - k + 1;
}

static char* convert_hit_to_string(const ubyte_t * seq, int len)
{
	char* hit_pattern = calloc(len + 1, sizeof(char));
	for (int j = 0; j < len; j++)
		hit_pattern[j] = "ACGTN"[seq[j]];

	hit_pattern[len] = '\0';
	return hit_pattern;
}

static uint64_t* get_pos_from_sa_interval(const bwaidx_t* idx, uint64_t k, uint64_t l, input_query* info_seq, uint64_t* numer_forward,
											uint64_t* numer_rc)
{
	//Current index of the SA, k<=t<=l
	uint64_t t;

	//Is SA[t]
	uint64_t pos;

	//Define if pos is the start position of the input patter inside the reference genome or the r.c of the input pattern
	bool strand;

	//Max number of occurrences of current hit found inside the reference
	uint64_t max_occ = l - k + 1;

	//Same meaning of numer_forward and numer_rc
	uint64_t n_f = 0;
	uint64_t n_rc = 0;

	/*
	 * Store all valid positions. The first n_f positions are relative to the forward of the current hit,
	 * the last n_rc to the r.c.
	 */
	uint64_t* set_pos = (uint64_t*) malloc(max_occ * sizeof(uint64_t));

	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr, "Max_occ: %" PRIu64 "\n", max_occ);

	for (t = k; t <= l; ++t)
	{
		//Calculate SA[t](base 1)
		pos = bwa_sa2pos(idx->bns, idx->bwt, t, info_seq->len, &strand);

		if (verbose_bound_backtracking_search > 2)
		{
			if (pos == ULLONG_MAX)
				fprintf(stderr, "Position return by bwa_sa2pos is not valid\n");
			else
			{
				fprintf(stderr, "Position return by bwa_sa2pos: %" PRIu64 "\n", pos);
				fprintf(stderr, "Strand: %" PRIu8 "\n", strand);
			}
		}
		//Check if pos is admissible
		if (pos != ULLONG_MAX)
		{
			if (strand == 0)
			{
				set_pos[n_f] = pos;
				n_f++;
			}
			else if ((strand == 1) && (info_seq->type_output == ALLOW_REV_COMP))
			{
				set_pos[max_occ - n_rc - 1] = pos;
				n_rc++;
			}
		}
	}

	//Number of valid positions found
	uint64_t realSize = n_f + n_rc;

	if (verbose_bound_backtracking_search > 2)
	{
		fprintf(stderr, "\nn_rc: %"PRIu64 "\n", n_rc);
		fprintf(stderr, "n_f: %" PRIu64"\n", n_f);
	}
	if ((n_rc == n_f) && (n_f == 0))
	{
		if (verbose_bound_backtracking_search > 2)
			fprintf(stderr, "No hits have been added\n");
		return 0;
	}
	//Compact set_pos
	if (n_f + n_rc < max_occ)
	{
		for (int i = 0; i < max_occ - realSize; i++)
			set_pos[n_f + i] = set_pos[max_occ - i - 1];
		set_pos = realloc(set_pos, realSize * sizeof(uint64_t));
	}
	if (verbose_bound_backtracking_search > 2)
	{
		fprintf(stderr, "The positions store inside set_pos are: \n");
		for (int i = 0; i < realSize; i++)
			printf("%" PRIu64 "\n", set_pos[i]);
	}

	*(numer_forward) = n_f;
	*(numer_rc) = n_rc;

	return set_pos;
}

//Same of gap_destroy_stack in bwtgap.c
static void destroy_stack(stack_t *stack)
{
	int i;
	for (i = 0; i != stack->n_sub_stacks; ++i)
		free(stack->stacks[i].start_substack);
	free(stack->stacks);
	free(stack);
}

static void get_string_to_pos(const bwaidx_t *idx, input_query* info_seq, search_result* sr)
{
	int64_t len = 0;

	if (verbose_bound_backtracking_search > 2)
	{
		fprintf(stderr, "\nStart get_string_to_pos\n");
		fprintf(stderr, "beg: %"PRIu64 "\n", sr->positions_to_ref[0]);
		fprintf(stderr, "end: %"PRIu64 "\n", sr->positions_to_ref[0] + info_seq->len);
	}

	//Get string that start in position position_to_ref[0]
	uint8_t* hit_code = bns_get_seq(idx->bns->l_pac, idx->pac, sr->positions_to_ref[0], sr->positions_to_ref[0] + info_seq->len, &len);

	if (sr->n_mismatches > 0)
	{
		size_t num_different_pos_found = 0;
		sr->different_positions = calloc(sr->n_mismatches, sizeof(uint32_t));

		if (info_seq->len == len) //this must always be true..
		{
			//Check where the mismatches are
			for (uint32_t j = 0; j < info_seq->len; ++j)
			{
				if ((sr->is_rev_comp == true) && (3 - hit_code[info_seq->len - j - 1] != 3 - info_seq->seq[info_seq->len - j - 1]))
				{
					sr->different_positions[num_different_pos_found] = j + 1; //base-1
					num_different_pos_found++;
				}
				else if ((sr->is_rev_comp == false) && (hit_code[j] != 3 - info_seq->seq[info_seq->len - j - 1]))
				{
					sr->different_positions[num_different_pos_found] = j + 1; //base-1
					num_different_pos_found++;
				}
			}
		}
		else //this must always be false..
		{
			printf("Impossible to extract the patter from the reference");
			exit(EXIT_FAILURE);
		}
	}
	else
		sr->different_positions = NULL;

	sr->hit = convert_hit_to_string(hit_code, len);
}

//In the next version of software create a method for add a entry into rs.
static void add_entry_to_result(const bwaidx_t* idx, search_result** rs, uint8_t n_mm, uint64_t n_f, uint64_t n_revC, uint64_t* pos,
								input_query* info_seq, int* aln_n)
{
	if (verbose_bound_backtracking_search > 2)
		fprintf(stderr, "-> Inside add_entry_to_result\n");

	if (n_f > 0)
	{
		if (verbose_bound_backtracking_search > 2)
			fprintf(stderr, "-Forward\n");

		rs[*(aln_n)]->positions_to_ref = malloc(n_f * sizeof(uint64_t));

		for (int i = 0; i < n_f; i++)
			rs[*(aln_n)]->positions_to_ref[i] = pos[i];

		qsort(rs[*(aln_n)]->positions_to_ref, n_f, sizeof(*rs[*(aln_n)]->positions_to_ref), comp);
		rs[*(aln_n)]->is_rev_comp = false;
		rs[*(aln_n)]->num_occur = n_f;
		rs[*(aln_n)]->n_mismatches = n_mm;
		get_string_to_pos(idx, info_seq, rs[*(aln_n)]);
		*(aln_n) = *(aln_n) + 1;

		if (info_seq->type_search == ARBITRARY_HIT)
			return;
	}

	if (n_revC > 0)
	{
		if (verbose_bound_backtracking_search > 2)
			fprintf(stderr, "-Reverse complement\n");

		rs[*(aln_n)]->positions_to_ref = malloc(n_revC * sizeof(uint64_t));

		for (int i = n_f; i < n_f + n_revC; i++)
			rs[*(aln_n)]->positions_to_ref[i - n_f] = pos[i];

		qsort(rs[*(aln_n)]->positions_to_ref, n_revC, sizeof(*rs[*(aln_n)]->positions_to_ref), comp);
		rs[*(aln_n)]->is_rev_comp = true;
		rs[*(aln_n)]->num_occur = n_revC;
		rs[*(aln_n)]->n_mismatches = n_mm;
		get_string_to_pos(idx, info_seq, rs[*(aln_n)]);
		*(aln_n) = *(aln_n) + 1;
	}
}

static void print_pattern_to_search(const ubyte_t * seq, int len)
{
	for (int j = 0; j < len; j++)
		fprintf(stderr, "%c", "ACGTN"[seq[j]]);

	fprintf(stderr, "\n");
}

static void print_character(int i)
{
	char t;

	switch (i)
	{
		case 0:
			t = 'A';
			break;
		case 1:
			t = 'C';
			break;
		case 2:
			t = 'G';
			break;
		case 3:
			t = 'T';
			break;

	}
	fprintf(stderr, "%c\n", t);
}

static int comp(const void * elem1, const void * elem2)
{
	if (*(uint64_t*) elem1 == *(uint64_t*) elem2)
		return 0;
	return *(uint64_t*) elem1 < *(uint64_t*) elem2 ? -1 : 1;
}
//To improve
void controlParam(const bwaidx_t *idx, const char* pattern_input, const uint8_t type_search, const uint8_t type_output)
{
	if (idx == 0)
	{
		fprintf(stderr, "Miss index.\n");
		exit(EXIT_FAILURE);
	}
	else if (pattern_input == 0)
	{
		fprintf(stderr, "Miss pattern to search.\n");
		exit(EXIT_FAILURE);
	}
	else if ((type_search != ARBITRARY_HIT) && (type_search != ALL_HITS))
	{
		fprintf(stderr, "type_search not legal.\n");
		exit(EXIT_FAILURE);
	}
	else if ((type_output != ALLOW_REV_COMP) && (type_output != NO_REV_COMP))
	{
		fprintf(stderr, "type_search not legal.\n");
		exit(EXIT_FAILURE);
	}
}
