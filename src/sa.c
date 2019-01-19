/* The MIT License

 Copyright (c) 2008 Genome Research Ltd (GRL).

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

#include "sa.h"
#include "kvec.h"

static inline uint64_t bwt_invPsi(const bwt_t *bwt, uint64_t k);

/*
 * It's the same of the original method present in bwt.c
 * This is used both in backtracking_search.c and fmdindex_build.c
 */
uint64_t bwt_sa(const bwt_t *bwt, uint64_t k)
{
	uint64_t sa = 0, mask = bwt->sa_intv - 1;
	while (k & mask)
	{
		++sa;
		k = bwt_invPsi(bwt, k);
	}

	/* without setting bwt->sa[0] = -1, the following line should be
	 changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k / bwt->sa_intv];
}

/*
 * Same of bwt_cal_sa present in bwt.c.
 * This method is used both in fmdindex_load.c and fmdindex_build.c
 * NOTE: bwt->bwt and bwt->occ must be precalculated
 */
void bwt_cal_sa(bwt_t *bwt, int intv)
{
	uint64_t isa, sa, i; // S(isa) = sa
	int intv_round = intv;

	kv_roundup32(intv_round);
	xassert(intv_round == intv, "SA sample interval is not a power of 2.");
	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa)
		free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (uint64_t*) calloc(bwt->n_sa, sizeof(uint64_t));
	// calculate SA value
	isa = 0;
	sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i)
	{
		if (isa % intv == 0)
			bwt->sa[isa / intv] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if (isa % intv == 0)
		bwt->sa[isa / intv] = sa;
	bwt->sa[0] = (uint64_t) -1; // before this line, bwt->sa[0] = bwt->seq_len
}

//It's the same method present in bwt.c
static inline uint64_t bwt_invPsi(const bwt_t *bwt, uint64_t k) // compute inverse CSA
{
	uint64_t x = k - (k > bwt->primary);
	x = bwt_B0(bwt, x);
	x = bwt->L2[x] + bwt_occ(bwt, k, x);
	return k == bwt->primary ? 0 : x;
}
