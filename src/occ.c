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
#include <stdio.h>
#include <string.h>
#include "occ.h"

static void bwt_occ4(const bwt_t *bwt, uint64_t k, uint64_t cnt[4]);
static inline int __occ_aux(uint64_t y, int c);

/*
 * Same of original method present in bwt.c.
 * This is used both in sa.c and backtracking_search.c files
 */
uint64_t bwt_occ(const bwt_t *bwt, uint64_t k, ubyte_t c)
{
	uint64_t n;
	uint32_t *p, *end;

	if (k == bwt->seq_len)
		return bwt->L2[c + 1] - bwt->L2[c];
	if (k == (uint64_t) (-1))
		return 0;
	k -= (k >= bwt->primary); // because $ is not in bwt

	// retrieve Occ at k/OCC_INTERVAL
	n = ((uint64_t*) (p = bwt_occ_intv(bwt, k)))[c];
	p += sizeof(uint64_t); // jump to the start of the first BWT cell

	// calculate Occ up to the last k/32
	end = p + (((k >> 5) - ((k & ~OCC_INTV_MASK) >> 5)) << 1);
	for (; p < end; p += 2)
		n += __occ_aux((uint64_t) p[0] << 32 | p[1], c);

	// calculate Occ
	n += __occ_aux(((uint64_t) p[0] << 32 | p[1]) & ~((1ull << ((~k & 31) << 1)) - 1), c);
	if (c == 0)
		n -= ~k & 31; // corrected for the masked bits

	return n;
}

/*
 * It's the same present in bwt.c file
 * NOTE: an analogy to bwt_occ() but more efficient, requiring k <= l
 */
void bwt_2occ(const bwt_t *bwt, uint64_t k, uint64_t l, ubyte_t c, uint64_t *ok, uint64_t *ol)
{
	uint64_t _k, _l;
	_k = (k >= bwt->primary) ? k - 1 : k;
	_l = (l >= bwt->primary) ? l - 1 : l;
	if (_l / OCC_INTERVAL != _k / OCC_INTERVAL || k == (uint64_t) (-1) || l == (uint64_t) (-1))
	{
		*ok = bwt_occ(bwt, k, c);
		*ol = bwt_occ(bwt, l, c);
	}
	else
	{
		uint64_t m, n, i, j;
		uint32_t *p;
		if (k >= bwt->primary)
			--k;
		if (l >= bwt->primary)
			--l;
		n = ((uint64_t*) (p = bwt_occ_intv(bwt, k)))[c];
		p += sizeof(uint64_t);
		// calculate *ok
		j = k >> 5 << 5;
		for (i = k / OCC_INTERVAL * OCC_INTERVAL; i < j; i += 32, p += 2)
			n += __occ_aux((uint64_t) p[0] << 32 | p[1], c);
		m = n;
		n += __occ_aux(((uint64_t) p[0] << 32 | p[1]) & ~((1ull << ((~k & 31) << 1)) - 1), c);
		if (c == 0)
			n -= ~k & 31; // corrected for the masked bits
		*ok = n;
		// calculate *ol
		j = l >> 5 << 5;
		for (; i < j; i += 32, p += 2)
			m += __occ_aux((uint64_t) p[0] << 32 | p[1], c);
		m += __occ_aux(((uint64_t) p[0] << 32 | p[1]) & ~((1ull << ((~l & 31) << 1)) - 1), c);
		if (c == 0)
			m -= ~l & 31; // corrected for the masked bits
		*ol = m;
	}
}


/*
 * It's the same of bwt_2occ4 present in bwt.c file
 * NOTE:an analogy to bwt_occ4() but more efficient, requiring k <= l
 */
void bwt_2occ4(const bwt_t *bwt, uint64_t k, uint64_t l, uint64_t cntk[4], uint64_t cntl[4])
{
	uint64_t _k, _l;
	_k = k - (k >= bwt->primary);
	_l = l - (l >= bwt->primary);
	if (_l >> OCC_INTV_SHIFT != _k >> OCC_INTV_SHIFT || k == (uint64_t) (-1)
	|| l == (uint64_t) (-1))
	{
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	}
	else
	{
		uint64_t x, y;
		uint32_t *p, tmp, *endk, *endl;
		k -= (k >= bwt->primary); // because $ is not in bwt
		l -= (l >= bwt->primary);
		p = bwt_occ_intv(bwt, k);
		memcpy(cntk, p, 4 * sizeof(uint64_t));
		p += sizeof(uint64_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
		// prepare cntk[]
		endk = p + ((k >> 4) - ((k & ~OCC_INTV_MASK) >> 4));
		endl = p + ((l >> 4) - ((l & ~OCC_INTV_MASK) >> 4));
		for (x = 0; p < endk; ++p)
			x += __occ_aux4(bwt, *p);
		y = x;
		tmp = *p & ~((1U << ((~k & 15) << 1)) - 1);
		x += __occ_aux4(bwt, tmp) - (~k & 15);
		// calculate cntl[] and finalize cntk[]
		for (; p < endl; ++p)
			y += __occ_aux4(bwt, *p);
		tmp = *p & ~((1U << ((~l & 15) << 1)) - 1);
		y += __occ_aux4(bwt, tmp) - (~l & 15);
		memcpy(cntl, cntk, 4 * sizeof(uint64_t));
		cntk[0] += x & 0xff;
		cntk[1] += x >> 8 & 0xff;
		cntk[2] += x >> 16 & 0xff;
		cntk[3] += x >> 24;
		cntl[0] += y & 0xff;
		cntl[1] += y >> 8 & 0xff;
		cntl[2] += y >> 16 & 0xff;
		cntl[3] += y >> 24;
	}
}
//It's the same method present in bwt.c file.
static void bwt_occ4(const bwt_t *bwt, uint64_t k, uint64_t cnt[4])
{
	uint64_t x;
	uint32_t *p, tmp, *end;
	if (k == (uint64_t) (-1))
	{
		memset(cnt, 0, 4 * sizeof(uint64_t));
		return;
	}
	k -= (k >= bwt->primary); // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	memcpy(cnt, p, 4 * sizeof(uint64_t));
	p += sizeof(uint64_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
	end = p + ((k >> 4) - ((k & ~OCC_INTV_MASK) >> 4)); // this is the end point of the following loop
	for (x = 0; p < end; ++p)
		x += __occ_aux4(bwt, *p);
	tmp = *p & ~((1U << ((~k & 15) << 1)) - 1);
	x += __occ_aux4(bwt, tmp) - (~k & 15);
	cnt[0] += x & 0xff;
	cnt[1] += x >> 8 & 0xff;
	cnt[2] += x >> 16 & 0xff;
	cnt[3] += x >> 24;
}

//It's the same method present in bwt.c file.
static inline int __occ_aux(uint64_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = ((c & 2) ? y : ~y) >> 1 & ((c & 1) ? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}
