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
#ifndef OCC_H
#define OCC_H

#include <inttypes.h>

#ifndef BWA_UBYTE
#define BWA_UBYTE

typedef unsigned char ubyte_t;

#endif

#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))

// requirement: (OCC_INTERVAL%16 == 0); please DO NOT change this line because some part of the code assume OCC_INTERVAL=0x80
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

#ifndef BWT_T
typedef struct
{
	uint64_t primary; // S^{-1}(0), or the primary index of BWT
	uint64_t L2[5]; // C(), cumulative count
	uint64_t seq_len; // sequence length
	uint64_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	// occurance array, separated to two parts
	uint32_t cnt_table[256];
	// suffix array
	int sa_intv;
	uint64_t n_sa;
	uint64_t *sa;

} bwt_t;
#define BWT_T
#endif

#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))

#define __occ_aux4(bwt, b)											\
	((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
	 + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

#ifdef __cplusplus
extern "C"
{
#endif
	uint64_t bwt_occ(const bwt_t *bwt, uint64_t k, ubyte_t c);
	void bwt_2occ(const bwt_t *bwt, uint64_t k, uint64_t l, ubyte_t c, uint64_t *ok, uint64_t *ol);
	void bwt_2occ4(const bwt_t *bwt, uint64_t k, uint64_t l, uint64_t cntk[4], uint64_t cntl[4]);
#ifdef __cplusplus
}
#endif

#endif
