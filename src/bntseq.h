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

#ifndef BWT_BNTSEQ_H
#define BWT_BNTSEQ_H

#include <stdbool.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include "sa.h"

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef uint8_t ubyte_t;
#endif

#ifndef BNTANN1_T
typedef struct
{
		int64_t offset;
		int32_t len;
		int32_t n_ambs;
		uint32_t gi;
		int32_t is_alt;
		char *name, *anno;
} bntann1_t;
#define BNTANN1_T
#endif

#ifndef BNTAMB1_T
typedef struct
{
		int64_t offset;
		int32_t len;
		char amb;

} bntamb1_t;
#define BNTAMB1_T
#endif


#ifndef BNTSEQ_T
typedef struct
{
		int64_t l_pac;
		int32_t n_seqs;
		uint32_t seed;
		bntann1_t *anns; // n_seqs elements
		int32_t n_holes;
		bntamb1_t *ambs; // n_holes elements
		FILE *fp_pac;
} bntseq_t;
#define BNTSEQ_T
#endif

extern unsigned char nst_nt4_table[256];

#ifdef __cplusplus
extern "C" {
#endif

	uint64_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, uint64_t sapos, int ref_len, bool *strand);
	int64_t bns_fasta2bntseq(gzFile fp_fa, const char *prefix, int for_only);
	uint8_t *bns_get_seq(int64_t l_pac, const uint8_t *pac, int64_t beg, int64_t end, int64_t *len);

#ifdef __cplusplus
}
#endif

//It's the same method present in bntseq.h in original version
static inline int64_t bns_depos(const bntseq_t *bns, int64_t pos, bool *is_rev)
{
	return (*is_rev = (pos >= bns->l_pac))? (bns->l_pac<<1) - 1 - pos : pos;
}

#endif
