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
#ifndef AC_KVEC_H
#define AC_KVEC_H

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "utils.h"

#ifndef BNTANN1_T
#define BNTANN1_T

typedef struct
{
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	int32_t is_alt;
	char *name, *anno;
} bntann1_t;

#endif

#ifndef BNTAMB1_T
#define BNTAMB1_T

typedef struct
{
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

#endif

#ifndef BNTSEQ_T
#define BNTSEQ_T

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

#endif

#ifndef BWT_T
#define BWT_T

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

#endif

#ifdef __cplusplus
extern "C"
{
#endif

	char *bwa_idx_infer_prefix(const char *hint);
	bwt_t *bwa_idx_load_bwt(const char *hint);
	bntseq_t *bns_restore(const char *prefix);
	bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename);
	bwt_t *bwt_restore_bwt(const char *fn);
	void bwt_destroy(bwt_t *bwt);
#ifdef __cplusplus
}
#endif
#endif
