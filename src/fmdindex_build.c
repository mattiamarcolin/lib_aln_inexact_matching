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

#include <string.h>
#include "fmdindex_load.h"
#include "sa.h"
#include "bntseq.h"
#include "rle.h"
#include "rope.h"
#include "utils.h"

#ifdef _DIVBWT
#include "divsufsort.h"
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

int is_bwt(ubyte_t *T, int n);
void bwt_bwtgen2(const char *fn_pac, const char *fn_bwt, int block_size);

static void bwt_dump_bwt(const char *fn, const bwt_t *bwt);
static void bwt_dump_sa(const char *fn, const bwt_t *bwt);
static int64_t bwa_seq_len(const char *fn_pac);
static bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is);
static int bwa_pac2bwt(int argc, char *argv[]);
static void bwt_bwtupdate_core(bwt_t *bwt);
static int bwa_bwtupdate(int argc, char *argv[]);
static int bwa_bwt2sa(int argc, char *argv[]);

int fmd_idx_build(const char *fa, const char *prefix, int algo_type, int block_size)
{
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

	char *str, *str2, *str3;
	int64_t l_pac;

	str = (char*) calloc(strlen(prefix) + 10, 1);
	str2 = (char*) calloc(strlen(prefix) + 10, 1);
	str3 = (char*) calloc(strlen(prefix) + 10, 1);

	// nucleotide indexing
	gzFile fp = xzopen(fa, "r");

	l_pac = bns_fasta2bntseq(fp, prefix, 0);
	err_gzclose(fp);

	if (algo_type == 0)
		algo_type = l_pac > 50000000 ? 2 : 3; // set the algorithm for generating BWT
	{
		strcpy(str, prefix);
		strcat(str, ".pac");
		strcpy(str2, prefix);
		strcat(str2, ".bwt");

		if (algo_type == 2)
			bwt_bwtgen2(str, str2, block_size);
		else if (algo_type == 1 || algo_type == 3)
		{
			bwt_t *bwt;
			bwt = bwt_pac2bwt(str, algo_type == 3);
			bwt_dump_bwt(str2, bwt);
			bwt_destroy(bwt);
		}

	}
	{
		bwt_t *bwt;
		strcpy(str, prefix);
		strcat(str, ".bwt");
		bwt = bwt_restore_bwt(str);
		bwt_bwtupdate_core(bwt);
		bwt_dump_bwt(str, bwt);
		bwt_destroy(bwt);
	}
	{
		gzFile fp = xzopen(fa, "r");
		l_pac = bns_fasta2bntseq(fp, prefix, 1);
		err_gzclose(fp);
	}
	{
		bwt_t *bwt;
		strcpy(str, prefix);
		strcat(str, ".bwt");
		strcpy(str3, prefix);
		strcat(str3, ".sa");
		bwt = bwt_restore_bwt(str);
		bwt_cal_sa(bwt, 32);
		bwt_dump_sa(str3, bwt);
		bwt_destroy(bwt);
	}
	free(str3);
	free(str2);
	free(str);
	return 0;
}

static void bwt_dump_bwt(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	err_fwrite(&bwt->primary, sizeof(uint64_t), 1, fp);
	err_fwrite(bwt->L2 + 1, sizeof(uint64_t), 4, fp);
	err_fwrite(bwt->bwt, 4, bwt->bwt_size, fp);
	err_fflush(fp);
	err_fclose(fp);
}
static void bwt_dump_sa(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	err_fwrite(&bwt->primary, sizeof(uint64_t), 1, fp);
	err_fwrite(bwt->L2 + 1, sizeof(uint64_t), 4, fp);
	err_fwrite(&bwt->sa_intv, sizeof(uint64_t), 1, fp);
	err_fwrite(&bwt->seq_len, sizeof(uint64_t), 1, fp);
	err_fwrite(bwt->sa + 1, sizeof(uint64_t), bwt->n_sa - 1, fp);
	err_fflush(fp);
	err_fclose(fp);
}

static int64_t bwa_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	ubyte_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int) c;
}

static bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is)
{
	bwt_t *bwt;
	ubyte_t *buf, *buf2;
	int64_t i, pac_size;
	FILE *fp;

	// initialization
	bwt = (bwt_t*) calloc(1, sizeof(bwt_t));
	bwt->seq_len = bwa_seq_len(fn_pac);
	bwt->bwt_size = (bwt->seq_len + 15) >> 4;
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (bwt->seq_len >> 2) + ((bwt->seq_len & 3) == 0 ? 0 : 1);
	buf2 = (ubyte_t*) calloc(pac_size, 1);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	memset(bwt->L2, 0, 5 * 4);
	buf = (ubyte_t*) calloc(bwt->seq_len + 1, 1);
	for (i = 0; i < bwt->seq_len; ++i)
	{
		buf[i] = buf2[i >> 2] >> ((3 - (i & 3)) << 1) & 3;
		++bwt->L2[1 + buf[i]];
	}
	for (i = 2; i <= 4; ++i)
		bwt->L2[i] += bwt->L2[i - 1];
	free(buf2);

	// Burrows-Wheeler Transform
	if (use_is)
	{
		bwt->primary = is_bwt(buf, bwt->seq_len);
	}
	else
	{
		rope_t *r;
		int64_t x;
		rpitr_t itr;
		const uint8_t *blk;

		r = rope_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN);
		for (i = bwt->seq_len - 1, x = 0; i >= 0; --i)
		{
			int c = buf[i] + 1;
			x = rope_insert_run(r, x, c, 1, 0) + 1;
			while (--c >= 0)
				x += r->c[c];
		}
		bwt->primary = x;
		rope_itr_first(r, &itr);
		x = 0;
		while ((blk = rope_itr_next_block(&itr)) != 0)
		{
			const uint8_t *q = blk + 2, *end = blk + 2 + *rle_nptr(blk);
			while (q < end)
			{
				int c = 0;
				int64_t l;
				rle_dec1(q, c, l);
				for (i = 0; i < l; ++i)
					buf[x++] = c - 1;
			}
		}
		rope_destroy(r);
	}
	bwt->bwt = (uint32_t*) calloc(bwt->bwt_size, 4);
	for (i = 0; i < bwt->seq_len; ++i)
		bwt->bwt[i >> 4] |= buf[i] << ((15 - (i & 15)) << 1);
	free(buf);
	return bwt;
}

static int bwa_pac2bwt(int argc, char *argv[]) // the "pac2bwt" command; IMPORTANT: bwt generated at this step CANNOT be used with BWA. bwtupdate is required!
{
	bwt_t *bwt;
	int c, use_is = 1;
	while ((c = getopt(argc, argv, "d")) >= 0)
	{
		switch (c)
		{
			case 'd':
				use_is = 0;
				break;
			default:
				return 1;
		}
	}
	if (optind + 2 > argc)
	{
		fprintf(stderr, "Usage: bwa pac2bwt [-d] <in.pac> <out.bwt>\n");
		return 1;
	}
	bwt = bwt_pac2bwt(argv[optind], use_is);
	bwt_dump_bwt(argv[optind + 1], bwt);
	bwt_destroy(bwt);
	return 0;
}

static void bwt_bwtupdate_core(bwt_t *bwt)
{
	uint64_t i, k, c[4], n_occ;
	uint32_t *buf;

	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	bwt->bwt_size += n_occ * sizeof(uint64_t); // the new size
	buf = (uint32_t*) calloc(bwt->bwt_size, 4); // will be the new bwt
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i)
	{
		if (i % OCC_INTERVAL == 0)
		{
			memcpy(buf + k, c, sizeof(uint64_t) * 4);
			k += sizeof(uint64_t); // in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4)
		}
		if (i % 16 == 0)
			buf[k++] = bwt->bwt[i / 16]; // 16 == sizeof(uint32_t)/2
		++c[bwt_B00(bwt, i)];
	}
	// the last element
	memcpy(buf + k, c, sizeof(uint64_t) * 4);
	xassert(k + sizeof(uint64_t) == bwt->bwt_size, "inconsistent bwt_size");
	// update bwt
	free(bwt->bwt);
	bwt->bwt = buf;
}

static int bwa_bwtupdate(int argc, char *argv[]) // the "bwtupdate" command
{
	bwt_t *bwt;
	if (argc != 2)
	{
		fprintf(stderr, "Usage: bwa bwtupdate <the.bwt>\n");
		return 1;
	}
	bwt = bwt_restore_bwt(argv[1]);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(argv[1], bwt);
	bwt_destroy(bwt);
	return 0;
}

static int bwa_bwt2sa(int argc, char *argv[]) // the "bwt2sa" command
{
	bwt_t *bwt;
	int c, sa_intv = 32;
	while ((c = getopt(argc, argv, "i:")) >= 0)
	{
		switch (c)
		{
			case 'i':
				sa_intv = atoi(optarg);
				break;
			default:
				return 1;
		}
	}
	if (optind + 2 > argc)
	{
		fprintf(stderr, "Usage: bwa bwt2sa [-i %d] <in.bwt> <out.sa>\n", sa_intv);
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	bwt_cal_sa(bwt, sa_intv);
	bwt_dump_sa(argv[optind + 1], bwt);
	bwt_destroy(bwt);
	return 0;
}
