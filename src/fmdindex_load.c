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

#include "fmdindex_load.h"

#include <string.h>
#include <errno.h>

#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_MAP_INIT_STR(str, int)

static void bwt_gen_cnt_table(bwt_t *bwt);
static void bwt_restore_sa(const char *fn, bwt_t *bwt);

char *bwa_idx_infer_prefix(const char *hint)
{
	char *prefix;
	int l_hint;
	FILE *fp;
	l_hint = strlen(hint);

	prefix = malloc(l_hint + 3 + 4 + 1);
	strcpy(prefix, hint);
	strcpy(prefix + l_hint, ".64.bwt");

	if ((fp = fopen(prefix, "rb")) != 0)
	{
		fclose(fp);
		prefix[l_hint + 3] = 0;
		return prefix;
	}
	else
	{
		strcpy(prefix + l_hint, ".bwt");
		if ((fp = fopen(prefix, "rb")) == 0)
		{
			free(prefix);
			return 0;
		}
		else
		{
			fclose(fp);
			prefix[l_hint] = 0;
			return prefix;
		}
	}
}

uint64_t fread_fix(FILE *fp, uint64_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	const int bufsize = 0x1000000; // 16M block
	uint64_t offset = 0;
	while (size)
	{
		int x = bufsize < size ? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0)
			break;
		size -= x;
		offset += x;
	}
	return offset;
}

bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*) calloc(1, sizeof(bwt_t));
	fp = xopen(fn, "rb");
	err_fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (err_ftell(fp) - sizeof(uint64_t) * 5) >> 2;
	bwt->bwt = (uint32_t*) calloc(bwt->bwt_size, 4);
	err_fseek(fp, 0, SEEK_SET);
	err_fread_noeof(&bwt->primary, sizeof(uint64_t), 1, fp);
	err_fread_noeof(bwt->L2 + 1, sizeof(uint64_t), 4, fp);
	fread_fix(fp, bwt->bwt_size << 2, bwt->bwt);
	bwt->seq_len = bwt->L2[4];
	err_fclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}

bwt_t *bwa_idx_load_bwt(const char *hint)
{
	char *tmp, *prefix;
	bwt_t *bwt;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0)
	{
		fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	tmp = calloc(strlen(prefix) + 5, 1);
	strcat(strcpy(tmp, prefix), ".bwt"); // FM-index
	bwt = bwt_restore_bwt(tmp);
	strcat(strcpy(tmp, prefix), ".sa");  // partial suffix array (SA)
	bwt_restore_sa(tmp, bwt);
	free(tmp);
	free(prefix);
	return bwt;
}
bntseq_t *bns_restore(const char *prefix)
{
	char ann_filename[1024], amb_filename[1024], pac_filename[1024], alt_filename[1024];
	FILE *fp;
	bntseq_t *bns;
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");
	bns = bns_restore_core(ann_filename, amb_filename, pac_filename);
	if (bns == 0)
		return 0;
	if ((fp = fopen(strcat(strcpy(alt_filename, prefix), ".alt"), "r")) != 0)
	{
		// read .alt file if present
		char str[1024];
		khash_t(str) *h;
		int c, i, absent;
		khint_t k;
		h = kh_init(str);
		for (i = 0; i < bns->n_seqs; ++i)
		{
			k = kh_put(str, h, bns->anns[i].name, &absent);
			kh_val(h, k) = i;
		}
		i = 0;
		while ((c = fgetc(fp)) != EOF)
		{
			if (c == '\t' || c == '\n' || c == '\r')
			{
				str[i] = 0;
				if (str[0] != '@')
				{
					k = kh_get(str, h, str);
					if (k != kh_end(h))
						bns->anns[kh_val(h, k)].is_alt = 1;
				}
				while (c != '\n' && c != EOF)
					c = fgetc(fp);
				i = 0;
			}
			else
				str[i++] = c; // FIXME: potential segfault here
		}
		kh_destroy(str, h);
		fclose(fp);
	}
	return bns;
}

bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[8192];
	FILE *fp;
	const char *fname;
	bntseq_t *bns;
	long long xx;
	int i;
	int scanres;
	bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = xopen(fname = ann_filename, "r");
		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		if (scanres != 3)
			goto badread;
		bns->l_pac = xx;
		bns->anns = (bntann1_t*) calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i)
		{
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			scanres = fscanf(fp, "%u%s", &p->gi, str);
			if (scanres != 2)
				goto badread;
			p->name = strdup(str);
			// read fasta comments 
			while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF)
				*q++ = c;
			while (c != '\n' && c != EOF)
				c = fgetc(fp);
			if (c == EOF)
			{
				scanres = EOF;
				goto badread;
			}
			*q = 0;
			if (q - str > 1 && strcmp(str, " (null)") != 0)
				p->anno = strdup(str + 1); // skip leading space
			else
				p->anno = strdup("");
			// read the rest
			scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			if (scanres != 3)
				goto badread;
			p->offset = xx;
		}
		err_fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		fp = xopen(fname = amb_filename, "r");
		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		if (scanres != 3)
			goto badread;
		l_pac = xx;
		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = bns->n_holes ? (bntamb1_t*) calloc(bns->n_holes, sizeof(bntamb1_t)) : 0;
		for (i = 0; i < bns->n_holes; ++i)
		{
			bntamb1_t *p = bns->ambs + i;
			scanres = fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			if (scanres != 3)
				goto badread;
			p->offset = xx;
			p->amb = str[0];
		}
		err_fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = xopen(pac_filename, "rb");
	}
	return bns;

	badread: if (EOF == scanres)
	{
		err_fatal(__func__, "Error reading %s : %s\n", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
	}
	err_fatal(__func__, "Parse error reading %s\n", fname);
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0)
		return;
	free(bwt->sa);
	free(bwt->bwt);
	free(bwt);
}

static void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i)
	{
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i & 3) == j) + ((i >> 2 & 3) == j) + ((i >> 4 & 3) == j) + (i >> 6 == j))
			<< (j << 3);
		bwt->cnt_table[i] = x;
	}
}

static void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	uint64_t primary;

	fp = xopen(fn, "rb");
	err_fread_noeof(&primary, sizeof(uint64_t), 1, fp);
	xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	err_fread_noeof(skipped, sizeof(uint64_t), 4, fp); // skip
	err_fread_noeof(&bwt->sa_intv, sizeof(uint64_t), 1, fp);
	err_fread_noeof(&primary, sizeof(uint64_t), 1, fp);
	xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (uint64_t*) calloc(bwt->n_sa, sizeof(uint64_t));
	bwt->sa[0] = -1;

	fread_fix(fp, sizeof(uint64_t) * (bwt->n_sa - 1), bwt->sa + 1);
	err_fclose(fp);
}
