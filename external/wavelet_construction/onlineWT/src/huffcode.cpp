#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alphabet.hpp"
#include "arrayutil.hpp"
#include "binheap.hpp"
#include "bitarray.hpp"
#include "bytearray.hpp"
#include "cocadautil.hpp"
#include "cstringutil.hpp"
#include "dynstr.hpp"
#include "huffcode.hpp"
#include "mathutil.hpp"
#include "strstream.hpp"

struct _hufftnode {
    size_t     chr_rank;
    hufftnode *left;
    hufftnode *right;
    byte_t       *ab_mask;
};


struct _huffcode {
    alphabet  *ab;
    size_t        size;
    hufftnode *tree;
    bincode   *code;
};


typedef struct __nodefreq {
    size_t node;
    size_t freq;
} _nodefreq;


int _nodefreq_comp_fn(void *p1, void *p2)
{
    _nodefreq *nf1, *nf2;
    nf1 = *(_nodefreq **)p1;
    nf2 = *(_nodefreq **)p2;
    if (nf1->freq < nf2->freq)
        return -1;
    else if (nf1->freq == nf2->freq)
        return 0;
    else
        return +1;
}


void _fill_code_table(huffcode *hcode, hufftnode *node, size_t code_len,
                      char *code)
{
    if (hufftnode_is_leaf(node)) {
        hcode->code[node->chr_rank].rawcode = bitarr_new_from_str(code, code_len);
        hcode->code[node->chr_rank].code_len = code_len;
        //printf("code of %c = %s\n",ab_char(hcode->ab, node->chr_rank), code);
    }
    else {
        code[code_len] = '0';
        _fill_code_table(hcode, hufftnode_left(node), code_len+1, code);
        code[code_len] = '1';
        _fill_code_table(hcode, hufftnode_right(node), code_len+1, code);
        code[code_len] = '\0';
    }
}


huffcode *huffcode_new(alphabet *ab, size_t freqs[])
{
    huffcode *hcode;
    size_t next, ab_bytesize;
    binheap *nfheap;

    hcode = NEW(huffcode);
    hcode->ab = ab;
    hcode->size = ab_size(ab);

    ab_bytesize = (size_t)multceil(hcode->size, BYTESIZE);

    hcode->tree = NEW_ARRAY(hufftnode, 2*hcode->size-1);
    for (size_t i = 0; i < hcode->size; i++) {
        hcode->tree[i].chr_rank = i;
        hcode->tree[i].left     = &hcode->tree[i];
        hcode->tree[i].right    = &hcode->tree[i];
        hcode->tree[i].ab_mask  = bytearr_new(ab_bytesize);
        bitarr_set_bit(hcode->tree[i].ab_mask, i, 1);
    }

    nfheap = binheap_new(&_nodefreq_comp_fn, sizeof(_nodefreq *), MIN_HEAP);

    for (size_t i=0; i<hcode->size; i++) {
        _nodefreq *new_nf = NEW(_nodefreq);
        new_nf->node = i;
        new_nf->freq = freqs[i];
        binheap_push(nfheap, new_nf);
    }
    next = hcode->size;
    while (binheap_size(nfheap)>1) {
        _nodefreq *smallest, *snd_smallest, *new_nf;
        smallest = (_nodefreq*)binheap_pop(nfheap);
        snd_smallest = (_nodefreq*)binheap_pop(nfheap);
        new_nf = NEW(_nodefreq);
        new_nf->node = next;
        new_nf->freq = smallest->freq + snd_smallest->freq;
        binheap_push(nfheap, new_nf);
        hcode->tree[next].left = &hcode->tree[smallest->node];
        hcode->tree[next].right = &hcode->tree[snd_smallest->node];
        hcode->tree[next].ab_mask = bytearr_new(ab_bytesize);
        hcode->tree[next].chr_rank = hcode->size;
        bitarr_or( hcode->tree[next].ab_mask, 
                   hcode->tree[smallest->node].ab_mask, hcode->size );
        bitarr_or( hcode->tree[next].ab_mask, 
                   hcode->tree[snd_smallest->node].ab_mask, hcode->size );
        FREE(smallest);
        FREE(snd_smallest);
        next++;
    }
    //assert(next==(2*hcode->size-1));

    hcode->code = NEW_ARRAY(bincode, hcode->size);
    char *chrcode = cstr_new(hcode->size);
    _fill_code_table(hcode, huffcode_tree(hcode), 0, chrcode);
    FREE(chrcode);
    return hcode;
}


huffcode *huffcode_new_from_str(alphabet *ab, char *src)
{
    huffcode *hcode;
    size_t *freqs = ab_count(ab, src, strlen(src));
    hcode = huffcode_new(ab, freqs);
    FREE(freqs);
    return hcode;
}


huffcode *huffcode_new_from_stream(alphabet *ab, strstream *sst)
{
    huffcode *hcode;
    size_t *freqs = ab_count_stream(ab, sst);
    hcode = huffcode_new(ab, freqs);
    FREE(freqs);
    return hcode;
}


huffcode *huffcode_new_online_from_stream(strstream *sst)
{
    strstream_reset(sst);
    size_t all_len = (UCHAR_MAX+1);
    size_t *all_freqs = NEW_ARRAY(size_t, all_len);
    size_t ab_len = 0;
    FILL_ARRAY(all_freqs, 0,all_len ,0);
    for (int c; (c=strstream_getc(sst))!=EOF;) {
        ab_len += (all_freqs[c]==0)?1:0;
        all_freqs[c] += 1;
    }
    char *abstr = cstr_new(ab_len);
    size_t *ab_freqs = NEW_ARRAY(size_t, ab_len);
    for (size_t i=0, k=0; i<all_len; i++) {
        if (all_freqs[i]>0) {
            abstr[k]  = (char)i;
            ab_freqs[k] = all_freqs[i];
            k++;
        }
    }
    alphabet *ab = alphabet_new(ab_len, abstr);
    huffcode *hcode= huffcode_new(ab, ab_freqs);
    FREE(ab_freqs);
    return hcode;
}



void huffcode_free(huffcode *hcode)
{
    for (size_t i=0; i<hcode->size; i++) {
        if (hcode->code[i].rawcode)
            FREE(hcode->code[i].rawcode);
    }
    FREE(hcode->code);
    for (size_t i=0; i<(2*hcode->size)-1; i++) {
        if (hcode->tree[i].ab_mask)
            FREE(hcode->tree[i].ab_mask);
    }
    FREE(hcode->tree);
    FREE(hcode);
}


void _print_htree(huffcode *hc, hufftnode *node, size_t level, char *code)
{
    char *space = cstr_new(4*level);
    cstr_fill(space, 0, 4*level, ' ');
    if (hufftnode_is_leaf(node)) {
        printf("%s[%p code=%s chr=%c]\n", space, node, code, ab_char(hc->ab,
                node->chr_rank));
        //bytearr_print(hufftnode_ab_mask(node), (size_t)mult_ceil(ab_size(hc->ab), BYTESIZE), 4, space);
    }
    else {
        printf("%s[%p code=%s]\n", space, node, code);
        //bytearr_print(hufftnode_ab_mask(node), (size_t)mult_ceil(ab_size(hc->ab), BYTESIZE), 4, space);
        char *ccode = cstr_new(level+1);
        strcpy(ccode, code);
        ccode[level] = '0';
        _print_htree(hc, hufftnode_left(node), level+1, ccode);
        ccode[level] = '1';
        _print_htree(hc, hufftnode_right(node), level+1, ccode);
    }
    FREE(space);
}


void huffcode_print(huffcode *hcode)
{
    printf("huffcode@%p\n",(void *)hcode);
    _print_htree(hcode, huffcode_tree(hcode), 0, "");
    printf ("codes:\n");
    for (size_t i=0; i<hcode->size; i++) {
        printf("%c: len=%zu raw=",ab_char(hcode->ab, i), hcode->code[i].code_len);
        bytearr_print( hcode->code[i].rawcode,
                       (size_t)multceil(hcode->code[i].code_len, BYTESIZE),
                       (size_t)multceil(hcode->code[i].code_len, BYTESIZE),
                       "");
    }
}


bincode huffcode_encode(huffcode *hcode, char *str, size_t len)
{
    size_t ba_size = MAX(4*BYTESIZE,
                         ((size_t)(len*3/BYTESIZE))*BYTESIZE);//GOT TO BE A MULTIPLE OF BYTESIZE
    byte_t *code = NEW_ARRAY(byte_t, ba_size/BYTESIZE);
    size_t code_size = 0;
    bincode charcode;
    for (size_t i=0; i<len; i++) {
        charcode = hcode->code[ab_rank(hcode->ab, str[i])];
        if (code_size+charcode.code_len >= ba_size) {
            ba_size *= 2;
            code = (byte_t*)realloc(code, ba_size/BYTESIZE);
        }
        bitarr_write(code, code_size, charcode.rawcode, 0, charcode.code_len);
        code_size += charcode.code_len;
    }
    code = (byte_t*)realloc(code, (size_t)multceil(code_size, BYTESIZE));
    bincode ret = {.rawcode=code, .code_len=code_size};
    return ret;
}


char *huffcode_decode(huffcode *hcode, byte_t *code, size_t code_len)
{
    size_t smaxlen = MAX(16, (size_t)(code_len/(size_t)log2(hcode->size)));
    char *ret = cstr_new(smaxlen);
    size_t slen=0;
    hufftnode *cur = huffcode_tree(hcode);
    size_t i=0;
    while (i<code_len) {
        if (bitarr_get_bit(code, i)) {
            cur = cur->right;
        }
        else {
            cur = cur->left;
        }
        if (hufftnode_is_leaf(cur)) {
            if (slen==smaxlen) {
                smaxlen += MAX(16, (size_t)((code_len-i)/(size_t)log2(hcode->size)));
                ret = (char*)realloc(ret, smaxlen*sizeof(char)+1);
                memset(ret+(slen*sizeof(char)), '\0', (smaxlen-slen)*sizeof(char)+1);
                //printf("smaxlen=%zu\n",smaxlen);
            }
            ret[slen++] = ab_char(hcode->ab, cur->chr_rank);
            cur = huffcode_tree(hcode);
        }
        i++;
    }
    //printf("smaxlen=%zu slen=%zu\n",smaxlen, slen);
    ret = (char*)realloc(ret, slen*sizeof(char)+1);
    return ret;
}


hufftnode *huffcode_tree(huffcode *code)
{
    return code->tree+(2*code->size)-2;
}

alphabet *huffcode_ab(huffcode *code)
{
    return code->ab;
}


bool hufftnode_is_leaf(hufftnode *node)
{
    return node->left==node->right;
}


hufftnode *hufftnode_left(hufftnode *node)
{
    return node->left;
}


hufftnode *hufftnode_right(hufftnode *node)
{
    return node->right;
}


byte_t *hufftnode_ab_mask(hufftnode *node)
{
    return node->ab_mask;
}

size_t hufftnode_char_rank(hufftnode *node)
{
    return node->chr_rank;
}
