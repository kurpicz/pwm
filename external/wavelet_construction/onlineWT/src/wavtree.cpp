#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "alphabet.hpp"
#include "arrayutil.hpp"
#include "bitsandbytes.hpp"
#include "bitvector.hpp"
#include "bytearray.hpp"
#include "cocadautil.hpp"
#include "csrsbitarray.hpp"
#include "cstringutil.hpp"
#include "dynstr.hpp"
#include "hashmap.hpp"
#include "huffcode.hpp"
#include "mathutil.hpp"
#include "wavtree.hpp"


/****************************************************************************
 * Convenience constants & hash functions                                   *
 ****************************************************************************/

typedef uintmax_t charcode;
static const charcode NULL_CODE = UINTMAX_MAX;
static const size_t  CHARCODE_BIT = sizeof(intmax_t)*BYTESIZE;
static const byte_t LEFT  = 0x0;
static const byte_t RIGHT = 0x1;


static size_t _charcode_hash(void *c) {
    return (size_t)(*((charcode *)c)); 
}


static bool _charcode_equals(void *c1, void *c2) {
    return ((*((charcode *)c1))==(*((charcode *)c2))); 
}


static inline void _hashmap_set_charcode_to_char(hashmap *hm, charcode u, char c)
{
    charcode *pu = NEW(charcode);
    *pu = u;
    char *pc = NEW(char);
    *pc = c;
    hashmap_set(hm, pu, pc);
}


/****************************************************************************
 * Template wavelet tree apparatus used for construction                    *
 ****************************************************************************/

typedef struct _tmp_wtnode {
    byte_t    nxt_chd;  
    bitvector *bv;
    struct _tmp_wtnode *chd[2];
}
tmp_wtnode;


typedef struct _tmp_wavtree {
    wtshape     shape;
    huffcode   *hcode;
    tmp_wtnode *tmp_root;
    size_t      nchars;
    size_t     *chcode_table;
} 
tmp_wavtree;


static tmp_wtnode *tmp_wtnode_new(size_t size)
{
    tmp_wtnode *node;
    node = NEW(tmp_wtnode);
    node->bv = bitvector_new(size);
    node->nxt_chd = LEFT;
    node->chd[LEFT] = NULL;
    node->chd[RIGHT] = NULL;
    return node;
}


static void tmp_wtnode_free(tmp_wtnode *node)
{
    if (node == NULL) return;
    tmp_wtnode_free(node->chd[LEFT]);
    tmp_wtnode_free(node->chd[RIGHT]);
    FREE(node);
}


static tmp_wavtree *tmp_wavtree_new(wtshape shape)
{
    tmp_wavtree *tmp_wt = NEW(tmp_wavtree);
    tmp_wt->shape = shape;
    tmp_wt->tmp_root = tmp_wtnode_new(0);
    tmp_wt->nchars = 0;
    tmp_wt->chcode_table = NEW_ARRAY(size_t, (UCHAR_MAX+1));
    FILL_ARRAY(tmp_wt->chcode_table, 0, (UCHAR_MAX+1), NULL_CODE);
    tmp_wt->hcode = NULL;
    return tmp_wt;
}


static tmp_wtnode *_tmp_wavtree_init_bal( alphabet *ab, size_t l, size_t r,  
                                          charcode *chcode_table, short depth )
{
    // assert about depth's data type
    assert(depth < CHARCODE_BIT); 

    // no need for an explicit node for an unary alphabet
    if (r-l<=1) {
    	chcode_table[(size_t)ab_char(ab, l)] = l;
        //printf("--> added bal chcode[%c]=%lu\n", ab_char(ab, l), (unsigned long)l );
        return NULL;
    }

    tmp_wtnode *node = tmp_wtnode_new(0);

    // partition alphabet
    size_t mid = (l+r)/2;

    node->chd[LEFT] = _tmp_wavtree_init_bal(ab, l, mid, chcode_table, depth+1);
    node->chd[RIGHT] = _tmp_wavtree_init_bal(ab, mid, r, chcode_table, depth+1);
    return node;
}


static void tmp_wavtree_init_bal(tmp_wavtree *tmp_wt, alphabet *ab)
{
    tmp_wt->nchars = ab_size(ab);
    tmp_wtnode_free(tmp_wt->tmp_root);
    tmp_wt->tmp_root = _tmp_wavtree_init_bal( ab, 0, tmp_wt->nchars, 
                                              tmp_wt->chcode_table, 0);
}


static tmp_wtnode *_tmp_wavtree_init_huff( hufftnode *htnode, huffcode *hcode,
                                           charcode *chcode_table, short depth, 
                                           charcode chcode )
{ 
    // assert about depth's data type
    assert(depth < CHARCODE_BIT); 

    // no need for an explicit node for an unary alphabet
    if (hufftnode_is_leaf(htnode)) {
        chcode_table[ (size_t)ab_char( huffcode_ab(hcode), 
                                    hufftnode_char_rank(htnode)) ] = chcode;
        //printf("--> added huff chcode[%c]=%lu\n", ab_char( huffcode_ab(hcode), hufftnode_char_rank(htnode)), (unsigned long)chcode );
        return NULL;
    }

    tmp_wtnode *node = tmp_wtnode_new(0);
    node->chd[LEFT]  = _tmp_wavtree_init_huff( hufftnode_left(htnode), hcode, 
                                                 chcode_table, depth+1, chcode);
    node->chd[RIGHT] = _tmp_wavtree_init_huff( hufftnode_right(htnode), hcode, 
                                                 chcode_table, depth+1, 
                                                 chcode | ((charcode)1 << depth));
    return node;
}


static void tmp_wavtree_init_huff(tmp_wavtree *tmp_wt, huffcode *hcode)
{
    tmp_wt->nchars = ab_size(huffcode_ab(hcode));
    tmp_wtnode_free(tmp_wt->tmp_root);
    tmp_wt->tmp_root = _tmp_wavtree_init_huff( huffcode_tree(hcode), hcode, 
                                               tmp_wt->chcode_table, 0, 0 );
    //huffcode_print(hcode);
}


static void tmp_wavtree_free(tmp_wavtree *tmp_wt)
{
    tmp_wtnode_free(tmp_wt->tmp_root);
    if (tmp_wt->hcode!=NULL) 
        huffcode_free(tmp_wt->hcode);
    FREE(tmp_wt->chcode_table);
    FREE(tmp_wt);
}


static void tmp_wavtree_app_char(tmp_wtnode *root, charcode chcode)
{
    byte_t bit;
    while (root != NULL) {
        bit = chcode%2;
        bitvec_append(root->bv, bit);
        root = root->chd[bit];
        chcode >>= 1;
    }
}


static void tmp_wavtree_app_new_char(tmp_wtnode *root, charcode chcode)
{
    tmp_wtnode *node=root, *parent=NULL;
    while (true) {
        if (node!=NULL) {
            bitvec_append(node->bv, node->nxt_chd);
            parent = node;
            node = node->chd[node->nxt_chd];
            parent->nxt_chd ^= 0x1;
        }
        else {
            // case first & second symbols: 'ignore'
            if (chcode<2) {
                return;
            }
            byte_t bit = parent->nxt_chd^0x1;
            size_t q = bitvec_count(parent->bv, bit);
            tmp_wtnode *new_node = tmp_wtnode_new(q);
            bitvec_append_n(new_node->bv, q-1, 0x0);
            bitvec_append(new_node->bv, 0x1);
            parent->chd[bit] = new_node;
            return;
        }
    }
}


static void tmp_wavtree_fill( tmp_wavtree *tmp_wt, strstream *sst, 
                               bool online ) 
{
    size_t old_cap = 0;
    if (online) {
        charcode chcode;
        strstream_reset(sst);
        for (int c; (c=strstream_getc(sst))!=EOF;) {
            if (bitvector_cap(tmp_wt->tmp_root->bv) != old_cap) {
                old_cap = bitvector_cap(tmp_wt->tmp_root->bv);
                //bitvec_print_cap(tmp_wt->tmp_root->bv);
            }
            if ((chcode=tmp_wt->chcode_table[c]) != NULL_CODE) {
                tmp_wavtree_app_char(tmp_wt->tmp_root, chcode);
            } else {
                chcode = (tmp_wt->nchars++);
                tmp_wt->chcode_table[c] = chcode;
                tmp_wavtree_app_new_char(tmp_wt->tmp_root, chcode);
            }
        }
    } 
    else {
        strstream_reset(sst);
        for (int c; (c=strstream_getc(sst))!=EOF;) {
            tmp_wavtree_app_char(tmp_wt->tmp_root, tmp_wt->chcode_table[c]);
        }
    }
}

/*

static void tmp_wavtree_app_n_char(tmp_wtnode *root, size_t n, charcode chcode)
{
    while (root != NULL) {
        byte_t bit = chcode%2;
        bitvec_append_n(root->bv, n, bit);
        root = root->chd[bit];
        chcode >>= 1;
    }
}


static void tmp_wavtree_app_new_n_char(tmp_wtnode *root, size_t n, charcode chcode)
{
    tmp_wtnode *node=root, *parent=NULL;
    while (true) {
        if (node!=NULL) {
            bitvec_append_n(node->bv, n, node->nxt_chd);
            parent = node;
            node = node->chd[node->nxt_chd];
                parent->nxt_chd ^= 0x1;
        }
        else {
            // case first & second symbols: 'ignore'
            if (chcode<2) {
                return;
            }
            byte_t bit = parent->nxt_chd ^0x1;
            size_t q = bitvec_count(parent->bv, bit);
            tmp_wtnode *new_node = tmp_wtnode_new(q);
            bitvec_append_n(new_node->bv, q-n, 0x0);
            bitvec_append_n(new_node->bv, n, 0x1);
            parent->chd[bit] = new_node;
            //parent->nxt_chd ^= 0x1;
            return;
        }
    }
}




static void tmp_wavtree_fill_mult( tmp_wavtree *tmp_wt, strstream *sst, 
                               bool online ) 
{
    if (online) {
        charcode chcode;
        int lastc=EOF, c;
        size_t k=0;
        strstream_reset(sst);
        while ( (c=strstream_getc(sst))!=EOF ) {
            if (c==lastc) {
                k++; 
                continue;
            }
            if (lastc != EOF) {
                if ((chcode=tmp_wt->chcode_table[lastc]) != NULL_CODE) {
                    tmp_wavtree_app_n_char(tmp_wt->tmp_root, k, chcode);
                } else {
                    chcode = (tmp_wt->nchars++);
                    tmp_wt->chcode_table[(char)lastc] = chcode;
                    tmp_wavtree_app_new_n_char(tmp_wt->tmp_root, k, chcode);
                }
            }
            k = 1;
            lastc = c;
        }
        if (lastc != EOF) {
            if ((chcode=tmp_wt->chcode_table[lastc]) != NULL_CODE) {
                tmp_wavtree_app_n_char(tmp_wt->tmp_root, k, chcode);
            } else {
                chcode = (tmp_wt->nchars++);
                tmp_wt->chcode_table[(char)lastc] = chcode;
                tmp_wavtree_app_new_n_char(tmp_wt->tmp_root, k, chcode);
            }
        }
    }
    else {
        charcode chcode;
        int lastc=EOF, c;
        size_t k=0;
        strstream_reset(sst);
        while ( (c=strstream_getc(sst))!=EOF ) {
            if (c==lastc) {
                k++; 
                continue;
            }
            if (lastc != EOF) {
                chcode=tmp_wt->chcode_table[lastc];
                tmp_wavtree_app_n_char(tmp_wt->tmp_root, k, chcode);
            }
            k = 1;
            lastc = c;
        }
        if (lastc != EOF) {
            chcode=tmp_wt->chcode_table[lastc];
            tmp_wavtree_app_n_char(tmp_wt->tmp_root, k, chcode);
        }
    }
}
*/
 
void tmp_wtnode_print(tmp_wtnode *node, size_t depth)
{
    size_t i;
    char *margin = cstr_new(2*depth+2);
    for (i=0; i<2*depth; i++) {
        margin[i] = ' ';
    }
    margin[i++] = '|';
    margin[i++] = ' ';
    margin[i++] = '\0';
    if (node==NULL) {
        printf ("%s@tmp_wtree_node NULL\n",margin);
    }
    else {
        printf ("%s@tmp_wtree_node %p\n",margin, node);
        printf ("%ssize: %zu\n", margin, bitvec_len(node->bv));
        printf ("%snxt_chd: %c\n",margin, node->nxt_chd?'1':'0');
        printf ("%sbits: \n", margin);
        bitvec_print(node->bv, 8);
        tmp_wtnode_print(node->chd[LEFT], depth+1);
        tmp_wtnode_print(node->chd[RIGHT], depth+1);
    }
}


void tmp_wavtree_print(tmp_wavtree *tmp_wt)
{
    printf ("tmp_wavtree@%p\n",tmp_wt);
    tmp_wtnode_print(tmp_wt->tmp_root, 0);
}


/****************************************************************************
 * Wavelet tree ADT and construction                                        *
 ****************************************************************************/

typedef struct _wtnode {
    csrsbitarray *ba;
    struct _wtnode *left;
    struct _wtnode *right;
}
wtnode;


struct _wavtree {
    wtshape     shape;
    wtnode     *root;
    alphabet   *ab;
    charcode   *chcode_table;
    hashmap    *char_table;
};


static wtnode *_wavtree_build(tmp_wtnode *temp_node)
{
    if (temp_node == NULL)
        return NULL;
    wtnode *node = NEW(wtnode);
    size_t nbits = (bitvec_len(temp_node->bv));
    byte_t *raw_bits = bitvec_detach(temp_node->bv);
    raw_bits = (byte_t *)realloc(raw_bits, (size_t)multceil(nbits, BYTESIZE));
    node->ba = csrsbitarr_new(raw_bits, nbits);
    node->left  = _wavtree_build(temp_node->chd[LEFT] );
    node->right = _wavtree_build(temp_node->chd[RIGHT]);
    return node;
}



static wavtree *wavtree_build(tmp_wavtree *tmp_wt)
{
    wavtree *wt = NEW(wavtree);
    wt->shape = tmp_wt->shape;
    wt->root = _wavtree_build(tmp_wt->tmp_root);
    wt->char_table = hashmap_new(_charcode_hash, _charcode_equals);
    wt->chcode_table = NEW_ARRAY(charcode, (UCHAR_MAX+1));
    FILL_ARRAY(wt->chcode_table, 0, (UCHAR_MAX+1), NULL_CODE);
    char *ab_symb = cstr_new(tmp_wt->nchars);
    
    for (int c=0, n=(UCHAR_MAX+1), k=0; c<n; c++) {
        if (tmp_wt->chcode_table[c] != NULL_CODE) {
            ab_symb[k++] = (char)c;
            wt->chcode_table[c] = tmp_wt->chcode_table[c];
            _hashmap_set_charcode_to_char(wt->char_table, tmp_wt->chcode_table[c], (char)c);
        }
    }

    wt->ab = alphabet_new(tmp_wt->nchars, ab_symb);
    
    return wt;
}


wavtree *wavtree_new (alphabet *ab, char *str, size_t len, wtshape shape)
{
    wavtree *wt = NULL;
    tmp_wavtree *tmp_wt;
    strstream *sst;
    switch (shape) {
    case WT_BALANCED:
        tmp_wt = tmp_wavtree_new(shape);
        sst = strstream_open_str(str, len);
        tmp_wavtree_init_bal(tmp_wt, ab);
        tmp_wavtree_fill(tmp_wt, sst, false);
        strstream_close(sst);
        wt = wavtree_build(tmp_wt);
        tmp_wavtree_free(tmp_wt);
        break;
    case WT_HUFFMAN:
        tmp_wt = tmp_wavtree_new(shape);
        tmp_wt->hcode = huffcode_new_from_str(ab, str);
        tmp_wavtree_init_huff(tmp_wt, tmp_wt->hcode);
        sst = strstream_open_str(str, len);
        tmp_wavtree_fill(tmp_wt, sst, false);//huffman online mode not supported
        strstream_close(sst);
        wt = wavtree_build(tmp_wt);
        tmp_wavtree_free(tmp_wt);
        break;
    }
    return wt;
}


wavtree *wavtree_new_from_stream (alphabet *ab, strstream *sst, wtshape shape)
{
    wavtree *wt = NULL;
    tmp_wavtree *tmp_wt;
    switch (shape) {
    case WT_BALANCED:
        tmp_wt = tmp_wavtree_new(shape);
        tmp_wavtree_init_bal(tmp_wt, ab);
        tmp_wavtree_fill(tmp_wt, sst, false);
        wt = wavtree_build(tmp_wt);
        tmp_wavtree_free(tmp_wt);
        break;
    case WT_HUFFMAN:
        tmp_wt = tmp_wavtree_new(shape);
        tmp_wt->hcode = huffcode_new_from_stream(ab, sst);
        tmp_wavtree_init_huff(tmp_wt, tmp_wt->hcode);
        tmp_wavtree_fill(tmp_wt, sst, false);//huffman online mode not supported
        wt = wavtree_build(tmp_wt);
        tmp_wavtree_free(tmp_wt);
        break;
    }
    return wt;
}


wavtree *wavtree_new_online_from_stream (strstream *sst, wtshape shape)
{
    wavtree *wt = NULL;
    tmp_wavtree *tmp_wt;
    switch (shape) {
    case WT_BALANCED:
        tmp_wt = tmp_wavtree_new(shape);
        tmp_wavtree_fill(tmp_wt, sst, true);
        wt = wavtree_build(tmp_wt);
        tmp_wavtree_free(tmp_wt);
        break;
    case WT_HUFFMAN:
        tmp_wt = tmp_wavtree_new(shape);
        tmp_wt->hcode = huffcode_new_online_from_stream(sst);
        tmp_wavtree_init_huff(tmp_wt, tmp_wt->hcode);
        tmp_wavtree_fill(tmp_wt, sst, false);//huffman online mode not supported
        wt = wavtree_build(tmp_wt);
        tmp_wavtree_free(tmp_wt);
        break;
    }
    return wt;
}


void _wtnode_free(wtnode *node)
{
    if (node==NULL) return;
    csrsbitarr_free(node->ba, true);
    _wtnode_free(node->left);
    _wtnode_free(node->right);
    FREE(node);
}


void wavtree_free(wavtree *wt)
{
    _wtnode_free(wt->root);
    FREE(wt);
}


/****************************************************************************
 * Wavelet tree operations                                                  *
 ****************************************************************************/

size_t wavtree_rank_pos(wavtree *wt, size_t pos)
{
    wtnode *node = wt->root;
    size_t rank = pos+1;
    while (node!=NULL) {
        if ( csrsbitarr_get(node->ba, rank-1) == 0 ) {
            rank = csrsbitarr_rank0(node->ba, rank-1);
            node = node->left;
        }
        else {
            rank = csrsbitarr_rank1(node->ba, rank-1);
            node = node->right;
        }
    }
    return rank;
}


static size_t _wavtree_rank_bal(wavtree *wt, size_t pos, char c)
{
    wtnode *node = wt->root;
    charcode chcode = wt->chcode_table[(size_t)c];
    size_t rank = pos+1;
    byte_t bit;
    while (rank>0 && node!=NULL) {
        bit = (byte_t)(chcode&1);
        rank = csrsbitarr_rank(node->ba, rank-1, bit);
        node = (bit==0) ? node->left: node->right;
        chcode>>=1;
    }
    return rank;
}

static size_t _wavtree_rank_huff(wavtree *wt, size_t pos, char c)
{
    wtnode *node = wt->root;
    charcode chcode = wt->chcode_table[(size_t)c];
    size_t rank = pos+1;
    byte_t bit;
    while (rank>0 && node!=NULL) {
        bit = (byte_t)(chcode&1);
        rank = csrsbitarr_rank(node->ba, rank-1, bit);
        node = (bit==0) ? node->left: node->right;
        chcode>>=1;
    }
    return rank;
}

typedef size_t (*_rank_func)(wavtree *wt, size_t pos, char c);


static _rank_func wavtree_rank_func[2] = {
    _wavtree_rank_bal,
    _wavtree_rank_huff
};


size_t wavtree_rank(wavtree *wt, size_t pos, char c){
    return wavtree_rank_func[wt->shape](wt, pos, c);
}


static size_t _wtree_select(wtnode *node, charcode chcode, size_t rank)
{
    if (node == NULL)
        return rank;

    bool bit = ((chcode & 1) == 1);
    wtnode *son = bit ? node->right : node->left;
    size_t sel = _wtree_select(son, chcode >> 1, rank);
    sel = bit ? csrsbitarr_select1(node->ba, sel) 
                : csrsbitarr_select0(node->ba, sel);
    return sel+1;
}


static size_t _wavtree_select_bal(wavtree *wt, char c, size_t rank)
{
    charcode chcode = wt->chcode_table[(size_t)c];
    return _wtree_select(wt->root, chcode, rank) - 1;
}


static size_t _wavtree_select_huff(wavtree *wt, char c, size_t rank)
{
    charcode chcode = wt->chcode_table[(size_t)c];
    return _wtree_select(wt->root, chcode, rank) - 1;
}


typedef size_t (*_select_func)(wavtree *wt, char c, size_t rank);


static _select_func wavtree_select_func[2] = {
    _wavtree_select_bal,
    _wavtree_select_huff,
};


size_t wavtree_select(wavtree *wt, char c, size_t rank) {
    return wavtree_select_func[wt->shape](wt, c, rank);
}


static charcode _wtree_char(wtnode *node, size_t pos)
{
    charcode chcode = 0;
    short depth = 0;
    while (node != NULL) {
        if (csrsbitarr_get(node->ba, pos) == 0) {
            pos = csrsbitarr_rank0(node->ba, pos) - 1;
            node = node->left;
        }
        else {
            pos = csrsbitarr_rank1(node->ba, pos) - 1;
            chcode |= 1ULL << depth;
            node = node->right;
        }
        depth++;
    }
    return chcode;
}


static char wavtree_char_bal(wavtree *wt, size_t pos)
{
    charcode chcode = _wtree_char(wt->root, pos);
    return ab_char(wt->ab, chcode);
}


static char wavtree_char_huff(wavtree *wt, size_t pos)
{
    charcode chcode = _wtree_char(wt->root, pos);
    return *((char *)hashmap_get(wt->char_table, &chcode));
}


typedef char (*_char_func)(wavtree *wt, size_t pos);


static _char_func wavtree_char_func[2] = {
    wavtree_char_bal,
    wavtree_char_huff,    
};


char wavtree_char(wavtree *wt, size_t pos)
{
    return wavtree_char_func[wt->shape](wt, pos);
}

// Print

void _wt_node_print(wtnode *node, size_t depth)
{
    size_t i;
    char *margin = cstr_new(2*depth+2);
    for (i=0; i<2*depth; i++) {
        margin[i] = ' ';
    }
    margin[i++] = '|';
    margin[i++] = ' ';
    margin[i++] = '\0';
    if (node==NULL) {
        printf ("%s@wtree_node NULL\n",margin);
    }
    else {
        printf ("%s@wtree_node %p\n",margin, node);
        printf ("%sbits: \n", margin);
        bytearr_print(csrsbitarr_data(node->ba),
                      multceil(csrsbitarr_len(node->ba), BYTESIZE),
                      multceil(csrsbitarr_len(node->ba), BYTESIZE),
                      margin);
        _wt_node_print(node->left, depth+1);
        _wt_node_print(node->right, depth+1);
    }
}

void wavtree_print(wavtree *wt)
{
    printf ("wavelet_tree@%p\n",wt);
    _wt_node_print(wt->root, 0);
    printf ("char codes:\n");
    size_t chcode_bits = sizeof(size_t)*BYTESIZE;
    char *chcode = cstr_new(chcode_bits);
    for (int c=0, n=(UCHAR_MAX+1); c<n; c++) {
        if (wt->chcode_table[c] != NULL_CODE) {
            sizet_to_cstr(chcode, wt->chcode_table[c], 'b');
            cstr_revert(chcode, chcode_bits);
            printf ("  %c : %s\n", (char)c, chcode);
        }
    }
}


