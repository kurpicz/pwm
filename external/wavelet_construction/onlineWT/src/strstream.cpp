#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cocadautil.hpp"
#include "strstream.hpp"


typedef enum {
    SSTR_STR = 0,
    SSTR_FILE = 1,
} sstream_type;


struct _strstream {
    sstream_type type;
    union {
        FILE *file;
        char *str;
    } src;
    size_t pos;
    size_t slen;
};


strstream *strstream_open_str(char *str, size_t slen)
{
    strstream *sst;
    sst = NEW(strstream);
    sst->type = SSTR_STR;
    sst->src.str = str;
    sst->pos = 0;
    sst->slen = slen;
    return sst;
}


strstream *strstream_open_file(char *filename)
{
    strstream *sst;
    sst = NEW(strstream);
    sst->type = SSTR_FILE;
    sst->src.file = fopen(filename, "r");
    sst->pos = 0;
    return sst;
}

void strstream_reset(strstream *sst)
{
    switch (sst->type) {
    case SSTR_STR:
        sst->pos = 0;
        break;
    case SSTR_FILE:
        rewind(sst->src.file);
        break;
    }
}

bool strstream_end(strstream *sst)
{
    switch (sst->type) {
    case SSTR_STR:
        return sst->slen <= sst->pos;
        break;
    case SSTR_FILE:
        return feof(sst->src.file);
        break;
    default:
        return true;
    }
}

int strstream_getc(strstream *sst)
{
    switch (sst->type) {
    case SSTR_STR:
        if (sst->pos>=sst->slen)
            return EOF;
        else
            return (int)sst->src.str[sst->pos++];
        break;
    case SSTR_FILE:
        return fgetc(sst->src.file);
        break;
    default:
        return '\0';
    }
}

size_t strstream_reads(strstream *sst, char *dest, size_t n)
{
    size_t nread;
    switch (sst->type) {
    case SSTR_STR:
        nread = MIN(n, (sst->pos<sst->slen)?(sst->slen-sst->pos):0);
        strncpy(dest, sst->src.str+sst->pos, nread);
        sst->pos += nread;
        //dest[nread] = '\0';
        return nread;
        break;
    case SSTR_FILE:
        return fread(dest, 1, n, sst->src.file);
        break;
    default:
        return 0;
    }
}

void strstream_close(strstream *sst)
{
    switch (sst->type) {
    case SSTR_STR:
        break;
    case SSTR_FILE:
        fclose(sst->src.file);
        break;
    }
    FREE(sst);
}
