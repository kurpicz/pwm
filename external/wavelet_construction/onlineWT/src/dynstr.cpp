#include <string.h>

#include "cocadautil.hpp"
#include "cstringutil.hpp"
#include "dynstr.hpp"


typedef struct _dynstr {
    char *str;
    size_t len;
    size_t capacity;
}
dynstr;

dynstr *dstr_new(size_t init_capacity)
{
    dynstr *ret;
    ret = NEW(dynstr);
    ret->capacity = init_capacity;
    ret->len = 0;
    ret->str = cstr_new(ret->capacity);
    return ret;
}

dynstr *dstr_new_from_str(char *src)
{
    dynstr *ret;
    ret = NEW(dynstr);
    ret->len = ret->capacity = strlen(src);
    ret->str = cstr_new(ret->capacity);
    strcpy(ret->str, src);
    return ret;
}

void dstr_free(dynstr *dstr)
{
    FREE(dstr->str);
    FREE(dstr);
}

size_t dstr_len(dynstr *dstr)
{
    return dstr->len;
}

size_t dstr_capacity(dynstr *dstr)
{
    return dstr->capacity;
}

char dstr_get(dynstr *dstr, size_t pos)
{
    return dstr->str[pos];
}

void dstr_set(dynstr *dstr, size_t pos, char c)
{
    dstr->str[pos] = c;
}

void _double(dynstr *dstr)
{
    dstr->capacity = MAX(1, 2*dstr->capacity);
    dstr->str = (char*)realloc(dstr->str, (dstr->capacity+1)*sizeof(char));
    cstr_fill(dstr->str, dstr->len, dstr->capacity, '\0');
}

void dstr_append(dynstr *dstr, char *suff)
{
    size_t slen = strlen(suff);
    while (dstr->len + slen >= dstr->capacity) {
        _double(dstr);
    }
    strcpy(dstr->str+dstr->len, suff);
    dstr->len += slen;
    dstr->str[dstr->len] = '\0';

}


void dstr_append_char(dynstr *dstr, char c)
{
    if (dstr->len == dstr->capacity) {
        _double(dstr);
    }
    dstr->str[dstr->len] = c;
    dstr->len++;
    dstr->str[dstr->len] = '\0';
}

const char *dstr_as_str(dynstr *dstr)
{
    return dstr->str;
}


char *dstr_detach(dynstr *dstr)
{
    char *str = dstr->str;
    str = (char*)realloc(str, (dstr->len+1));
    FREE(dstr);
    return str;
}

