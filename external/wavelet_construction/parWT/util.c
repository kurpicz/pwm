#include <stdio.h>
#include <stdlib.h>
#include "util.h"

symbol* read_text_from_file(const char* fn, unsigned long* n) {

  FILE* fp = fopen(fn, "r");
  size_t read;

  if (!fp) {
    fprintf(stderr, "Error opening file \"%s\".\n", fn);
    exit(-1);
  }

  fseek(fp, 0L, SEEK_END);
  *n = (unsigned long)ftell(fp);

  symbol* t;
  t = malloc(*n);

  *n = *n / sizeof(symbol); /* Number of symbols*/

  fseek(fp, 0L, SEEK_SET);

  read = fread(t, sizeof(symbol), *n, fp);
  if(read != *n){
    fprintf(stderr, "Error reading file \"%s\".\n", fn);
    exit(-1);
  }

  fclose(fp);

  return t;

}
