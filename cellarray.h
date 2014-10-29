#define __CELLARRAY_H

#include <stdio.h>
#include <stdlib.h>

//barebones cellaray. Align data for better cpu performance
#define BIN_REFINE_FACTOR     1
#define ALIGNMENT 64

//cellarray does not "need" to be aligned but the individual elements do.
typedef struct{
  DOUBLE *x  __attribute__((aligned(ALIGNMENT)));
  DOUBLE *y  __attribute__((aligned(ALIGNMENT)));
  DOUBLE *z  __attribute__((aligned(ALIGNMENT)));
  int nelements;
}cellarray __attribute__((aligned(ALIGNMENT)));