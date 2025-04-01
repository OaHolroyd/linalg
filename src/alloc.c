#include "alloc.h"

#include <stdlib.h>

double **malloc_d2d(const int ni, const int nj) {
  /* allocate row memory */
  double **p_2arr = malloc(ni * sizeof(double *));
  if (!p_2arr) {
    return NULL;
  }

  /* allocate main memory */
  double *mem = malloc(ni * nj * sizeof(double));
  if (!mem) {
    free(p_2arr);
    return NULL;
  }

  /* match rows to memory */
  for (int i = 0; i < ni; i++) {
    p_2arr[i] = &(mem[i * nj]);
  } // i end

  return p_2arr;
}

void internal_free_2d(void **arr) {
  if (arr) {
    free(arr[0]); // free the main memory
    free(arr); // free the row pointers
  }
}
