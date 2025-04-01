#include "alloc.h"

#include <stdlib.h>

double **malloc_d2d(const int n, const int m) {
  /* allocate row memory */
  double **rows = malloc(n * sizeof(double *));
  if (!rows) {
    return NULL;
  }

  /* allocate main memory */
  double *mem = malloc(n * m * sizeof(double));
  if (!mem) {
    free(rows);
    return NULL;
  }

  /* match rows to memory */
  for (int i = 0; i < n; i++) {
    rows[i] = &(mem[i * m]);
  } // i end

  return rows;
}

double **calloc_d2d(int n, int m) {
  /* allocate row memory */
  double **rows = malloc(n * sizeof(double *));
  if (!rows) {
    return NULL;
  }

  /* allocate main memory */
  double *mem = calloc(n * m, sizeof(double));
  if (!mem) {
    free(rows);
    return NULL;
  }

  /* match rows to memory */
  for (int i = 0; i < n; i++) {
    rows[i] = &(mem[i * m]);
  } // i end

  return rows;
}

void internal_free_2d(void **arr) {
  if (arr) {
    free(arr[0]); // free the main memory
    free(arr); // free the row pointers
  }
}
