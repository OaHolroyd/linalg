#include "io.h"

void mat_fprintf(
    FILE *stream, const char *fmt, const double *A, const int n, const int m
) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m - 1; j++) {
      fprintf(stream, fmt, A[i * m + j]);
      fprintf(stream, " ");
    }
    fprintf(stream, fmt, A[i * m + m - 1]);
    fprintf(stream, "\n");
  }
}

void mat_fprint(FILE *stream, const double *A, const int n, const int m) {
  mat_fprintf(stream, "%4g", A, n, m);
}

void mat_printf(const char *fmt, const double *A, const int n, const int m) {
  mat_fprintf(stdout, fmt, A, n, m);
}

void mat_print(const double *A, const int n, const int m) {
  mat_fprint(stdout, A, n, m);
}

int mat_outputf(
    const char *filename, const char *fmt, const double *A, const int n,
    const int m
) {
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    return 1;
  }

  mat_fprintf(fp, fmt, A, n, m);

  fclose(fp);

  return 0;
}

int mat_output(
    const char *filename, const double *A, const int n, const int m
) {
  return mat_outputf(filename, "%.8f", A, n, m);
}
