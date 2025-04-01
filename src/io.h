#ifndef IO_H
#define IO_H

#include <stdio.h>

/**
 * Print a matrix to a file stream in a specified format.
 *
 * A single space will be added between entries, so this should not be included
 * in the format string. It is suggested that a fixed-width format is used, e.g.
 * "%5.1lf".
 *
 * @param stream output stream
 * @param fmt format string
 * @param A pointer to the flattened matrix data
 * @param n number of rows
 * @param m number of columns
 */
void mat_fprintf(FILE *stream, const char *fmt, const double *A, int n, int m);

/**
 * Print a matrix to a file stream using the default format.
 *
 * @param stream output stream
 * @param A pointer to the flattened matrix data
 * @param n number of rows
 * @param m number of columns
 */
void mat_fprint(FILE *stream, const double *A, int n, int m);

/**
 * Print a matrix to stdout in a specified format.
 *
 * A single space will be added between entries, so this should not be included
 * in the format string. It is suggested that a fixed-width format is used, e.g.
 * "%5.1lf".
 *
 * @param fmt format string
 * @param A pointer to the flattened matrix data
 * @param n number of rows
 * @param m number of columns
 */
void mat_printf(const char *fmt, const double *A, int n, int m);

/**
 * Print a matrix to stdout using the default format.
 *
 * @param A pointer to the flattened matrix data
 * @param n number of rows
 * @param m number of columns
 */
void mat_print(const double *A, int n, int m);

/**
 * Output a matrix to a text file in a specified format.
 *
 * A single space will be added between entries, so this should not be included
 * in the format string. It is suggested that a fixed-width format is used, e.g.
 * "%5.1lf".
 *
 * @param filename name of the file to write to
 * @param fmt format string
 * @param A pointer to the flattened matrix data
 * @param n number of rows
 * @param m number of columns
 * @return 0 on success, 1 on failure
 */
int mat_outputf(
    const char *filename, const char *fmt, const double *A, int n, int m
);

/**
 * Output a matrix to a text file using the default format.
 *
 * @param filename name of the file to write to
 * @param A pointer to the flattened matrix data
 * @param n number of rows
 * @param m number of columns
 * @return 0 on success, 1 on failure
 */
int mat_output(const char *filename, const double *A, int n, int m);

#endif // IO_H
