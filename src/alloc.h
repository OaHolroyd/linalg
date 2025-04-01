#ifndef ALLOC_H
#define ALLOC_H

/**
 * Allocates a 2D 'array' of doubles. Must be freed with free_2d.
 *
 * This actually allocates a 1D array of double-pointers into a contiguous block
 * of double-memory, so it is not an array in the C sense. This does mean that
 * you can index into it like a 2D array, e.g.
 *   arr[i][j] = 0.0;
 * It also means that arr[0] points to the first row AND the flattened array, so
 *   arr[0][m] == arr[1][0]
 * Although this second property is often convenient (e.g. for LAPACK functions
 * which expect a pointer to a contiguous block of memory), it does mean that
 * row-overflows are not detectable, unless it is the last row:
 *   arr[0][m] = 0.0; // will not error, but may be unexpected
 *   arr[n-1][m] = 0.0; // will error, as this is the last row
 *
 * If performance is absolutely critical, it may be faster to use a 1D array and
 * do the indexing manually:
 *   arr_1d = arr[0];
 *   arr_1d[i * m + j] = 0.0; // This is equivalent to arr[i][j]
 *
 * @param n number of rows
 * @param m number of columns
 * @return pointer to a 2D array of size n x m
 */
double **malloc_d2d(int n, int m);

/**
 * Frees a 2D array of any type allocated with malloc_X2d.
 *
 * This frees both the pointers and the contiguous block of memory.
 *
 * @param arr pointer to the 2D array to free.
 */
#define free_2d(arr) internal_free_2d((void **)arr)
void internal_free_2d(void **arr);

#endif // ALLOC_H
