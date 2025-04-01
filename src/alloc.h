#ifndef ALLOC_H
#define ALLOC_H

/**
 * Allocates a 2D array of doubles with the given dimensions, setting up row
 * pointers to the correct addresses. Must be freed with free_2d.
 *
 * @param ni number of rows
 * @param nj number of columns
 * @return pointer to the allocated 2D array
 */
double **malloc_d2d(int ni, int nj);

/**
 * Frees a 2D array of any type allocated with malloc_X2d.
 *
 * @param arr pointer to the 2D array to free.
 */
#define free_2d(arr) internal_free_2d((void **)arr)
void internal_free_2d(void **arr);

#endif // ALLOC_H
