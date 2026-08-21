#ifndef HEAP_H
#define HEAP_H

#include <stdbool.h>

#include "optics_types.h"
#include "helpers.h"

/*
 * BinaryHeap
 * ----------
 * Min-heap of point indices ordered by reach_distance (ties broken by
 * orig_idx via the less() comparator).  Used as the OPTICS seed set.
 *
 * Fields
 * ------
 * arr:      Array of point indices stored in heap order.
 * size:     Current number of elements in the heap.
 * capacity: Allocated length of arr (== total number of points + 1).
 * points:   Pointer into the shared Point array (for reach_distance lookup).
 * pos:      pos[p_idx] = current heap position of point p_idx, -1 if absent.
 */
typedef struct
{
    int    *arr;
    int     size;
    int     capacity;
    Point  *points;
    int    *pos;
} BinaryHeap;

void heap_init        (BinaryHeap *h, Point *pts, int n);
bool heap_is_empty    (const BinaryHeap *h);
void heap_insert      (BinaryHeap *h, int p_idx);
void heap_decrease_key(BinaryHeap *h, int p_idx);
int  heap_extract_min (BinaryHeap *h);
void heap_free        (BinaryHeap *h);

#endif
