#include "heap.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * Internal helpers
 * ========================================================================= */

static void heap_swap(BinaryHeap *h, int i, int j)
{
    int p_i = h->arr[i];
    int p_j = h->arr[j];

    h->arr[i] = p_j;
    h->arr[j] = p_i;

    h->pos[p_i] = j;
    h->pos[p_j] = i;
}

/*
 * heap_sift_up
 *
 * Restore the min-heap property upward from position i after an insertion or
 * key decrease.  Compares using reach_distance with orig_idx as the
 * tiebreaker (via the less() helper).
 */
static void heap_sift_up(BinaryHeap *h, int i)
{
    while (i > 0)
    {
        int parent = (i - 1) / 2;
        int ci = h->arr[i];
        int pi = h->arr[parent];

        if (less(h->points,
                 h->points[ci].reach_distance, ci,
                 h->points[pi].reach_distance, pi))
        {
            heap_swap(h, i, parent);
            i = parent;
        }
        else
        {
            break;
        }
    }
}

/*
 * heap_sift_down
 *
 * Restore the min-heap property downward from position i after the root is
 * replaced during extraction.
 */
static void heap_sift_down(BinaryHeap *h, int i)
{
    for (;;)
    {
        int smallest = i;
        int left     = 2 * i + 1;
        int right    = 2 * i + 2;

        if (left < h->size)
        {
            int sl = h->arr[smallest];
            int li = h->arr[left];
            if (less(h->points,
                     h->points[li].reach_distance, li,
                     h->points[sl].reach_distance, sl))
                smallest = left;
        }

        if (right < h->size)
        {
            int sm = h->arr[smallest];
            int ri = h->arr[right];
            if (less(h->points,
                     h->points[ri].reach_distance, ri,
                     h->points[sm].reach_distance, sm))
                smallest = right;
        }

        if (smallest == i) break;

        heap_swap(h, i, smallest);
        i = smallest;
    }
}

/* ============================================================================
 * Public API
 * ========================================================================= */

void heap_init(BinaryHeap *h, Point *pts, int n)
{
    h->size     = 0;
    h->capacity = n + 1;
    h->points   = pts;
    h->arr      = (int *)malloc(h->capacity * sizeof(int));
    h->pos      = (int *)malloc((size_t)n * sizeof(int));

    if (!h->arr || !h->pos)
    {
        free(h->arr);
        free(h->pos);
        fprintf(stderr, "heap_init: out of memory\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++)
        h->pos[i] = -1;
}

bool heap_is_empty(const BinaryHeap *h)
{
    return h->size == 0;
}

/*
 * heap_insert
 *
 * Add point p_idx to the heap.  The point's reach_distance must already be
 * set in h->points[p_idx] before calling.  p_idx must not already be present
 * (check h->pos[p_idx] == -1 beforehand).
 */
void heap_insert(BinaryHeap *h, int p_idx)
{
    int i = h->size;
    h->arr[i]     = p_idx;
    h->pos[p_idx] = i;
    h->size++;
    heap_sift_up(h, i);
}

/*
 * heap_decrease_key
 *
 * Update the heap after a point's reach_distance has been lowered.
 * The caller must write the new value to h->points[p_idx].reach_distance
 * BEFORE calling this function; new_reach is provided only for symmetry and
 * is not used directly here.
 */
void heap_decrease_key(BinaryHeap *h, int p_idx, double new_reach)
{
    (void)new_reach;          /* already committed to points[] by caller */
    heap_sift_up(h, h->pos[p_idx]);
}

/*
 * heap_extract_min
 *
 * Remove and return the index of the point with the smallest reach_distance
 * (ties broken by orig_idx).  The returned point's pos entry is set to -1.
 */
int heap_extract_min(BinaryHeap *h)
{
    int min_p = h->arr[0];

    h->size--;
    if (h->size > 0)
    {
        h->arr[0]          = h->arr[h->size];
        h->pos[h->arr[0]]  = 0;
        heap_sift_down(h, 0);
    }

    h->pos[min_p] = -1;
    return min_p;
}

void heap_free(BinaryHeap *h)
{
    free(h->arr);
    free(h->pos);
    h->arr = NULL;
    h->pos = NULL;
}
