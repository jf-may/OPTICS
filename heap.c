#include "heap.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * Internal helpers
 * ========================================================================= */

static void heap_swap(BinaryHeap *heap, int i, int j)
{
    int p_i = heap->arr[i];
    int p_j = heap->arr[j];

    heap->arr[i] = p_j;
    heap->arr[j] = p_i;

    heap->pos[p_i] = j;
    heap->pos[p_j] = i;
}

/*
 * heap_sift_up
 *
 * Restore the min-heap property upward from position i after an insertion or
 * key decrease. Compares using reach_distance.
 */
static void heap_sift_up(BinaryHeap *heap, int i)
{
    while (i > 0)
    {
        int parent = (i - 1) / 2;
        int children_idx = heap->arr[i];
        int parent_idx = heap->arr[parent];

        if (less(heap->points[children_idx].reach_distance, children_idx,
                 heap->points[parent_idx].reach_distance, parent_idx))
        {
            heap_swap(heap, i, parent);
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
static void heap_sift_down(BinaryHeap *heap, int i)
{
    for (;;)
    {
        int smallest = i;
        int left     = 2 * i + 1;
        int right    = 2 * i + 2;

        if (left < heap->size)
        {
            int smallest_idx = heap->arr[smallest];
            int left_idx = heap->arr[left];
            if (less(heap->points[left_idx].reach_distance, left_idx,
                     heap->points[smallest_idx].reach_distance, smallest_idx))
                smallest = left;
        }

        if (right < heap->size)
        {
            int smallest_idx = heap->arr[smallest];
            int right_idx = heap->arr[right];
            if (less(heap->points[right_idx].reach_distance, right_idx,
                     heap->points[smallest_idx].reach_distance, smallest_idx))
                smallest = right;
        }

        if (smallest == i) break;

        heap_swap(heap, i, smallest);
        i = smallest;
    }
}

/* ============================================================================
 * Public API
 * ========================================================================= */

void heap_init(BinaryHeap *heap, Point *pts, int n)
{
    heap->size     = 0;
    heap->capacity = n + 1;
    heap->points   = pts;
    heap->arr      = (int *)malloc(heap->capacity * sizeof(int));
    heap->pos      = (int *)malloc((size_t)n * sizeof(int));

    if (!heap->arr || !heap->pos)
    {
        free(heap->arr);
        free(heap->pos);
        fprintf(stderr, "heap_init: out of memory\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++)
        heap->pos[i] = -1;
}

bool heap_is_empty(const BinaryHeap *heap)
{
    return heap->size == 0;
}

/*
 * heap_insert
 *
 * Add point p_idx to the heap.
 */
void heap_insert(BinaryHeap *heap, int p_idx)
{
    int i = heap->size;
    heap->arr[i]     = p_idx;
    heap->pos[p_idx] = i;
    heap->size++;
    heap_sift_up(heap, i);
}

/*
 * heap_decrease_key
 *
 * Update the heap after a point's reach_distance has been lowered.
 */
void heap_decrease_key(BinaryHeap *heap, int p_idx)
{
    heap_sift_up(heap, heap->pos[p_idx]);
}

/*
 * heap_extract_min
 *
 * Remove and return the index of the point with the smallest reach_distance.
 * The returned point's pos entry is set to -1.
 */
int heap_extract_min(BinaryHeap *heap)
{
    int min_p = heap->arr[0];

    heap->size--;
    if (heap->size > 0)
    {
        heap->arr[0]            = heap->arr[heap->size];
        heap->pos[heap->arr[0]] = 0;
        heap_sift_down(heap, 0);
    }

    heap->pos[min_p] = -1;
    return min_p;
}

void heap_free(BinaryHeap *heap)
{
    free(heap->arr);
    free(heap->pos);
    heap->arr = NULL;
    heap->pos = NULL;
}
