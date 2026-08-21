#include "rtree.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * Internal: R-Tree Builder Utilities
 * ========================================================================= */

/*
 * allocate_node_memory
 * --------------------
 * High-performance memory allocator for an R-Tree node. Instead of allocating
 * memory for every individual bounding box (which causes severe fragmentation),
 * it allocates one contiguous block of memory for the entire node, and maps
 * the internal BoundingBox pointers to their correct offsets within the block.
 */
static void allocate_node_memory(RTreeNode *node, int dim)
{
    int total_doubles = MAX_ENTRIES * 2 * dim;
    node->coords_block = (double *)malloc(total_doubles * sizeof(double));

    if (!node->coords_block)
    {
        fprintf(stderr, "R-Tree Error: Memory allocation failed for node "
                        "coords_block.\n");
        exit(EXIT_FAILURE);
    }

    // Map the pointers of each BoundingBox directly into the contiguous block
    for (int i = 0; i < MAX_ENTRIES; i++)
    {
        // min_coords starts at the beginning of this entry's designated space
        node->box[i].min_coords = node->coords_block + (i * 2 * dim);

        // max_coords starts immediately after min_coords
        node->box[i].max_coords = node->coords_block + (i * 2 * dim) + dim;
    }
}

/*
 * quicksort_indices
 * -----------------
 * A robust, production-grade QuickSort implementation for spatially sorting points.
 * Uses Hoare's partition, middle-element pivot selection, and tail-call
 * optimization to completely eliminate the risk of stack overflow.
 */
static void quicksort_indices(Point *points, int *indices, int left_limit,
                              int right_limit, int axis)
{
    while (left_limit < right_limit)
    {
        // Pivot selection
        int mid = left_limit + (right_limit - left_limit) / 2;
        double pivot = points[indices[mid]].coords[axis];

        // Hoare partition scheme
        int i = left_limit - 1;
        int j = right_limit + 1;

        while (true)
        {
            do i++; while (points[indices[i]].coords[axis] < pivot);
            do j--; while (points[indices[j]].coords[axis] > pivot);

            if (i >= j) break;

            int tmp = indices[i];
            indices[i] = indices[j];
            indices[j] = tmp;
        }

        // Always recurse on the smaller partition, loop on the larger one.
        if ((j - left_limit) < (right_limit - (j + 1)))
        {
            quicksort_indices(points, indices, left_limit, j, axis);
            left_limit = j + 1;
        }
        else
        {
            quicksort_indices(points, indices, j + 1, right_limit, axis);
            right_limit = j;
        }
    }
}

/*
 * build_rtree_recursive
 * ---------------------
 * Recursively divides the dataset into spatial chunks to build a perfectly
 * balanced R-Tree. Returns a pointer to the created node.
 */
static RTreeNode *build_rtree_recursive(Point *points, int size, int dim,
                                        int *indices)
{
    // Allocate the node structure
    RTreeNode *node = (RTreeNode *)calloc(1, sizeof(RTreeNode));
    if (!node)
    {
        fprintf(stderr, "R-Tree Error: Memory allocation failed for "
                        "RTreeNode.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate the contiguous coordinate memory block for this node
    allocate_node_memory(node, dim);

    /* ------------------------------------------------------------------------
     * BASE CASE: Small enough to be a leaf node
     * ----------------------------------------------------------------------*/
    if (size <= MAX_ENTRIES)
    {
        node->is_leaf = true;
        node->count = size;

        for (int i = 0; i < size; i++)
        {
            int p_idx = indices[i];
            node->point_indices[i] = p_idx;

            // Define the bounding box of a single point (min == max)
            for (int d = 0; d < dim; d++)
            {
                node->box[i].min_coords[d] = points[p_idx].coords[d];
                node->box[i].max_coords[d] = points[p_idx].coords[d];
            }
        }
        return node;
    }

    /* ------------------------------------------------------------------------
     * RECURSIVE CASE: Internal node containing smaller child nodes
     * ----------------------------------------------------------------------*/
    node->is_leaf = false;

    // Find the axis with the widest spatial spread
    int best_axis = 0;
    double max_spread = -1.0;

    for (int d = 0; d < dim; d++)
    {
        double min_val = points[indices[0]].coords[d];
        double max_val = min_val;

        for (int i = 1; i < size; i++)
        {
            double val = points[indices[i]].coords[d];
            if (val < min_val) min_val = val;
            if (val > max_val) max_val = val;
        }

        double spread = max_val - min_val;
        if (spread > max_spread)
        {
            max_spread = spread;
            best_axis = d;
        }
    }

    // Sort the points along the widest axis
    quicksort_indices(points, indices, 0, size - 1, best_axis);

    // Divide remaining points into children.
    // Adding (MAX_ENTRIES - 1) before dividing ensures we round up.
    int chunk_size = (size + MAX_ENTRIES - 1) / MAX_ENTRIES;

    int child_idx = 0;
    for (int i = 0; i < size; i += chunk_size)
    {
        int current_chunk_size = chunk_size;
        if (i + current_chunk_size > size) current_chunk_size = size - i;

        // Build the child recursively
        RTreeNode *child = build_rtree_recursive(points, current_chunk_size,
                                                 dim, &indices[i]);
        node->children[child_idx] = child;

        // Wrap this new child in a tightly fitting bounding box
        for (int d = 0; d < dim; d++)
        {
            // Seed the bounds with the child's first box
            node->box[child_idx].min_coords[d] = child->box[0].min_coords[d];
            node->box[child_idx].max_coords[d] = child->box[0].max_coords[d];

            // Expand the bounds to encompass all remaining boxes in the child
            for (int k = 1; k < child->count; k++)
            {
                if (child->box[k].min_coords[d] < node->box[child_idx].min_coords[d])
                    node->box[child_idx].min_coords[d] = child->box[k].min_coords[d];

                if (child->box[k].max_coords[d] > node->box[child_idx].max_coords[d])
                    node->box[child_idx].max_coords[d] = child->box[k].max_coords[d];
            }
        }
        child_idx++;
    }

    node->count = child_idx;
    return node;
}

/* ============================================================================
 * Public Helpers
 * ========================================================================= */

/*
 * min_dist_sq_point_box
 * ---------------------
 * Calculates the shortest squared distance from a point to an AABB.
 * Returns 0 if the point is physically inside the box along all dimensions.
 */
double min_dist_sq_point_box(const Point *p, const BoundingBox *b,
                             const int dim)
{
    double dist_sq = 0.0;

    for (int i = 0; i < dim; i++)
    {
        if (p->coords[i] < b->min_coords[i])
        {
            double diff = b->min_coords[i] - p->coords[i];
            dist_sq += diff * diff;
        }
        else if (p->coords[i] > b->max_coords[i])
        {
            double diff = p->coords[i] - b->max_coords[i];
            dist_sq += diff * diff;
        }
    }

    return dist_sq;
}

/* ============================================================================
 * Public API: R-Tree
 * ========================================================================= */

/*
 * build_rtree
 * -----------
 * Public wrapper to kick off the Top-Down bulk loading process.
 * Validates inputs and prepares the index array.
 */
RTreeNode* build_rtree(Point *points, const int size, const int dim)
{
    if (!points || size < 1)
    {
        fprintf(stderr, "R-Tree Error: Invalid point array or size "
                        "provided.\n");
        return NULL;
    }

    int *indices = (int *)malloc((size_t)size * sizeof(int));
    if (!indices)
    {
        fprintf(stderr, "R-Tree Error: Memory allocation failed for indices "
                        "array.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < size; i++) indices[i] = i;

    RTreeNode *root = build_rtree_recursive(points, size, dim, indices);

    free(indices);
    return root;
}

/*
 * free_rtree
 * ----------
 * Recursively traverses the tree and frees all allocated memory.
 */
void free_rtree(RTreeNode *node)
{
    if (!node) return;

    if (!node->is_leaf)
    {
        for (int i = 0; i < node->count; i++) free_rtree(node->children[i]);
    }

    if (node->coords_block) free(node->coords_block);

    free(node);
}
