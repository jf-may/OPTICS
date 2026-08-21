#ifndef RTREE_H
#define RTREE_H

#include <stdbool.h>
#include "optics_types.h"

#define MAX_ENTRIES 32

/*
 * BoundingBox
 * -----------
 * Represents an Axis-Aligned Bounding Box (AABB) in N-dimensional space.
 *
 * Fields
 * ------
 * min_coords:
 *     Pointer to an array of size 'dim' holding the minimum bounds across all axes.
 * max_coords:
 *     Pointer to an array of size 'dim' holding the maximum bounds across all axes.
 */
typedef struct
{
    double *min_coords;
    double *max_coords;
} BoundingBox;

/*
 * RTreeNode
 * ---------
 * A node within the R-tree spatial index hierarchy. Used to accelerate
 * epsilon-neighborhood range queries from O(n) to average O(log n).
 *
 * Fields
 * ------
 * is_leaf:
 *     True if this node is at the bottom of the tree and contains direct
 *     references to data points. False if it is an internal node.
 * count:
 *     The number of valid entries currently stored in this node (from 0 up
 *     to MAX_ENTRIES).
 * coords_block:
 *     A single contiguous block of memory holding all min/max coordinates for
 *     every bounding box in this node. Prevents severe memory fragmentation.
 * box:
 *     Array of BoundingBoxes. For internal nodes, box[i] bounds the corresponding
 *     child node. For leaf nodes, box[i] bounds the corresponding data point.
 * children:
 *     (Union field) Array of pointers to sub-nodes. Only safely accessible
 *     when is_leaf is false.
 * point_indices:
 *     (Union field) Array of integers referencing the original index of points
 *     in the main dataset. Only safely accessible when is_leaf is true.
 */
typedef struct RTreeNode
{
    bool is_leaf;
    int count;
    double *coords_block;
    BoundingBox box[MAX_ENTRIES];
    union
    {
        struct RTreeNode *children[MAX_ENTRIES];
        int point_indices[MAX_ENTRIES];
    };
} RTreeNode;

/* Computes the minimum squared distance from a point to a bounding box */
double min_dist_sq_point_box(const Point *p, const BoundingBox *b,
                             const int dim);

/* Builds a perfectly balanced R-Tree using a Top-Down Bulk Loading strategy */
RTreeNode* build_rtree(Point *points, const int size, const int dim);

/* Recursively frees the R-Tree and all associated memory blocks */
void free_rtree(RTreeNode *node);

#endif
