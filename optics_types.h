#ifndef OPTICS_TYPES_H
#define OPTICS_TYPES_H

#include <stdbool.h>

/*
 * Point
 * -----
 * Internal representation of one input sample used by the OPTICS pipeline.
 *
 * Fields
 * ------
 * coords:
 *     Pointer to the first coordinate of this point inside the contiguous
 *     coordinate storage block.
 * processed:
 *     True once the point has been removed from the OPTICS seed structure and
 *     assigned to the final ordering.
 * core_distance:
 *     Core distance of the point. Set to INFINITY when undefined (fewer than
 *     minPts neighbors within eps).
 * reach_distance:
 *     Reachability distance of the point. Set to INFINITY when undefined.
 */
typedef struct
{
    double *coords;
    bool    processed;
    double  core_distance;
    double  reach_distance;
} Point;

/*
 * ClusterOrdering
 * ---------------
 * The final output of the OPTICS algorithm, representing the hierarchical
 * clustering structure as a set of parallel arrays.
 *
 * Fields
 * ------
 * ordering:
 *     Array of point indices representing the sequence in which points were
 *     processed. This sequence forms the x-axis of an OPTICS reachability plot.
 * core:
 *     Array of core distances for each point, ordered identically to the
 *     'ordering' array.
 * reach:
 *     Array of reachability distances for each point, ordered identically to
 *     the 'ordering' array. This forms the y-axis of a reachability plot.
 */
typedef struct
{
    int *ordering;
    double *core;
    double *reach;
} ClusterOrdering;

#endif
