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
 * dim:
 *     Number of coordinates / feature dimensions.
 * orig_idx:
 *     Original position of this point in the input array.
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
    int     dim;
    int     orig_idx;
    bool    processed;
    double  core_distance;
    double  reach_distance;
} Point;

typedef struct
{
    int *ordering;
    double *core;
    double *reach;
} OpticsResult;

#endif
