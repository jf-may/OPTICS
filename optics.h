/*
 * Ordering Points To Identify the Clustering Structure (OPTICS)
 *
 * OPTICS is a density-based clustering algorithm closely related to DBSCAN.
 * It identifies core samples and expands clusters from them, but unlike
 * DBSCAN, it preserves the full reachability hierarchy across a range of
 * neighborhood radii.
 *
 * This implementation uses an R-Tree spatial index to accelerate neighborhood
 * searches, bringing the overall time complexity down from O(n^2) to an average
 * of O(n log n) for spatial data.
 *
 * Author: Juarez Ferriol May.
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef OPTICS_H
#define OPTICS_H

#include "optics_types.h"
#include "helpers.h"
#include "heap.h"
#include "rtree.h"

/*
 * run_optics
 * ----------
 * Perform OPTICS clustering and extract an ordered list of points together
 * with their reachability and core distances.
 *
 * Parameters
 * ----------
 * points : Point *
 *     Input array of n points.  coords, dim, and orig_idx must be populated
 *     by the caller; processed, core_distance, and reach_distance are
 *     initialised internally.
 *
 * size : int
 *     Number of points.
 *
 * eps : double
 *     Maximum distance between two samples for one to be considered in the
 *     neighborhood of the other. Pass INFINITY to consider all points as
 *     potential neighbors.
 *
 * minPts : int
 *     Minimum number of points (including the query point itself) required
 *     in an eps-neighborhood for a point to be classified as a core point.
 *
 * Returns
 * -------
 * ClusterOrdering : struct
 *     Contains dynamically allocated arrays (ordering, core, reach).
 *     The caller assumes ownership of these arrays and must free them
 *     using free_cluster_ordering().
 *
 * References
 * ----------
 * [1] Ankerst et al., "OPTICS: Ordering Points to Identify the Clustering
 *     Structure", SIGMOD 1999. https://doi.org/10.1145/304181.304187
 *
 * [2] Scikit-learn OPTICS documentation and source.
 */
ClusterOrdering run_optics(Point *points, const int size, const double eps,
                           const int minPts);

/*
 * free_cluster_ordering
 * ---------------------
 * Safely releases the memory arrays dynamically allocated inside the
 * ClusterOrdering result struct.
 */
void free_cluster_ordering(ClusterOrdering *res);

#endif
