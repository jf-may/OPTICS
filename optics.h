/*
 * Ordering Points To Identify the Clustering Structure (OPTICS)
 *
 * OPTICS is a density-based clustering algorithm closely related to DBSCAN.
 * It identifies core samples and expands clusters from them, but—unlike
 * DBSCAN—preserves the full reachability hierarchy across a range of
 * neighborhood radii.
 *
 * This implementation uses brute-force O(n) neighbor search, making the
 * overall complexity O(n^2).  A spatial index (R*-tree, kd-tree) would
 * reduce this to O(n log n) on average and is the main planned improvement.
 *
 * See the run_optics() docstring below for the public interface.
 */

/*
 * Author: Juarez Ferriol May.
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef OPTICS_H
#define OPTICS_H

#include "optics_types.h"
#include "helpers.h"
#include "heap.h"

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
 *     neighborhood of the other.  Pass INFINITY to consider all points as
 *     potential neighbors (full hierarchy, slowest).
 *
 * minPts : int
 *     Minimum number of points (including the query point itself) required
 *     in an eps-neighborhood for a point to be classified as a core point.
 *
 * ordering_out : int **
 *     Output: dynamically allocated array of length size holding point
 *     indices in OPTICS traversal order.  Caller must free().
 *
 * reach_out : double **
 *     Output: reachability distances parallel to ordering_out.  The first
 *     point of each new cluster run has reach_distance == INFINITY.
 *     Caller must free().
 *
 * core_out : double **
 *     Output: core distances parallel to ordering_out.  Non-core points
 *     have INFINITY.  Caller must free().
 *
 * References
 * ----------
 * [1] Ankerst et al., "OPTICS: Ordering Points to Identify the Clustering
 *     Structure", SIGMOD 1999. https://doi.org/10.1145/304181.304187
 * [2] Scikit-learn OPTICS documentation and source.
 */
OpticsResult run_optics(Point *points, const int size, const double eps,
                        const int minPts);

/*
 * free_optics_result
 * ------------------
 * Releases the memory arrays allocated inside OpticsResult.
 */
void free_optics_result(OpticsResult *res);

#endif
