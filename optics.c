#include "optics.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * Internal: neighbor utilities
 * ========================================================================= */

/*
 * sort_neighbors
 *
 * Insertion-sort of a parallel (neigh_ids, distances) array by the OPTICS
 * ordering (ascending distance, ties broken by orig_idx).
 *
 * TODO: replace with a faster sort for large n.
 */
static void sort_neighbors(const Point *points, int *neigh_ids,
                           double *distances, const int count)
{
    for (int i = 1; i < count; i++)
    {
        int    key_id   = neigh_ids[i];
        double key_dist = distances[i];
        int    j        = i - 1;

        while (j >= 0 &&
               less(points, key_dist, key_id, distances[j], neigh_ids[j]))
        {
            neigh_ids[j + 1] = neigh_ids[j];
            distances[j + 1] = distances[j];
            j--;
        }

        neigh_ids[j + 1] = key_id;
        distances[j + 1] = key_dist;
    }
}

/*
 * get_neighbors
 *
 * Collect all points within epsilon of point p_idx (including p_idx itself at
 * distance 0), sort them by the OPTICS ordering, and store results in the
 * caller-supplied arrays.
 *
 * Returns the number of neighbors found.
 *
 * TODO: replace brute-force O(n) scan with an R*-tree or kd-tree for O(log n).
 */
static int get_neighbors(Point *points, const int size, const int p_idx,
                         const double epsilon, int *neighbors_out,
                         double *distances_out)
{
    Point *p = &points[p_idx];
    int count = 0;

    for (int i = 0; i < size; i++)
    {
        double dist = euclidean_distance(p, &points[i]);
        if (dist <= epsilon)
        {
            neighbors_out[count] = i;
            distances_out[count] = dist;
            count++;
        }
    }

    sort_neighbors(points, neighbors_out, distances_out, count);
    return count;
}

/* ============================================================================
 * Internal: OPTICS core steps
 * ========================================================================= */

/*
 * compute_core_distance
 *
 * OPTICS definition: the smallest epsilon' such that N_{epsilon'}(p) >= min_pts.
 * Equivalently, the distance to the min_pts-th nearest neighbor (p counts as
 * its own neighbor at distance 0).
 *
 * Returns INFINITY if p has fewer than min_pts neighbors within epsilon.
 */
static double compute_core_distance(Point *points, const int size, const int p_idx,
                                    const double epsilon, const int min_pts,
                                    int *neighbors_out, double *distances_out)
{
    int count = get_neighbors(points, size, p_idx, epsilon, neighbors_out,
                              distances_out);
    if (count < min_pts) return INFINITY;
    return distances_out[min_pts - 1];
}

/*
 * update_seeds
 *
 * After processing center_idx as a core point, update the seed heap with
 * improved reachability distances for all unprocessed neighbors.
 *
 * OPTICS reachability rule:
 *   reach_dist(o, p) = max(core_dist(p), dist(p, o))
 */
static void update_seeds(Point *points, const int size, const int center_idx,
                         BinaryHeap *seeds, const double epsilon,
                         int *neighbors, double *distances)
{
    Point *center    = &points[center_idx];
    double core_dist = center->core_distance;

    int count = get_neighbors(points, size, center_idx, epsilon, neighbors,
                              distances);

    for (int i = 0; i < count; i++)
    {
        int    n_idx   = neighbors[i];
        Point *n_point = &points[n_idx];

        if (n_point->processed) continue;

        double new_reach = fmax(core_dist, distances[i]);

        if (n_point->reach_distance == INFINITY)
        {
            n_point->reach_distance = new_reach;
            heap_insert(seeds, n_idx);
        }
        else if (new_reach < n_point->reach_distance)
        {
            n_point->reach_distance = new_reach;
            heap_decrease_key(seeds, n_idx, new_reach);
        }
    }
}

/* ============================================================================
 * Public API
 * ========================================================================= */

OpticsResult run_optics(Point *points, const int size, const double epsilon,
                        const int min_pts)
{
    if (!points || size < 1 || epsilon < 0.0 || min_pts < 1)
    {
        fprintf(stderr, "run_optics: invalid input parameters\n");
        exit(EXIT_FAILURE);
    }

    /*
     * Initialize points. Coords, dim, orig_idx and size come from csv parser
     * or manually. TODO: The user should only give a coords array either
     * through a csv or manually and we should get everything else.
     */
    for (int i = 0; i < size; i++)
    {
        points[i].processed      = false;
        points[i].core_distance  = INFINITY;
        points[i].reach_distance = INFINITY;
    }

    /*
     * Initialize heap. Not needed yet but possible error point so we call it
     * before allocating any more memory so we only have to think about freeing
     * the heap and not everything else if something goes wrong.
     */
    BinaryHeap seeds;
    heap_init(&seeds, points, size);

    /*
     * Allocate all the outputs (ordering, core_vals, reach_vals) and the
     * intermediate arrays (neighbors, distances) and check for correct
     * memory allocation.
     */
    int    *ordering  = (int    *)malloc((size_t)size * sizeof(int));
    double *core_vals = (double *)malloc((size_t)size * sizeof(double));
    double *reach_vals = (double *)malloc((size_t)size * sizeof(double));
    int    *neighbors = (int    *)malloc((size_t)size * sizeof(int));
    double *distances = (double *)malloc((size_t)size * sizeof(double));

    if (!ordering || !core_vals || !reach_vals || !neighbors || !distances)
    {
        free(ordering);
        free(core_vals);
        free(reach_vals);
        free(neighbors);
        free(distances);
        heap_free(&seeds);
        fprintf(stderr, "run_optics: out of memory\n");
        exit(EXIT_FAILURE);
    }

    /* Index of processed points */
    int order_idx = 0;

    /* For each unprocessed point i */
    for (int i = 0; i < size; i++)
    {
        if (points[i].processed) continue;

        /*
         * First pass: heap empty, no point processed yet. Compute everything
         * and append to the order.
         */
        Point *p    = &points[i];
        p->processed    = true;
        p->core_distance = compute_core_distance(points, size, i, epsilon, min_pts,
                                                  neighbors, distances);

        ordering[order_idx]   = i;
        reach_vals[order_idx] = p->reach_distance;
        core_vals[order_idx]  = p->core_distance;
        order_idx++;

        /*
         * If the last point processed was a core point:
         */
        if (p->core_distance != INFINITY)
        {
            /*
             * Then populate the heap with its neighbors, sorted from closest
             * to farthest.
             */
            update_seeds(points, size, i, &seeds, epsilon, neighbors, distances);

            /*
             * As long as the heap has points, process them first, from closest
             * to farthest.
             */
            while (!heap_is_empty(&seeds))
            {
                int    closest_idx = heap_extract_min(&seeds);
                Point *closest     = &points[closest_idx];

                closest->processed    = true;
                closest->core_distance = compute_core_distance(points, size,
                                                               closest_idx,
                                                               epsilon, min_pts,
                                                               neighbors,
                                                               distances);

                ordering[order_idx]   = closest_idx;
                reach_vals[order_idx] = closest->reach_distance;
                core_vals[order_idx]  = closest->core_distance;
                order_idx++;

                /*
                 * If the last point processed was a core point, populate the
                 * heap with its neighbors. TODO: I don't actually know how
                 * this interacts with the previous points in the heap.
                 */
                if (closest->core_distance != INFINITY)
                {
                    update_seeds(points, size, closest_idx, &seeds, epsilon,
                                 neighbors, distances);
                }
            }
        }
    }

    free(neighbors);
    free(distances);
    heap_free(&seeds);

    OpticsResult result;
    result.ordering = ordering;
    result.reach    = reach_vals;
    result.core     = core_vals;

    return result;
}

void free_optics_result(OpticsResult *res)
{
    if (!res) return;

    free(res->ordering);
    free(res->reach);
    free(res->core);

    res->ordering = NULL;
    res->reach    = NULL;
    res->core     = NULL;
}
