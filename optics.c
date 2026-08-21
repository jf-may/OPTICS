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
 * Performs an R-Tree range query to collect all points within 'epsilon' of
 * point 'p_idx'. Results are sorted by OPTICS ordering rules and stored
 * in the caller-supplied output arrays.
 *
 * Returns the number of valid neighbors found.
 */
static int get_neighbors(Point *points, const int p_idx, const double epsilon,
                         const RTreeNode *root, int *neighbors_out,
                         double *distances_out)
{
    if (!root) return 0;

    Point *p = &points[p_idx];
    int count = 0;
    double eps_sq = epsilon * epsilon;

    int stack_capacity = 256;
    const RTreeNode **stack = malloc(stack_capacity * sizeof(RTreeNode *));

    if (!stack)
    {
        fprintf(stderr, "OPTICS Error: Initial stack allocation failed in "
                        "get_neighbors.\n");
        exit(EXIT_FAILURE);
    }

    int top = 0;

    stack[top++] = root;

    while (top > 0)
    {
        const RTreeNode *node = stack[--top];

        if (node->is_leaf)
        {
            for (int i = 0; i < node->count; i++)
            {
                int node_p_idx = node->point_indices[i];
                double dist_sq = euclidean_distance_sq(p, &points[node_p_idx]);

                if (dist_sq <= eps_sq)
                {
                    neighbors_out[count] = node_p_idx;
                    distances_out[count] = sqrt(dist_sq);
                    count++;
                }
            }
        }
        else
        {
            for (int i = 0; i < node->count; i++)
            {
                double box_dist_sq = min_dist_sq_point_box(p, &node->box[i]);

                if (box_dist_sq <= eps_sq)
                {
                    if (top >= stack_capacity)
                    {
                        stack_capacity *= 2;
                        const RTreeNode **new_stack = realloc(stack,
                                                              stack_capacity *
                                                              sizeof(RTreeNode *));

                        if (!new_stack)
                        {
                            fprintf(stderr, "OPTICS Error: Stack reallocation "
                                            "failed in get_neighbors.\n");
                            free(stack);
                            exit(EXIT_FAILURE);
                        }
                        stack = new_stack;
                    }

                    stack[top++] = node->children[i];
                }
            }
        }
    }

    free(stack);

    sort_neighbors(points, neighbors_out, distances_out, count);
    return count;
}

/* ============================================================================
 * Internal: OPTICS core steps
 * ========================================================================= */

/*
 * compute_core_distance
 *
 * The distance to the min_pts-th nearest neighbor (p counts as
 * its own neighbor at distance 0).
 *
 * Returns INFINITY if p has fewer than min_pts neighbors within epsilon.
 */
static double compute_core_distance(Point *points, const int p_idx,
                                    const double epsilon, const int min_pts,
                                    const RTreeNode *root, int *neighbors_out,
                                    double *distances_out)
{
    int count = get_neighbors(points, p_idx, epsilon, root, neighbors_out,
                              distances_out);
    if (count < min_pts) return INFINITY;
    return distances_out[min_pts - 1];
}

/*
 * update_seeds
 *
 * After processing a core point, we examine all its neighbors. If a neighbor
 * hasn't been processed yet, we calculate its reachability from this core
 * point. If this new reachability is smaller than its current one, we update
 * it and promote it in the priority queue.
 */
static void update_seeds(Point *points, const int core_idx,
                         BinaryHeap *seeds, const double epsilon,
                         const RTreeNode *root, int *neighbors,
                         double *distances)
{
    Point *core      = &points[core_idx];
    double core_dist = core->core_distance;

    int count = get_neighbors(points, core_idx, epsilon, root, neighbors,
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

ClusterOrdering run_optics(Point *points, const int size, const double epsilon,
                           const int min_pts)
{
    if (!points || size < 1 || epsilon < 0.0 || min_pts < 1)
    {
        fprintf(stderr, "run_optics error: Invalid input parameters "
                        "provided.\n");
        exit(EXIT_FAILURE);
    }

    /* Initialize points. */
    for (int i = 0; i < size; i++)
    {
        points[i].processed      = false;
        points[i].core_distance  = INFINITY;
        points[i].reach_distance = INFINITY;
    }

    /* Build the spatial index (R-Tree) for O(n log n) neighborhood queries. */
    RTreeNode *rtree_root = build_rtree(points, size);

    /* Initialize the priority queue */
    BinaryHeap seeds;
    heap_init(&seeds, points, size);

    /* Allocate memory for intermediate and output arrays */
    int    *neighbors  = (int    *)malloc((size_t)size * sizeof(int));
    double *distances  = (double *)malloc((size_t)size * sizeof(double));
    int    *ordering   = (int    *)malloc((size_t)size * sizeof(int));
    double *core_vals  = (double *)malloc((size_t)size * sizeof(double));
    double *reach_vals = (double *)malloc((size_t)size * sizeof(double));

    if (!neighbors || !distances || !ordering || !core_vals || !reach_vals)
    {
        free(neighbors);
        free(distances);
        free(ordering);
        free(core_vals);
        free(reach_vals);
        free_rtree(rtree_root);
        heap_free(&seeds);

        fprintf(stderr, "run_optics error: Out of memory during array "
                        "initialization.\n");
        exit(EXIT_FAILURE);
    }

    /* Index of processed points */
    int ordering_idx = 0;

    /* Main OPTICS iteration: for each unprocessed point i */
    for (int i = 0; i < size; i++)
    {
        if (points[i].processed) continue;

        Point *p         = &points[i];
        p->processed     = true;
        p->core_distance = compute_core_distance(points, i, epsilon, min_pts,
                                                 rtree_root, neighbors,
                                                 distances);

        ordering[ordering_idx]   = i;
        reach_vals[ordering_idx] = p->reach_distance;
        core_vals[ordering_idx]  = p->core_distance;
        ordering_idx++;

        /* If the point processed is a core point, we expand the cluster */
        if (p->core_distance != INFINITY)
        {
            update_seeds(points, i, &seeds, epsilon, rtree_root, neighbors,
                         distances);

            /* Exhaust the priority queue */
            while (!heap_is_empty(&seeds))
            {
                int    closest_idx = heap_extract_min(&seeds);
                Point *closest     = &points[closest_idx];

                closest->processed     = true;
                closest->core_distance = compute_core_distance(points,
                                                               closest_idx,
                                                               epsilon,
                                                               min_pts,
                                                               rtree_root,
                                                               neighbors,
                                                               distances);

                ordering[ordering_idx]   = closest_idx;
                reach_vals[ordering_idx] = closest->reach_distance;
                core_vals[ordering_idx]  = closest->core_distance;
                ordering_idx++;

                /*
                 * If the last point processed is also a core point, we find
                 * its neighbors. If any neighbor has a shorter reachability
                 * path, it its promoted.
                 */
                if (closest->core_distance != INFINITY)
                {
                    update_seeds(points, closest_idx, &seeds, epsilon,
                                 rtree_root, neighbors, distances);
                }
            }
        }
    }

    free(neighbors);
    free(distances);
    free_rtree(rtree_root);
    heap_free(&seeds);

    ClusterOrdering result;
    result.ordering = ordering;
    result.reach    = reach_vals;
    result.core     = core_vals;

    return result;
}

void save_cluster_ordering_to_csv(const ClusterOrdering *co, int size,
                                 const char *filename)
{
    if (!co || !co->ordering || !co->core || !co->reach)
    {
        fprintf(stderr, "save_cluster_ordering_to_csv error: Invalid "
                        "input pointer.\n");
        exit(EXIT_FAILURE);
    }

    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        fprintf(stderr, "save_cluster_ordering_to_csv error: Could not open "
                        "file '%s' for writing.\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Headers */
    fprintf(file, "ordering,core_distance,reachability_distance\n");

    for (int i = 0; i < size; i++)
    {
        fprintf(file, "%d,%f,%f\n", co->ordering[i], co->core[i],
                co->reach[i]);
    }

    fclose(file);
}

void free_cluster_ordering(ClusterOrdering *co)
{
    if (!co) return;

    free(co->ordering);
    free(co->reach);
    free(co->core);
}
