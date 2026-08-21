#ifndef DOPTICS_H
#define DOPTICS_H

#include "optics.h"

/*
 * Distributes the dataset using a k-d tree spatial partitioning strategy.
 * The master node (rank 0) partitions the data and scatters it to all nodes.
 */
void doptics_scatter_data(Point *global_points, const int global_size,
                          const int dim, const int rank, const int num_procs,
                          Point **local_points_out, int *local_size_out);

/*
 * Merges two cluster orderings across partition boundaries.
 */
ClusterOrdering doptics_merge(Point *pts1, int n1, ClusterOrdering *co1,
                              Point *pts2, int n2, ClusterOrdering *co2,
                              double epsilon, int min_pts, int dim);

#endif
