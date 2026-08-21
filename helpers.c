#include "helpers.h"

#include <math.h>
#include <stdbool.h>

#define DIST_EPS 1e-12

bool less(const Point *points,
          const double d1, const int idx1,
          const double d2, const int idx2)
{
    if (fabs(d1 - d2) > DIST_EPS) return d1 < d2;
    return points[idx1].orig_idx < points[idx2].orig_idx;
}

double euclidean_distance_sq(const Point *a, const Point *b, const int dim)
{
    double dist_sq = 0.0;
    for (int i = 0; i < dim; i++)
    {
        double diff = a->coords[i] - b->coords[i];
        dist_sq += diff * diff;
    }
    return dist_sq;
}
