#ifndef HELPERS_H
#define HELPERS_H

#include "optics_types.h"

/*
 * Compare two (distance, point-index) pairs using the OPTICS ordering.
 *
 * Returns true when (d1, idx1) precedes (d2, idx2), where:
 *   1. Smaller distance comes first.
 *   2. On ties, smaller index comes first.
 */
bool less(const double d1, const int idx1,
          const double d2, const int idx2);

/*
 * euclidean_distance_sq
 *
 * Calculates the squared distance between two points.
 */
double euclidean_distance_sq(const Point *a, const Point *b, const int dim);

#endif
