#ifndef HELPERS_H
#define HELPERS_H

#include "optics_types.h"

/*
 * Compare two (distance, point-index) pairs using the OPTICS ordering.
 *
 * Returns true when (d1, idx1) precedes (d2, idx2), where:
 *   1. Smaller distance comes first.
 *   2. On ties, smaller original-dataset index comes first.
 */
bool less(const Point *points,
          const double d1, const int idx1,
          const double d2, const int idx2);

/*
 * Euclidean distance between two points.
 * Aborts if the points have different dimensionality.
 */
double euclidean_distance(const Point *a, const Point *b);

#endif
