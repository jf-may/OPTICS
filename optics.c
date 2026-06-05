/* 
    OPTICS - Ordering Points To Identify the Clustering Structure
    Juarez Ferriol May. 2026.
    Usage: ./optics input.csv epsilon minPts output.csv
    Input CSV: one point per line, comma-separated doubles, no header,
        comments with #
    Output CSV: order_idx, orig_point_idx, reachability_dist, core_dist
        UNDEFINED distances are output as -1.0
*/

#include <stdio.h>      // printf, fprintf, stderr, fopen
#include <stdlib.h>     // exit
#include <string.h>     // strtok
#include <stdbool.h>    // bool
#include <math.h>       // INFINITY, fmax

#define UNDEFINED (-1.0) // Undefined for CSV output
#define INF (INFINITY)

typedef struct
{
    double *coords;         // Value array of dimension size
    int dim;                // Dimension
    int orig_idx;           // Initial index in input file
    bool processed;
    double core_distance;
    double reach_distance;
} Point;

typedef struct
{
    int *arr;       // Array of point indices
    int size;       // Current array size
    int capacity;   // Current allocated size
    Point *points;  // Access to reach_distance and orig_idx
    int *pos;       // pos[p_idx] = current heap index, -1 if absent
} BinaryHeap;

Point *load_csv(const char *filename, int *out_n, int *out_dim)
{
    // Open file
    FILE *f = fopen(filename, "r");
    if (!f)
    {
        fprintf(stderr, "Cannot open %s\n", filename);
        exit(1);
    }

    // First pass: count points and determine dimension
    char line[4096]; // POSSIBLE POINT OF FAILURE. LINE MAY BE BIGGER THAN 4096 CHARACTERS
    int n = 0;
    int dim = -1;

    // Read file in newline-terminated string chunks
    while (fgets(line, sizeof(line), f))
    {
        // Ignore empty lines and comments
        if (line[0] == "\n" || line[0] == "#") continue;

        // On start of data line, count one point
        n++;

        // On first data line, determine dimension by counting delimiters
        if (dim == -1)
        {
            // Divide line into tokens delimited by comma
            char *token = strtok(line, ",");
            
            dim = 0;
            while (token)
            {
                dim++;
                token = strtok(NULL, ","); // Continues from where it left off
            }
        }
    }

    // Go back to the start of the file
    rewind(f);

    if (n == 0 || dim <= 0)
    {
        fprintf(stderr, "Empty or invalid CSV\n");
        exit(1);
    }

    // Allocate memory
    double *coords_block = (double *)malloc(n * dim * sizeof(double));
    Point *points = (Point *)malloc(n * sizeof(Point));

    // Read file in newline-terminated string chunks
    int i = 0;
    while (fgets(line, sizeof(line), f))
    {
        // Ignore empty lines and comments
        if (line[0] == "\n" || line[0] == "#") continue;

        // Divide line into tokens delimited by comma
        char *token = strtok(line, ",");

        int j = 0;
        while (token && j < dim)
        {
            coords_block[i * dim + j] = atof(token); // Write values to coords
            token = strtok(NULL, ","); // Continues from where it left off
            j++;
        }

        if (j != dim)
        {
            fprintf(stderr, "Row %d has wrong number of columns\n", i);
            exit(1);
        }

        // points.coords[i] points to the i-th row of coords_block
        points[i].coords = coords_block + i * dim;
        points[i].dim = dim;
        points[i].orig_idx = i;
        i++;
    }

    // Close input file
    fclose(f);

    *out_n = n;
    *out_dim = dim;
    return points;
}

// Eucledian distance for arbitrary dimension
static double euclid_distance(const Point *a, const Point *b)
{
    double sum = 0.0;
    for (int j = 0; j < a->dim; j++)
    {
        double d = a->coords[j] - b->coords[j];
        sum += d * d;
    }
    return sqrt(sum);
}

static int get_neighbors(Point *points, int n, int p_idx, double eps, int *neigh_ids, double *distances)
{
    Point *p = &points[p_idx];
    int count = 0;

    // Brute force, go through every point
    for (int i = 0; i < n; i++)
    {
        double distance = euclid_distance(p, &points[i]);

        // If it is in the epsilon radius, add it as neighbor
        if (distance <= eps)
        {
            neigh_ids[count] = i;
            distances[count] = distance;
            count++;
        }
    }

    // Sort distances

    return count;
}

// Returns INF if not core, else the minPts-th smallest distance
static double compute_core_distance(Point *points, int n, int p_idx, double eps, int minPts)
{
    int *neighbors = (int *)malloc(n * sizeof(int));
    double *distances = (double *)malloc(n * sizeof(double));

    // Returns sorted neighbor ids and distances in an epsilon radius
    int count = get_neighbors(points, n, p_idx, eps, &neighbors, &distances);

    // If there are less than minPts neighbors, return INF
    if (count < minPts)
    {
        free(distances);
        return INF;
    }

    // Returns distance to minPts closest neighbor
    double core_distance = distances[minPts - 1];
    free(distances);
    return core_distance;
}

static void heap_init(BinaryHeap *h, Point *pts, int n)
{
    h->size = 0;
    h->capacity = n + 1;
    h->points = pts;
    h->arr = (int *)malloc(h->capacity * sizeof(int));
    h->pos = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) h->pos[i] = -1;
}

static bool heap_is_empty(const BinaryHeap *h)
{
    return h->size == 0;
}

static void heap_free(BinaryHeap *h)
{
    free(h->arr);
    free(h->pos);
}

static void heap_swap(BinaryHeap *h, int i, int j)
{
    int p_i = h->arr[i];
    int p_j = h->arr[j];
    h->arr[i] = p_j;
    h->arr[j] = p_i;
    h->pos[p_i] = j;
    h->pos[p_j] = i;
}

static void heap_bubble_up(BinaryHeap *h, int idx)
{
    while (idx > 0)
    {
        int parent = (idx - 1) / 2;
        int p_idx = h->arr[idx];
        int parent_idx = h->arr[parent];
        double p_reach = h->points[p_idx].reach_distance;
        double parent_reach = h->points[parent_idx].reach_distance;
        // Not necessary unless second case happens.
        // Definition made for code clarity. Probably should just not do it
        int p_oidx = h->points[p_idx].orig_idx;
        int parent_oidx = h->points[parent_idx].orig_idx;

        if (p_reach < parent_reach ||
            (p_reach == parent_reach && p_oidx < parent_oidx))
        {
            heap_swap(h, idx, parent);
            idx = parent;
        }
        else break;
    }
}

static void heap_bubble_down()
{

}

// Insert a point into the heap
static void heap_insert(BinaryHeap *h, int p_idx)
{
    // If the array is full, duplicate its size
    if (h->size >= h->capacity)
    {
        h->capacity *= 2;
        h->arr = (int *)realloc(h->arr, h->capacity * sizeof(int));
    }

    // Insert point last and bubble up until necessary
    h->arr[h->size] = p_idx;
    h->pos[p_idx] = h->size;
    heap_bubble_up(h, h->size);
    h->size++;
}

// Improve a point's reachability
static void heap_decrease_key(BinaryHeap *h, int p_idx, double new_reach)
{
    if (h->pos[p_idx] == -1) return; // Not in heap
    
    Point *p = &h->points[p_idx];

    // If the new reachability is better, update and bubble up
    if (new_reach < p->reach_distance)
    {
        p->reach_distance = new_reach;
        heap_bubble_up(h, h->pos[p_idx]);
    }
}

static int heap_extract_min(BinaryHeap *h)
{
    if (h->size == 0) return -1;

    // Save the smallest index and remove that element in the heap
    int min_idx = h->arr[0];
    h->pos[min_idx] = -1;
    h->size--;

    // Replace it with the last element in the heap and bubble it down
    if (h->size > 0)
    {
        h->arr[0] = h->arr[h->size];
        h->pos[h->arr[0]] = 0;
        heap_bubble_down(h, 0);
    }

    return min_idx;
}

static void update_seeds(Point *points, int n, int center_idx, BinaryHeap *seeds,
                  double eps, int minPts)
{
    Point *center = &points[center_idx];
    double core_dist = center->core_distance;

    int *neighbors = (int *)malloc(n * sizeof(int));
    double *distances = (double *)malloc(n * sizeof(double));
    int count = get_neighbors(points, n, center_idx, eps, &neighbors, &distances);

    // Iterate across all neighbors
    for (int i = 0; i < count; i++)
    {
        int n_idx = neighbors[i];
        Point *n = &points[n_idx];

        // Ignore already processed points
        if (n->processed) continue;

        double new_reach = fmax(core_dist, distances[i]);

        if (n->reach_distance == INF)
        {
            n->reach_distance = new_reach;
            heap_insert(seeds, n_idx);
        }
        else if (new_reach < n->reach_distance)
        {
            heap_decrease_key(seeds, n_idx, new_reach);
        }
    }
    free(neighbors);
    free(distances);
}

void run_optics(Point *points, int n, double eps, int minPts,
                int **ordering_out, double **reach_out, double **core_out)
{
    // Initialize points as unprocessed
    for (int i = 0; i < n; i++)
    {
        points[i].processed = false;
        points[i].core_distance = INF;
        points[i].reach_distance = INF;
    }

    // Initialize variables
    BinaryHeap seeds;
    heap_init(&seeds, points, n);

    int order_idx = 0;
    int *ordering = (int *)malloc(n * sizeof(int));
    double *reach_vals = (double *)malloc(n * sizeof(double));
    double *core_vals = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        // Ignore already processed points
        if (points[i].processed) continue;

        Point *p = &points[i];

        // Calculate core distance and mark as processed
        p->processed = true;
        p->core_distance = compute_core_distance(points, n, i, eps, minPts);

        // Save variables to return arrays
        ordering[order_idx] = i;
        reach_vals[order_idx] = p->reach_distance;
        core_vals[order_idx] = p->core_distance;
        order_idx++;

        // If it is core point
        if (p->core_distance != INF)
        {
            // Populate the heap with p's neighbors
            update_seeds(points, n, i, &seeds, eps, minPts);

            while (!heap_is_empty(&seeds))
            {
                // Extract the closest point and compute distances
                int q_idx = heap_extract_min(&seeds);
                Point *q = &points[q_idx];
                q->processed = true;
                q->core_distance = compute_core_distance(points, n, q_idx, eps, minPts);

                // Save variables to return arrays
                ordering[order_idx] = q_idx;
                reach_vals[order_idx] = q->reach_distance;
                core_vals[order_idx] = q->core_distance;
                order_idx++;

                if (q->core_distance != INF)
                {
                    update_seeds(points, n, q_idx, &seeds, eps, minPts);
                }
            }
        }
    }

    heap_free(&seeds);
    *ordering_out = ordering;
    *reach_out = reach_vals;
    *core_out = core_vals;
}

int main (int argc, char **argv)
{
    // Initialize variables
    if (argc != 5)
    {
        printf("Usage: input.csv epsilon minPts output.csv\n");
        return 1;
    }

    // Search radius for points (eps)
    double eps = atof(argv[2]);
    if (eps <= 0)
    {
        fprintf(stderr, "Invalid epsilon\n");
        return 1;
    }

    // Minimum points in epsilon radius for core point (minPts)
    int minPts = atoi(argv[3]);
    if (minPts <= 1)
    {
        fprintf(stderr, "Invalid minPts\n");
        return 1;
    }

    // Number of points (n) and dimensionality (dim)
    int n, dim;

    // Dynamic Point (struct) array from csv (points)
    Point *points = load_csv(argv[1], &n, &dim);

    // Run OPTICS
    int *ordering;
    double *reach_vals;
    double *core_vals;
    run_optics(points, n, eps, minPts, &ordering, &reach_vals, &core_vals);

    // Open output file
    FILE *out = fopen(argv[4], "w");
    if (!out)
    {
        fprintf(stderr, "Cannot open output %s\n", argv[4]);
        return 1;
    }

    // Write headers
    fprintf(out, "order_idx,original_point_idx,reachability_distance,core_distance\n");

    // Write values
    for (int i = 0; i < n; i++)
    {
        int idx = ordering[i];
        // r = UNDEFINED if reach_vals[i] == INF, else r = reach_vals [i]
        double r = (reach_vals[i] == INF ? UNDEFINED : reach_vals[i]);
        double c = (core_vals[i] == INF ? UNDEFINED : core_vals[i]);
        fprintf(out, "%d,%d,%.10f,%.10f\n", i, points[idx].orig_idx, r, c);
    }

    // Close output file
    fclose(out);

    // Free memory
    free(ordering);
    free(reach_vals);
    free(core_vals);
    free(points[0].coords); // points[0].coords points to coords_block
    free(points);

    printf("OPTICS completed. Ordering and distances written to %s\n", argv[4]);
    return 0;
}
