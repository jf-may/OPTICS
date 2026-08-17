#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "optics.h"
#include "csv_parser.h"

static void print_results(const int *ordering, const double *core,
                          const double *reach, const int n)
{
    printf("%-3s %-6s %-9s %-10s\n",
           "pos", "pt_idx", "core_dist", "reach_dist");
    printf("--- ------ --------- ----------\n");

    for (int i = 0; i < n; i++)
    {
        printf("%-3d %-6d ", i, ordering[i]);

        if (isinf(core[i]))   printf("%-9s ", "INF");
        else                  printf("%-9.4f ", core[i]);

        if (isinf(reach[i]))  printf("%-10s\n", "INF");
        else                  printf("%-10.4f\n", reach[i]);
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage ./optics {filename}.\n");
        return EXIT_FAILURE;
    }

    char *filename = argv[1];
    printf("\n=== CSV file: %s ===\n\n", filename);

    CsvData data = load_csv(filename);
    printf("Loaded %d points of dimension %d\n\n", data.n, data.dim);

    ClusterOrdering res = run_optics(data.pts, data.n, 1.0, 2);
    print_results(res.ordering, res.core, res.reach, data.n);

    free_cluster_ordering(&res);
    free_csv_data(&data);

    return EXIT_SUCCESS;
}
