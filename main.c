#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    if (argc < 3)
    {
        fprintf(stderr, "Usage ./optics {input filename} {output filename}.\n");
        return EXIT_FAILURE;
    }

    char *input_filename = argv[1];
    printf("\n=== Input CSV file: %s ===\n\n", input_filename);

    CsvData data = load_csv(input_filename);
    printf("Loaded %d points of dimension %d\n\n", data.n, data.dim);

    ClusterOrdering results = run_optics(data.pts, data.n, 5.0, 100);

    char *output_filename = argv[2];
    printf("\n=== Output CSV file: %s ===\n\n", output_filename);

    /* For debugging */
    print_results(results.ordering, results.core, results.reach, data.n);

    save_cluster_ordering_to_csv(&results, data.n, output_filename);

    free_cluster_ordering(&results);
    free_csv_data(&data);

    return EXIT_SUCCESS;
}
