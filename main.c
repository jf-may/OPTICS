#include <stdio.h>
#include <stdlib.h>

#include "optics.h"
#include "csv_parser.h"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s <input_csv> <output_csv> <epsilon> "
                        "<min_pts>.\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *input_filename  = argv[1];
    char *output_filename = argv[2];
    double epsilon        = atof(argv[3]);
    int min_pts           = atoi(argv[4]);

    printf("\n=== Input CSV file: %s ===\n\n", input_filename);
    CsvData data = load_csv(input_filename);
    printf("Loaded %d points of dimension %d\n", data.n, data.dim);

    ClusterOrdering results = run_optics(data.pts, data.n, epsilon, min_pts);

    printf("\n=== Output CSV file: %s ===\n\n", output_filename);
    save_cluster_ordering_to_csv(&results, data.n, output_filename);

    free_cluster_ordering(&results);
    free_csv_data(&data);

    return EXIT_SUCCESS;
}
