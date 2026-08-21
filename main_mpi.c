#include "optics.h"
#include "csv_parser.h"
#include "doptics.h"
#include "optics_types.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (argc < 5)
    {
        if (rank == 0)
        {
            fprintf(stderr, "Usage: %s <input_csv> <output_csv> <epsilon> "
                            "<min_pts>\n", argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    char *input_filename  = argv[1];
    char *output_filename = argv[2];
    double epsilon        = atof(argv[3]);
    int min_pts           = atoi(argv[4]);

    CsvData data = {0};
    int global_size = 0;
    int dim = 0;

    /* Master node reads data */
    if (rank == 0)
    {
        printf("\n=== Parallel MN_DOPTICS ===\n");
        CsvData data = load_csv(input_filename);
        printf("Master loaded %d points of dimension %d\n", data.n, data.dim);
        global_size = data.n;
        dim = data.dim;
    }

    /* Broadcast dimensions */
    MPI_Bcast(&global_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Data partitioning and distribution */
    Point *local_points = NULL;
    int local_size = 0;
    doptics_scatter_data(data.pts, global_size, dim, rank, num_procs,
                         &local_points, &local_size);

    /* Local OPTICS execution */
    ClusterOrdering local_co = run_optics(local_points, local_size, dim,
                                          epsilon, min_pts);

    /* Hierarchical merge */
    for (int i = 1; i < num_procs; i *= 2)
    {
        /* Receiver */
        if (rank % (2 * i) == 0)
        {
            int sender = rank + i;

            /* If the sender doesn't exist, skip to the next merge level */
            if (sender >= num_procs) continue;

            /* Receive remote points */
            int remote_size;
            MPI_Recv(&remote_size, 1, MPI_INT, sender, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            double *remote_coords = malloc(remote_size * dim * sizeof(double));
            MPI_Recv(remote_coords, remote_size * dim, MPI_DOUBLE, sender, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            /* Reconstruct remote points */
            Point *remote_points = malloc(remote_size * sizeof(Point));
            for (int i = 0; i < remote_size; i++)
            {
                remote_points[i].coords = remote_coords + (i * dim);
                /* Not needed anymore, set to safe values */
                remote_points[i].processed = true;
                remote_points[i].core_distance = INFINITY;
                remote_points[i].reach_distance = INFINITY;
            }

            /* Recieve remote cluster ordering */
            ClusterOrdering remote_co;
            remote_co.ordering = malloc(remote_size * sizeof(int));
            remote_co.core     = malloc(remote_size * sizeof(double));
            remote_co.reach    = malloc(remote_size * sizeof(double));

            MPI_Recv(remote_co.ordering, remote_size, MPI_INT, sender, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(remote_co.core, remote_size, MPI_DOUBLE, sender, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(remote_co.reach, remote_size, MPI_DOUBLE, sender, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            /* Merge cluster orderings */
            ClusterOrdering merged_co = doptics_merge(local_points, local_size,
                                                      &local_co, remote_points,
                                                      remote_size, &remote_co,
                                                      epsilon, min_pts, dim);

            /* Merge point arrays */
            int merged_size = local_size + remote_size;
            Point *merged_points  = malloc(merged_size * sizeof(Point));
            double *merged_coords = malloc(merged_size * sizeof(double));

            for (int i = 0; i < local_size; i++)
            {
                for (int d = 0; d < dim; d++)
                {
                    merged_coords[i * dim + d] = local_points[i].coords[d];
                    merged_points[i] = local_points[i];
                    merged_points[i].coords = merged_coords + (i * dim);
                }
            }

            for (int i = 0; i < remote_size; i++)
            {
                for (int d = 0; d < dim; d++)
                {
                    merged_coords[(local_size + i) * dim + d] = remote_points[i].coords[d];
                    merged_points[local_size + i] = local_points[i];
                    merged_points[local_size + i].coords = merged_coords + ((local_size + i) * dim);
                }
            }

            /* Cleanup */
            free(local_points[0].coords);
            free(local_points);
            free_cluster_ordering(&local_co);
            free(remote_coords);
            free(remote_points[0].coords);
            free(remote_points);
            free_cluster_ordering(&remote_co);

            /* Promote merged structures */
            local_co = merged_co;
            local_points = merged_points;
            local_size = merged_size;
        }
        /* Sender */
        else if (rank % (2 * i) == i)
        {
            int receiver = rank - i;
            MPI_Send(&local_size, 1, MPI_INT, receiver, 0, MPI_COMM_WORLD);

            MPI_Send(local_points[0].coords, local_size * dim, MPI_DOUBLE,
                     receiver, 0, MPI_COMM_WORLD);

            MPI_Send(local_co.ordering, local_size, MPI_INT, receiver, 0,
                     MPI_COMM_WORLD);
            MPI_Send(local_co.core, local_size, MPI_DOUBLE, receiver, 0,
                     MPI_COMM_WORLD);
            MPI_Send(local_co.reach, local_size, MPI_DOUBLE, receiver, 0,
                     MPI_COMM_WORLD);

            break;
        }
    }

    if (rank == 0)
    {
        printf("Saving merged parallel cluster ordering to: %s\n",
               output_filename);
        save_cluster_ordering_to_csv(&local_co, local_size, output_filename);
    }

    free(local_points[0].coords);
    free(local_points);
    free_cluster_ordering(&local_co);
    if (rank == 0) free_csv_data(&data);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
