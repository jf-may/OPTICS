#define _POSIX_C_SOURCE 200809L

#include "csv_parser.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * Internal helpers
 * ========================================================================= */

/*
 * parse_csv_numeric_row
 *
 * Parse one CSV line into an array of doubles.
 *
 * Parameters
 * ----------
 * line       : null-terminated input line.
 * values_out : destination array; when NULL the fields are counted but not
 *              stored (first-pass mode).
 *
 * Returns
 * -------
 *  > 0 : number of fields parsed (data row).
 *    0 : blank or comment line.
 *   -1 : malformed row (empty field, trailing comma, unexpected character).
 */
static int parse_csv_numeric_row(const char *line, double *values_out)
{
    int count     = 0;
    const char *p = line;

    while (*p != '\0' && isspace((unsigned char) *p))
    {
        p++;
    }

    /* Blank or comment line */
    if (*p == '\0' || *p == '\n' || *p == '\r' || *p == '#')
        return 0;

    /* Leading comma → empty first field */
    if (*p == ',')
        return -1;

    for (;;)
    {
        char  *endptr;
        double val = strtod(p, &endptr);

        if (endptr == p)        /* no numeric token could be parsed */
            return -1;

        if (values_out != NULL)
            values_out[count] = val;
        count++;

        /* Advance past optional whitespace after the number */
        p = endptr;
        while (*p == ' ' || *p == '\t')
            p++;

        if (*p == ',')
        {
            p++;
            while (*p == ' ' || *p == '\t')  /* whitespace after separator */
                p++;
            if (*p == '\n' || *p == '\r' || *p == '\0')
                return -1;                   /* trailing comma */
        }
        else if (*p == '\n' || *p == '\r' || *p == '\0')
        {
            break;                           /* clean end of line */
        }
        else
        {
            return -1;                       /* unexpected character */
        }
    }

    return count;
}

/*
 * csv_fatal
 *
 * Release all held resources, print an error message, and abort.
 */
static void csv_fatal(FILE *f, char *line, double *coords_block,
                       Point *points, const char *msg, const char *filename)
{
    if (f) fclose(f);
    free(line);
    free(coords_block);
    free(points);
    fprintf(stderr, "%s: %s\n", filename, msg);
    exit(EXIT_FAILURE);
}

/* ============================================================================
 * Public API
 * ========================================================================= */

CsvData load_csv(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL)
    {
        fprintf(stderr, "Cannot open %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char   *line      = NULL;
    size_t  line_size = 0;
    int     n         = 0;
    int     dim       = -1;

    /* First pass: count valid data rows and infer dimensionality. */
    while (getline(&line, &line_size, f) != -1)
    {
        int cols    = parse_csv_numeric_row(line, NULL);

        if (cols == 0) continue;    /* blank or comment */

        if (cols < 0)
            csv_fatal(f, line, NULL, NULL, "malformed or empty row", filename);

        if (dim == -1)
        {
            dim = cols;          /* set from first valid row */
        }
        else if (cols != dim)
        {
            csv_fatal(f, line, NULL, NULL, "inconsistent number of columns",
                      filename);
        }

        n++;
    }

    if (n == 0)
        csv_fatal(f, line, NULL, NULL, "no data rows found", filename);

    rewind(f);

    /* Allocate storage: one contiguous coords block + point array.*/
    double *coords_block = (double *)malloc((size_t)n * dim * sizeof(double));
    Point *points        = (Point  *)malloc((size_t)n * sizeof(Point));

    if (!coords_block || !points)
        csv_fatal(f, line, coords_block, points, "out of memory", filename);

    /* Second pass: populate the Point array.*/
    int row = 0;
    while (getline(&line, &line_size, f) != -1)
    {
        double *coords_row = coords_block + (size_t)row * dim;
        int     cols       = parse_csv_numeric_row(line, coords_row);

        if (cols == 0) continue;

        points[row].coords         = coords_row;
        points[row].orig_idx       = row;
        points[row].processed      = false;
        points[row].core_distance  = INFINITY;
        points[row].reach_distance = INFINITY;
        row++;
    }

    free(line);
    fclose(f);

    CsvData result;
    result.pts = points;
    result.n = n;
    result.dim = dim;

    return result;
}

void free_csv_data(CsvData *data)
{
    if (!data) return;

    /* The first point's coords pointer is the start of the contiguous block */
    free(data->pts[0].coords);
    free(data->pts);
}
