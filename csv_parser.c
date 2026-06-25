#define _POSIX_C_SOURCE 200809L

#include "csv_parser.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * CSV parser for OPTICS input
 * ============================================================================
 *
 * Supported format:
 *   - One point per line.
 *   - Comma-separated floating-point values (parsed with strtod).
 *   - Blank lines ignored.
 *   - Lines whose first non-space character is '#' are treated as comments.
 *   - Every data row must have the same number of columns.
 *
 * Not supported:
 *   - Quoted fields.
 *   - Mixed types (all values must be numeric).
 *
 * Memory layout:
 *   A single contiguous coords block is allocated for all n*dim values.
 *   Each Point's coords pointer is an offset into that block.
 *   free_points() releases both allocations.
 * ========================================================================= */

/* -------------------------------------------------------------------------
 * Internal helpers
 * ---------------------------------------------------------------------- */

/*
 * skip_spaces
 *
 * Return a pointer to the first non-whitespace character in s.
 */
static char *skip_spaces(char *s)
{
    while (*s != '\0' && isspace((unsigned char)*s))
        s++;
    return s;
}

/*
 * parse_csv_numeric_row
 *
 * Parse one CSV line into an array of doubles.
 *
 * Parameters
 * ----------
 * line   : null-terminated input line.
 * values : destination array; when NULL the fields are counted but not
 *          stored (first-pass mode).
 * out_dim: set to the number of fields parsed on success.
 *
 * Returns
 * -------
 *  > 0 : number of fields parsed (data row).
 *    0 : blank or comment line — caller should skip.
 *   -1 : malformed row (empty field, trailing comma, unexpected character).
 */
static int parse_csv_numeric_row(char *line, double *values, int *out_dim)
{
    int   count = 0;
    char *p     = skip_spaces(line);

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

        if (values != NULL)
            values[count] = val;
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

    if (out_dim != NULL)
        *out_dim = count;

    return count;
}

/*
 * csv_fatal
 *
 * Release all held resources, print an error message, and abort.
 * Passing NULL for any pointer is safe (free(NULL) is a no-op).
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

/* -------------------------------------------------------------------------
 * Public API
 * ---------------------------------------------------------------------- */

/*
 * load_csv
 *
 * Two-pass loader: first pass counts rows and validates dimensionality
 * consistency; second pass populates the Point array.
 */
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

    /* ------------------------------------------------------------------
     * First pass: count valid data rows and infer dimensionality.
     * ---------------------------------------------------------------- */
    while (getline(&line, &line_size, f) != -1)
    {
        int row_dim = 0;
        int cols    = parse_csv_numeric_row(line, NULL, &row_dim);

        if (cols == 0) continue;    /* blank or comment */

        if (cols < 0)
            csv_fatal(f, line, NULL, NULL, "malformed row", filename);

        if (row_dim <= 0)
            csv_fatal(f, line, NULL, NULL, "empty data row", filename);

        if (dim == -1)
        {
            dim = row_dim;          /* set from first valid row */
        }
        else if (row_dim != dim)
        {
            csv_fatal(f, line, NULL, NULL,
                      "inconsistent number of columns", filename);
        }

        n++;
    }

    if (n == 0)
        csv_fatal(f, line, NULL, NULL, "no data rows found", filename);

    rewind(f);

    /* ------------------------------------------------------------------
     * Allocate storage: one contiguous coords block + point array.
     * ---------------------------------------------------------------- */
    double *coords_block = (double *)malloc((size_t)n * (size_t)dim *
                                             sizeof(double));
    Point  *points       = (Point  *)malloc((size_t)n * sizeof(Point));

    if (!coords_block || !points)
        csv_fatal(f, line, coords_block, points, "out of memory", filename);

    /* ------------------------------------------------------------------
     * Second pass: populate the Point array.
     * ---------------------------------------------------------------- */
    int row = 0;
    while (getline(&line, &line_size, f) != -1)
    {
        double *coords_row = coords_block + (size_t)row * (size_t)dim;
        int     row_dim    = 0;
        int     cols       = parse_csv_numeric_row(line, coords_row, &row_dim);

        if (cols == 0) continue;

        points[row].coords         = coords_row;
        points[row].dim            = dim;
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

    data->pts = NULL;
}
