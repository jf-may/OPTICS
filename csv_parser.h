#ifndef CSV_PARSER_H
#define CSV_PARSER_H

#include "optics_types.h"

typedef struct
{
    Point *pts;
    int n;
    int dim;
} CsvData;

/*
 * load_csv
 * --------
 * Load a numeric CSV file into an array of Point structures.
 *
 * Format:
 *   - One point per line.
 *   - Comma-separated floating-point values.
 *   - Blank lines and lines starting with '#' (after optional spaces) ignored.
 *   - All data rows must have the same number of columns.
 *
 * Parameters
 * ----------
 * filename : const char *    Input file path.
 *
 * Returns
 * -------
 * Dynamically allocated Point array.  Must be released with free_points().
 */
CsvData load_csv(const char *filename);

/*
 * free_csv_data
 * -----------
 * Releases memory returned by load_csv() (coords block + Point array).
 */
void free_csv_data(CsvData *data);

#endif
