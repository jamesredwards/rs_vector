#ifndef __RS_ROLLING_H_
#define __RS_ROLLING_H_

#include "rs_vector.h"

typedef struct rs_rolling {
        size_t window;
        size_t count;
        rs_vector *source_data;
        circular_array *window_data;
        rs_vector *mins;
        rs_vector *maxs;
        rs_vector *sums;
        rs_vector *means;
        rs_vector *variances;
        rs_vector *stddevs;
        rs_vector *skews;
        rs_vector *kurts;
} rs_rolling;

rs_rolling *rs_rolling_alloc(rs_vector *v, size_t window);
void rs_rolling_free(rs_rolling *r);
void rs_rolling_roll(rs_rolling *r, size_t start_index);
void rs_rolling_print(rs_rolling *r);

#endif