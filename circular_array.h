#ifndef __CIRCULAR_ARRAY_H_
#define __CIRCULAR_ARRAY_H_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
taken from https://embeddedartistry.com/blog/2017/4/6/circular-buffers-in-cc
*/

typedef struct circular_array {
        double *data;
        size_t head;
        size_t tail;
        size_t size;
        double min;
        double max;
        double sum;
        double mean;
        double M2;
        double M3;
        double M4;
} circular_array;

circular_array *circular_array_alloc(size_t size);
void circular_array_free(circular_array *ca);
int circular_array_reset(circular_array *ca);
int circular_array_put(circular_array *ca, double item);
int circular_array_get(circular_array *ca, double *out_value);
bool circular_array_is_empty(circular_array *ca);
bool circular_array_is_full(circular_array *ca);
void circular_array_update_stats_put(circular_array *ca, double item);
void circular_array_update_stats_pop(circular_array *ca, double item);

inline double circular_array_min(circular_array *ca) { return ca->min; }
inline double circular_array_max(circular_array *ca) { return ca->max; }
inline double circular_array_sum(circular_array *ca) { return ca->sum; }
inline double circular_array_mean(circular_array *ca) { return ca->mean; }

inline double circular_array_variance(circular_array *ca) {
        return (ca->M2 / (ca->size - 1.0));
}
inline double circular_array_stddev(circular_array *ca) {
        return sqrt(circular_array_variance(ca));
}
inline double circular_array_skewness(circular_array *ca) {
        double fac = pow(ca->size - 1.0, 1.5) / (ca->size + 0.0);
        return ((fac * ca->M3) / pow(ca->M2, 1.5));
}
inline double circular_array_kurtosis(circular_array *ca) {
        double fac = ((ca->size - 1.0) / ca->size) * (ca->size - 1.0);
        return ((fac * ca->M4) / (ca->M2 * ca->M2) - 3.0);
}

#endif