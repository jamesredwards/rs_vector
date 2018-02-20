#ifndef __RS_VECTOR_H_
#define __RS_VECTOR_H_

#include <gsl/gsl_rstat.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define CAPACITY_INCREASE_FACTOR 2

typedef struct rs_vector {
    double min;
    double max;
    double sum;
    double mean;;
    double M2;
    double M3;
    double M4;
    size_t count;
    size_t capacity;
    double *data;
} rs_vector;

rs_vector *rs_vector_alloc(size_t init_capacity);
void rs_vector_free(rs_vector *v);
rs_vector *rs_vector_reset(rs_vector *v);
int rs_vector_resize(rs_vector *v, size_t new_size);
int rs_vector_expand(rs_vector *v);
int rs_vector_contract(rs_vector *v);
int rs_vector_item_push(rs_vector *v, double item);
double rs_vector_item_pop(rs_vector *v);
void rs_vector_update(rs_vector *v, double item);
void rs_vector_update_remove(rs_vector *v, double item);

inline double rs_vector_min(rs_vector *v)
{
    return v->min;
}
inline double rs_vector_max(rs_vector *v)
{
    return v->max;
}
inline double rs_vector_sum(rs_vector *v)
{
    return v->sum;
}
inline double rs_vector_mean(rs_vector *v)
{
    return v->mean;
}
inline double rs_vector_variance(rs_vector *v)
{
    return (v->M2 / (v->count - 1.0));
}
inline double rs_vector_stddev(rs_vector *v)
{
    return sqrt(rs_vector_variance(v));
}
inline double rs_vector_skewness(rs_vector *v)
{
    double fac = pow(v->count - 1.0, 1.5) / (v->count + 0.0);
    return ((fac * v->M3) / pow(v->M2, 1.5));
}
inline double rs_vector_kurtosis(rs_vector *v)
{
    double fac = ((v->count - 1.0) / v->count) * (v->count - 1.0);
    return ((fac * v->M4) / (v->M2 * v->M2) - 3.0);
}

void rs_vector_print_stats(rs_vector *v);

static inline double rs_vector_get(rs_vector *v, size_t index)
{
    if (index > v->count - 1) {
        fprintf(stderr, "[rs_vector_get] index %zu out of range %zu", index, v->count);
        return 0.0;
    }
    return v->data[index];
}


#endif