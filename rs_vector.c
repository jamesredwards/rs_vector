#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include "rs_vector.h"

rs_vector *rs_vector_alloc(size_t init_capacity)
{
    if (init_capacity == 0) {
        init_capacity = 1;
    }
    rs_vector *v = malloc(sizeof(rs_vector));
    if (!v) {
        fprintf(stderr, "[rs_vector_alloc] malloc error\n");
        return NULL;
    }
    v->capacity = init_capacity;
    v->data = malloc(sizeof(double) * v->capacity);
    return v;
}

void rs_vector_free(rs_vector *v)
{
    if (v) {
        if (v->data) {
            free(v->data);
            v->data = NULL;
        }
        free(v);
        v = NULL;
    }
}

void rs_vector_reset(rs_vector *v)
{
    //zero everything out
    memset(v->data, 0, v->count * sizeof(double));
    v->mean = 0.0;
    v->M2 = 0.0;
    v->M3 = 0.0;
    v->M4 = 0.0;
    v->sum = 0.0;
    v->min = 0.0;
    v->max = 0.0;
}

int rs_vector_resize(rs_vector *v, size_t new_size)
{
    v->capacity = new_size;
    double *data = realloc(v->data, v->capacity * sizeof(double));
    if (data) {
        v->data = data;
        return 0;
    }
    return -1;
}

int rs_vector_item_push(rs_vector *v, double item)
{
    if (v->count == v->capacity - 1) {
        rs_vector_expand(v);
    }
    rs_vector_update(v, item);
    v->data[v->count] = item;
    return 0;
}

double rs_vector_item_pop(rs_vector *v)
{
    double item = v->data[v->count];
    v->data[v->count] = 0.0;
    //v->count--; will be (de)incremented in update
    rs_vector_update_remove(v, item);
    if (v->count == (v->capacity / CAPACITY_INCREASE_FACTOR) - 1) {
        rs_vector_contract(v);
    }
    return item;
}

double rs_vector_item_pop_index(rs_vector *v, size_t index)
{
    double item = v->data[index];
    v->data[index] = 0.0;
    //v->count--; will be (de)incremented in update
    rs_vector_update_remove(v, item);
    
    //TODO work out logic for contracting if popped from
    //start of array
    /*if (v->count == (v->capacity / CAPACITY_INCREASE_FACTOR) - 1)
    {
        rs_vector_contract(v);
    }*/
    return item;
}

int rs_vector_expand(rs_vector *v)
{
    size_t new_size = v->capacity * CAPACITY_INCREASE_FACTOR;
    int rc = rs_vector_resize(v, new_size);
    return rc;
}

int rs_vector_contract(rs_vector *v)
{
    size_t new_size = v->capacity / CAPACITY_INCREASE_FACTOR;
    int rc = rs_vector_resize(v, new_size);
    return rc;
}

void rs_vector_print_stats(rs_vector *v)
{
    printf("-------------------\nStats\n-------------------\n");
    printf("Min      : %7.6f\n", rs_vector_min(v));
    printf("Max      : %7.6f\n", rs_vector_max(v));
    printf("Sum      : %7.6f\n", rs_vector_sum(v));
    printf("Mean     : %7.6f\n", rs_vector_mean(v));
    printf("Variance : %7.6f\n", rs_vector_variance(v));
    printf("Std dev  : %7.6f\n", rs_vector_stddev(v));
    printf("Skewness : %7.6f\n", rs_vector_skewness(v));
    printf("Kurtosis : %7.6f\n", rs_vector_kurtosis(v));
}

void rs_vector_update(rs_vector *v, double item)
{
    v->sum += item;
    double delta = item - v->mean;
    double delta_n, delta_nsq, term1, n;

    /* update min and max */
    if (v->count == 0)
    {
        v->min = item;
        v->max = item;
    }
    else
    {
        if (item < v->min)
            v->min = item;
        if (item > v->max)
            v->max = item;
    }

    /* from GSL rstat and John D. Cook, MIT license
        http://www.johndcook.com/blog/skewness_kurtosis/ */
    n = (double)++(v->count);
    delta_n = delta / n;
    delta_nsq = delta_n * delta_n;
    term1 = delta * delta_n * (n - 1.0);
    v->mean += delta_n;
    v->M4 += term1 * delta_nsq * (n * n - 3.0 * n + 3.0) +
             6.0 * delta_nsq * v->M2 - 4.0 * delta_n * v->M3;
    v->M3 += term1 * delta_n * (n - 2.0) - 3.0 * delta_n * v->M2;
    v->M2 += term1;
}

/* need to use n+1 after decrement in term1, M3/M4 */
/* TODO algo to update min/max                     */
void rs_vector_update_remove(rs_vector *v, double item)
{
    v->sum -= item;
    double delta, delta_n, delta_nsq, term1;

    double n1 = (double)v->count;
    double n = (double)--v->count;
    delta = item - v->mean;
    delta_n = delta / n;
    delta_nsq = delta_n * delta_n;
    term1 = delta * delta_n * n1;
    v->mean -= delta_n;
    v->M2 -= term1;
    v->M3 -= term1 * delta_n * (n1 - 2.0) - 3.0 * delta_n * v->M2;
    v->M4 -= term1 * delta_nsq * (n1 * n1 - 3.0 * n1 + 3.0) +
             6.0 * delta_nsq * v->M2 - 4.0 * delta_n * v->M3;
}

rs_rolling *rs_rolling_alloc(rs_vector *v, size_t window)
{
    rs_rolling *r = malloc(sizeof(rs_rolling));
    if (!r)
    {
        fprintf(stderr, "[rs_rolling_alloc] malloc error\n");
        return NULL;
    }
    else 
    {
        r->v = v;
        r->window = window;
        r->data = rs_vector_alloc(1);
        r->mins = rs_vector_alloc(1);
        r->maxs = rs_vector_alloc(1);
        r->sums = rs_vector_alloc(1);
        r->means = rs_vector_alloc(1);
        r->variances = rs_vector_alloc(1);
        r->stddevs = rs_vector_alloc(1);
        r->skews = rs_vector_alloc(1);
        r->kurts = rs_vector_alloc(1);

        if (!r->data || !r->mins || !r->maxs || !r->sums || !r->means ||
            !r->variances || !r->stddevs || !r->skews || !r->kurts)
            {
                fprintf(stderr, "[rs_rolling_alloc] malloc error\n");
                exit(1);
            }
    }
    return r;
}
void rs_rolling_roll(rs_rolling *r, size_t start_index)
{
    if (r)
    {
        //initial calc to push data onto r->data
        for (size_t i = start_index; i < r->window; i++)
        {
            rs_vector_item_push(r->data, rs_vector_get(r->v, i));
        }
        rs_vector_item_push(r->sums, rs_vector_sum(r->v));
        rs_vector_item_push(r->mins, rs_vector_min(r->v));
        rs_vector_item_push(r->maxs, rs_vector_max(r->v));
        rs_vector_item_push(r->means, rs_vector_mean(r->v));
        rs_vector_item_push(r->variances, rs_vector_variance(r->v));
        rs_vector_item_push(r->stddevs, rs_vector_stddev(r->v));
        rs_vector_item_push(r->skews, rs_vector_skewness(r->v));
        rs_vector_item_push(r->kurts, rs_vector_kurtosis(r->v));

        //now roll baby, roll
        for (size_t i = r->window; i < r->v->count; i++)
        {
            rs_vector_item_pop_index(r->data, i - r->window);
            rs_vector_item_push(r->data, rs_vector_get(r->v, i));
            rs_vector_item_push(r->sums, rs_vector_sum(r->data));
            rs_vector_item_push(r->mins, rs_vector_min(r->data));
            rs_vector_item_push(r->maxs, rs_vector_max(r->data));
            rs_vector_item_push(r->means, rs_vector_mean(r->data));
            rs_vector_item_push(r->variances, rs_vector_variance(r->data));
            rs_vector_item_push(r->stddevs, rs_vector_stddev(r->data));
            rs_vector_item_push(r->skews, rs_vector_skewness(r->data));
            rs_vector_item_push(r->kurts, rs_vector_kurtosis(r->data));
        }
    }
}

void rs_rolling_free(rs_rolling *r)
{
    if (r)
    {
        if (r->data)
        {
            free(r->data);
            r->data = NULL;
        }
        if (r->mins) free(r->mins);
        if (r->maxs) free(r->maxs);
        if (r->sums) free(r->sums);
        if (r->means) free(r->means);
        if (r->variances) free(r->variances);
        if (r->stddevs) free(r->stddevs);
        if (r->skews) free(r->skews);
        if (r->kurts) free(r->kurts);
        free(r);
        r = NULL;
    }
}

int main(int argc, char *argv[])
{
    if (argc == 2) {
        //Test numeric stability with gsl random uniform
        rs_vector *v = rs_vector_alloc(1);

        const gsl_rng_type *T;
        gsl_rng *r;

        size_t n_vars = strtoul(argv[1], NULL, 0);

        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);

        for (size_t i = 0; i < n_vars; i++)
        {
            double u = gsl_rng_uniform(r);
            rs_vector_item_push(v, u);
        }

        rs_vector_print_stats(v);

        gsl_rng_free(r);

        rs_vector_free(v);
    }
    else {
        fprintf(stderr, "Please supply a number between 1 and 1,000,000,000\n");
        exit(1);
    }
}