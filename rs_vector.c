#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "circular_array.h"
#include "rs_rolling.h"
#include "rs_vector.h"

rs_vector *rs_vector_alloc(size_t init_capacity) {
        if (init_capacity == 0) {
                init_capacity = 1;
        }
        rs_vector *v = malloc(sizeof(rs_vector));
        if (!v) {
                fprintf(stderr, "[rs_vector_alloc] malloc error\n");
                return NULL;
        }
        v->capacity = init_capacity + 1;
        v->data = malloc(sizeof(double) * v->capacity);
        rs_vector_reset(v);
        return v;
}

void rs_vector_free(rs_vector *v) {
        if (v) {
                if (v->data) {
                        free(v->data);
                        v->data = NULL;
                }
                free(v);
                v = NULL;
        }
}

void rs_vector_reset(rs_vector *v) {
        // zero everything out
        memset(v->data, 0, v->capacity * sizeof(double));
        v->count = 0;
        v->mean = 0.0;
        v->M2 = 0.0;
        v->M3 = 0.0;
        v->M4 = 0.0;
        v->sum = 0.0;
        v->min = 0.0;
        v->max = 0.0;
}

int rs_vector_resize(rs_vector *v, size_t new_size) {

        static size_t call_ctr = 0;
        call_ctr++;
        
        v->capacity = new_size;
        double *data = realloc(v->data, v->capacity * sizeof(double));
        if (data) {
                v->data = data;
                return 0;
        }
        return -1;
}

int rs_vector_item_push(rs_vector *v, double item) {

        /*fprintf(stderr, 
                "[rs_vector_item_push] count: %zu, current cap: %zu\n", 
                v->count, v->capacity);*/
        static size_t func_call = 0;

        if (v->count == v->capacity - 1) {
                func_call++;
                rs_vector_expand(v);
        }
        v->data[v->count] = item;
        rs_vector_update(v, item);
        return 0;
}

double rs_vector_item_pop(rs_vector *v) {
        double item = v->data[v->count];
        v->data[v->count] = 0.0;
        // v->count--; will be (de)incremented in update
        rs_vector_update_remove(v, item);
        if (v->count == (v->capacity / CAPACITY_INCREASE_FACTOR) - 1) {
                rs_vector_contract(v);
        }
        return item;
}

int rs_vector_expand(rs_vector *v) {

        size_t new_size = v->capacity * CAPACITY_INCREASE_FACTOR;
        int rc = rs_vector_resize(v, new_size);

        return rc;
}

int rs_vector_contract(rs_vector *v) {

        size_t new_size = v->capacity / CAPACITY_INCREASE_FACTOR;
        int rc = rs_vector_resize(v, new_size);
        return rc;
}

void rs_vector_print_stats(rs_vector *v) {
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

void rs_vector_print(rs_vector *v)
{
        fprintf(stdout, "rs_vector [%.2f, ", rs_vector_get(v, 0));
        for (size_t i = 1; i < v->count - 1; i++) {
                fprintf(stdout, "%.2f, ", rs_vector_get(v, i));
        }
        fprintf(stdout, "%.2f]\n", rs_vector_get(v, v->count - 1));
}

void rs_vector_update(rs_vector *v, double item) {
        v->sum += item;
        double delta = item - v->mean;
        double delta_n, delta_nsq, term1, n;

        /* update min and max */
        if (v->count == 0) {
                v->min = item;
                v->max = item;
        } else {
                if (item < v->min) v->min = item;
                if (item > v->max) v->max = item;
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
void rs_vector_update_remove(rs_vector *v, double item) {
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

int main(int argc, char *argv[]) {

        if (argc > 1) {
                size_t n_vars = strtoul(argv[1], NULL, 0);
                size_t window = strtoul(argv[2], NULL, 0);

                
                rs_vector *v = rs_vector_alloc(1);

                for (size_t i = 0; i < n_vars; i++) {
                        rs_vector_item_push(v, (double)i+1);
                }
                //rs_vector_print(v);
                rs_rolling *r = rs_rolling_alloc(v, window);
                rs_rolling_roll(r, 0);
                rs_rolling_print(r);

                rs_rolling_free(r);
                rs_vector_free(v);
        }
}