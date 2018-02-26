#include "rs_vector.h"
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <string.h>
#include "circular_array.h"
#include "rs_rolling.h"

rs_vector *rs_vector_alloc(size_t init_capacity) {
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
        memset(v->data, 0, v->count * sizeof(double));
        v->mean = 0.0;
        v->M2 = 0.0;
        v->M3 = 0.0;
        v->M4 = 0.0;
        v->sum = 0.0;
        v->min = 0.0;
        v->max = 0.0;
}

int rs_vector_resize(rs_vector *v, size_t new_size) {
        // fprintf(stderr, "v->count %zu, v->capacity %zu, new_size %zu\n",
        // v->count, v->capacity, new_size);
        v->capacity = new_size;
        double *data = realloc(v->data, v->capacity * sizeof(double));
        if (data) {
                v->data = data;
                return 0;
        }
        return -1;
}

int rs_vector_item_push(rs_vector *v, double item) {
        if (v->count == v->capacity - 1) {
                rs_vector_expand(v);
        }
        rs_vector_update(v, item);
        v->data[v->count] = item;
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

double rs_vector_item_pop_index(rs_vector *v, size_t index) {
        double item = v->data[index];
        v->data[index] = 0.0;
        // v->count--; will be (de)incremented in update
        rs_vector_update_remove(v, item);

        // TODO work out logic for contracting if popped from
        // start of array
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
        /*if (argc == 2) {
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
            fprintf(stderr, "Please supply a number between 1 and
        1,000,000,000\n"); exit(1);
        }*/

        if (argc > 1) {
                size_t n_vars = strtoul(argv[1], NULL, 0);
                size_t window = strtoul(argv[2], NULL, 0);
                rs_vector *v = rs_vector_alloc(n_vars);

                for (size_t i = 5; i < n_vars + 5; i++) {
                        rs_vector_item_push(v, (double)i);
                }

                // fprintf(stderr, "rs_vector set up correctly\n");

                rs_rolling *r = rs_rolling_alloc(v, window);
                rs_rolling_roll(r, 0);
                rs_rolling_print(r);

                rs_rolling_free(r);
                rs_vector_free(v);
                // printf("%zu\n", sizeof(double));
                // printf("%zu\n", sizeof(rs_vector));
        }
}