#include "circular_array.h"

circular_array *circular_array_alloc(size_t size) {
        if (size < 1) {
                fprintf(stderr, "[circular_array_alloc] size must be > 0\n");
                exit(-1);
        }
        circular_array *ca = malloc(sizeof(circular_array));
        if (!ca) {
                fprintf(stderr, "[circular_array_alloc] malloc error\n");
                return NULL;
        }
        ca->size = size;
        ca->data = malloc(sizeof(double) * ca->size);
        circular_array_reset(ca);
        return ca;
}

void circular_array_free(circular_array *ca) {
        if (ca) {
                if (ca->data) {
                        free(ca->data);
                }
                free(ca);
                ca = NULL;
        }
}

int circular_array_reset(circular_array *ca) {
        int rc = -1;

        if (ca) {
                ca->head = 0;
                ca->tail = 0;
                memset(ca->data, 0.0, ca->size * sizeof(double));
                rc = 0;
        }
        return rc;
}

int circular_array_put(circular_array *ca, double item) {
        int rc = -1;

        if (ca) {
                ca->data[ca->head] = item;
                ca->head = (ca->head + 1) % ca->size;
                circular_array_update_stats_put(ca, item);

                if (ca->head == ca->tail) {
                        ca->tail = (ca->tail + 1) % ca->size;
                }
                rc = 0;
        }
        return rc;
}

int circular_array_get(circular_array *ca, double *out_value) {
        int rc = -1;

        if (ca && ca->data && !circular_array_is_empty(ca)) {
                *out_value = ca->data[ca->tail];
                circular_array_update_stats_pop(ca, *out_value);
                rc = 0;
        }
        return rc;
}

bool circular_array_is_empty(circular_array *ca) {
        return (ca->head == ca->tail);
}

bool circular_array_is_full(circular_array *ca) {
        return ((ca->head + 1) % ca->size) == ca->tail;
}

void circular_array_update_stats_put(circular_array *ca, double item)
{
        ca->sum += item;
        double delta = item - ca->mean;
        double delta_n, delta_nsq, term1, n;

        /* update min and max */
        if (ca->size == 0) {
                ca->min = item;
                ca->max = item;
        } else {
                if (item < ca->min)
                        ca->min = item;
                if (item > ca->max)
                        ca->max = item;
        }

        /* from GSL rstat and John D. Cook, MIT license
            http://www.johndcook.com/blog/skewness_kurtosis/ */
        n = (double)ca->size;
        delta_n = delta / n;
        delta_nsq = delta_n * delta_n;
        term1 = delta * delta_n * (n - 1.0);
        ca->mean += delta_n;
        ca->M4 += term1 * delta_nsq * (n * n - 3.0 * n + 3.0) +
                 6.0 * delta_nsq * ca->M2 - 4.0 * delta_n * ca->M3;
        ca->M3 += term1 * delta_n * (n - 2.0) - 3.0 * delta_n * ca->M2;
        ca->M2 += term1;
}

void circular_array_update_stats_pop(circular_array *ca, double item)
{
        ca->sum -= item;
        double delta, delta_n, delta_nsq, term1;

        double n1 = (double)ca->size;
        double n = n1 - 1.0;
        delta = item - ca->mean;
        delta_n = delta / n;
        delta_nsq = delta_n * delta_n;
        term1 = delta * delta_n * n1;
        ca->mean -= delta_n;
        ca->M2 -= term1;
        ca->M3 -= term1 * delta_n * (n1 - 2.0) - 3.0 * delta_n * ca->M2;
        ca->M4 -= term1 * delta_nsq * (n1 * n1 - 3.0 * n1 + 3.0) +
                 6.0 * delta_nsq * ca->M2 - 4.0 * delta_n * ca->M3;
}

void circular_array_print(circular_array *ca)
{
        fprintf(stdout, "circular_array [%.2f, ", ca->data[0]);
        for (size_t i = 1; i < ca->size - 1; i++) {
                fprintf(stdout, "%.2f, ", ca->data[i]);
        }
        fprintf(stdout, "%.2f]\n", ca->data[ca->size-1]);
}