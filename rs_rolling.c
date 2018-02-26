#include "rs_rolling.h"

rs_rolling *rs_rolling_alloc(rs_vector *v, size_t window) {
        rs_rolling *r = malloc(sizeof(rs_rolling));
        if (!r) {
                fprintf(stderr, "[rs_rolling_alloc] malloc error\n");
                return NULL;
        } else {
                size_t size = v->count - window;
                r->source_data = v;
                r->window = window;
                r->window_data = circular_array_alloc(window);
                r->mins = rs_vector_alloc(size);
                r->maxs = rs_vector_alloc(size);
                r->sums = rs_vector_alloc(size);
                r->means = rs_vector_alloc(size);
                r->variances = rs_vector_alloc(size);
                r->stddevs = rs_vector_alloc(size);
                r->skews = rs_vector_alloc(size);
                r->kurts = rs_vector_alloc(size);

                if (!r->window_data || !r->mins || !r->maxs || !r->sums ||
                    !r->means || !r->variances || !r->stddevs || !r->skews ||
                    !r->kurts) {
                        fprintf(stderr, "[rs_rolling_alloc] malloc error\n");
                        exit(1);
                }
        }
        return r;
}

void rs_rolling_roll(rs_rolling *r, size_t start_index) {
        if (r) {

                for (size_t i = start_index; i < (start_index + r->window); i++) {
                        circular_array_put(r->window_data, rs_vector_get(r->source_data, i));
                }
                rs_vector_item_push(r->sums, circular_array_sum(r->window_data));
                rs_vector_item_push(r->mins, circular_array_min(r->window_data));
                rs_vector_item_push(r->maxs, circular_array_max(r->window_data));
                rs_vector_item_push(r->means, circular_array_mean(r->window_data));
                rs_vector_item_push(r->variances, circular_array_variance(r->window_data));
                rs_vector_item_push(r->stddevs, circular_array_stddev(r->window_data));
                rs_vector_item_push(r->skews, circular_array_skewness(r->window_data));
                rs_vector_item_push(r->kurts, circular_array_kurtosis(r->window_data));
                r->count++;

                // now roll baby, roll
                size_t count = r->source_data->count;
                // fprintf(stderr, "r->v->count: %zu\n", count);

                double tmp = 0.0;

                for (size_t i = (start_index + r->window); i < count; i++) {
                        circular_array_get(r->window_data, &tmp);
                        circular_array_put(r->window_data, rs_vector_get(r->source_data, i));
                        rs_vector_item_push(r->sums, circular_array_sum(r->window_data));
                        rs_vector_item_push(r->mins, circular_array_min(r->window_data));
                        rs_vector_item_push(r->maxs, circular_array_max(r->window_data));
                        rs_vector_item_push(r->means, circular_array_mean(r->window_data));;
                        rs_vector_item_push(r->variances, circular_array_variance(r->window_data));
                        rs_vector_item_push(r->stddevs, circular_array_stddev(r->window_data));
                        rs_vector_item_push(r->skews, circular_array_skewness(r->window_data));
                        rs_vector_item_push(r->kurts, circular_array_kurtosis(r->window_data));
                        r->count++;
                }
                fprintf(stderr,
                        "r->v->count %zu r->window_data->count %zu r->count %zu "
                        "r->sums->count %zu\n",
                        r->source_data->count, r->window_data->size, r->count, r->sums->count);
        }
}

void rs_rolling_free(rs_rolling *r) {
        if (r) {
                if (r->window_data) {
                        circular_array_free(r->window_data);
                }
                if (r->mins) free(r->mins);
                if (r->maxs) free(r->maxs);
                if (r->sums) free(r->sums);
                if (r->means) free(r->means);
                if (r->variances) free(r->variances);
                if (r->stddevs) free(r->stddevs);
                if (r->skews) free(r->skews);
                if (r->kurts) free(r->kurts);
                // if (r->v) rs_vector_free(r->v);
                free(r);
                r = NULL;
        }
}

void rs_rolling_print(rs_rolling *r) {
        if (r) {
                fprintf(stdout, "Sum;Min;Max;Mean;Variance;Stddev;Skew;Kurt\n");

                for (size_t i = 0; i < r->count; i++) {
                        fprintf(stdout,
                                "%.6f;%.6f;%.6f;%.6f;%.6f;%.6f;%.6f;%.6f\n",
                                rs_vector_get(r->sums, i),
                                rs_vector_get(r->mins, i),
                                rs_vector_get(r->maxs, i),
                                rs_vector_get(r->means, i),
                                rs_vector_get(r->variances, i),
                                rs_vector_get(r->stddevs, i),
                                rs_vector_get(r->skews, i),
                                rs_vector_get(r->kurts, i));
                }
        }
}
