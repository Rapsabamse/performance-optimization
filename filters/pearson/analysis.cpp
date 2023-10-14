/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <pthread.h>

namespace Analysis {

pthread_mutex_t lock;
struct thread_data {
    unsigned int thread_id;
    std::vector<Vector>* datasets;
    std::vector<double>* result;
    unsigned int number_of_threads;
    unsigned int* result_index;
};

void* correlation_coefficients_par(void* thread_args)
{
    thread_data* my_data;
    my_data = (thread_data*) thread_args;

    std::vector<double> parResults;

    unsigned int size = my_data->datasets->size() / my_data->number_of_threads;
    unsigned int start_index = my_data->thread_id * size;
    unsigned int end_index = start_index + size;

    for (int sample1 = my_data->thread_id; sample1 < end_index - 1; sample1 += my_data->number_of_threads) {
        for (int sample2 = sample1 + 1; sample2 < end_index; sample2++) {
            double corr = pearson((*my_data->datasets)[sample1], (*my_data->datasets)[sample2]);
            parResults.insert(std::begin(*my_data->result) + (*my_data->result_index), corr);
            (*my_data->result_index)++;
        }
    }

    pthread_mutex_lock(&lock); // prevent race conditions between threads writing to result
    std::cout << "finished :)\n";
    pthread_mutex_unlock(&lock);

    pthread_exit(NULL);
}

std::vector<double> correlation_coefficients(std::vector<Vector> datasets, int MAX_THREADS)
{
    std::vector<double> result {};
    unsigned int result_index = 0;

    thread_data thread_data_array[MAX_THREADS];
    pthread_t p_threads[MAX_THREADS];

    pthread_mutex_init(&lock, NULL);

    for (int i = 0; i < MAX_THREADS; i++) {
        thread_data_array[i].thread_id = i;
        thread_data_array[i].datasets = &datasets;
        thread_data_array[i].result = &result;
        thread_data_array[i].number_of_threads = MAX_THREADS;
        thread_data_array[i].result_index = &result_index;

        pthread_create(&p_threads[i], NULL, correlation_coefficients_par, (void*) &thread_data_array[i]);
    }

    for (int i = 0; i < MAX_THREADS; i++) {
        pthread_join(p_threads[i], NULL);
    }

    return result;
}

double pearson(Vector vec1, Vector vec2)
{
    auto x_mean { vec1.mean() };
    auto y_mean { vec2.mean() };

    auto x_mm { vec1 - x_mean };
    auto y_mm { vec2 - y_mean };

    auto x_mag { x_mm.magnitude() };
    auto y_mag { y_mm.magnitude() };

    auto x_mm_over_x_mag { x_mm / x_mag };
    auto y_mm_over_y_mag { y_mm / y_mag };

    auto r { x_mm_over_x_mag.dot(y_mm_over_y_mag) };

    return std::max(std::min(r, 1.0), -1.0);
}
};