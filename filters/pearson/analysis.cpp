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
    int* syncVar;
    unsigned int number_of_threads;
};

void* correlation_coefficients_par(void* thread_args)
{
    thread_data* my_data;
    my_data = (thread_data*) thread_args;

    std::vector<double> parResults;

    unsigned int size = my_data->datasets->size() / my_data->number_of_threads;
    unsigned int start_index = my_data->thread_id * size;
    unsigned int end_index = start_index + size;
    //The last index wasnt read by original function, need to skip it for the last thread
    if(my_data->thread_id == (my_data->number_of_threads - 1)) end_index--;

    int b = 0;
    for (int sample1 { start_index }; sample1 < end_index; sample1 ++) {
        for (int sample2 = sample1 + 1; sample2 < my_data->datasets->size(); sample2++) {
            double corr = pearson((*my_data->datasets)[sample1], (*my_data->datasets)[sample2]);
            //my_data->result->at(result_i) = corr;
            parResults.push_back(corr);
            b++;
        }
    }

    int waiting = true;
    while(waiting){
        if(*my_data->syncVar == my_data->thread_id){
            for(auto i = 0; i < b; i++){
                my_data->result->push_back(parResults[i]);
            }
            *my_data->syncVar++;
            waiting = false;
        }
    }

    pthread_exit(NULL);
}

std::vector<double> correlation_coefficients(std::vector<Vector> datasets, int MAX_THREADS)
{
    std::vector<double> result {};

    unsigned int result_index = 0;

    thread_data thread_data_array[MAX_THREADS];
    pthread_t p_threads[MAX_THREADS];

    pthread_mutex_init(&lock, NULL);

    int syncVar = 0;
    for (int i = 0; i < MAX_THREADS; i++) {
        thread_data_array[i].thread_id = i;
        thread_data_array[i].datasets = &datasets;
        thread_data_array[i].result = &result;
        thread_data_array[i].syncVar = &syncVar;
        thread_data_array[i].number_of_threads = MAX_THREADS;

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