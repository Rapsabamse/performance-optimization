/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <iostream>
#include <cmath>
#include <pthread.h>

Vector::Vector()
    : size { 0 }
    , data { nullptr }
{
}

Vector::~Vector()
{
    if (data) {
        delete[] data;
    }

    size = 0;
}

Vector::Vector(unsigned size)
    : size { size }
    , data { new double[size] }
{
}

Vector::Vector(unsigned size, double* data)
    : size { size }
    , data { data }
{
}

struct thread_vector{
        int thread_id;
        int thread_amount;
        int size;
        double* data;
        double* dataOther;
};

void *threadVector(void * thread_arg){
    //Remake threadarg from void ptr to thread_data_blur struct so thread can use values
    struct thread_vector *my_data;
    my_data = (struct thread_vector *) thread_arg;

    for (auto i { my_data->thread_id }; i < my_data->size; my_data->thread_amount) {
        data[i] = other.data[i];
    }
    pthread_exit(NULL);
}

Vector::Vector(const Vector& other, int MAX_THREADS)
    : Vector { other.size }
{
    struct thread_vector thread_data_array[MAX_THREADS];
    pthread_t p_threads[MAX_THREADS];

    for (auto i { 0 }; i < size; i++) {
        thread_data_array[i].thread_id = i;
        thread_data_array[i].thread_amount = MAX_THREADS;
        thread_data_array[i].size = size;
        thread_data_array[i].data = data;
        thread_data_array[i].dataOther = other.data;

        //create threads and run threadSblurX, thread_data_array is passed as a parameter
        pthread_create(
            &p_threads[i],
            NULL,
            threadVector,
            (void*) &thread_data_array[i]
        );
    }
}

unsigned Vector::get_size() const
{
    return size;
}

double* Vector::get_data()
{
    return data;
}

double Vector::operator[](unsigned i) const
{
    return data[i];
}

double& Vector::operator[](unsigned i)
{
    return data[i];
}

double Vector::mean() const
{
    double sum { 0 };

    for (auto i { 0 }; i < size; i++) {
        sum += data[i];
    }

    return sum / static_cast<double>(size);
}

double Vector::magnitude() const
{
    auto dot_prod { dot(*this) };
    return std::sqrt(dot_prod);
}

Vector Vector::operator/(double div)
{
    auto result { *this };

    for (auto i { 0 }; i < size; i++) {
        result[i] /= div;
    }

    return result;
}

Vector Vector::operator-(double sub)
{
    auto result { *this };

    for (auto i { 0 }; i < size; i++) {
        result[i] -= sub;
    }

    return result;
}

double Vector::dot(Vector rhs) const
{
    double result { 0 };

    for (auto i { 0 }; i < size; i++) {
        result += data[i] * rhs[i];
    }

    return result;
}