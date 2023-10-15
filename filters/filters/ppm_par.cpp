/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "ppm.hpp"
#include <pthread.h>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>

namespace PPM {

int Reader::get_file_size(std::string filename)
{
    std::ifstream f;
    f.open(filename);

    if (!f) { return -1; }

    f.seekg(0, std::istream::end);
    std::streamoff len = f.tellg();
    f.seekg(0);
    return len;
}

void Reader::fill(std::string filename)
{
    std::ifstream f {};

    f.open(filename);

    if (!f) {
        stream.setstate(std::ios::failbit);
        return;
    }

    stream << f.rdbuf();

    f.close();
}

std::string Reader::get_magic_number()
{
    std::string magic {};

    std::getline(stream, magic);

    return magic;
}

std::pair<unsigned, unsigned> Reader::get_dimensions()
{
    std::string line {};

    while (std::getline(stream, line) && line[0] == '#')
        ;

    std::regex regex { "^(\\d+) (\\d+)$" };
    std::smatch matches {};

    std::regex_match(line, matches, regex);

    if (matches.ready()) {
        return { std::stoul(matches[1]), std::stoul(matches[2]) };
    } else {
        return { 0, 0 };
    }
}

unsigned Reader::get_color_max()
{
    std::string line {};

    std::getline(stream, line);

    std::regex regex { "^(\\d+)$" };
    std::smatch matches {};

    std::regex_match(line, matches, regex);

    if (matches.ready()) {
        return std::stoul(matches[1]);
    } else {
        return 0;
    }
}

struct thread_data {
    std::string filename;
    char* R;
    char* G;
    char* B;
    unsigned int size;
    unsigned int start_pos;
    unsigned int stream_length;
    unsigned int thread_id;
};

void* read_rgb_par(void* thread_args)
{
    thread_data* my_data;
    my_data = (thread_data*) thread_args;
    std::string line;
    std::string res;

    std::ifstream f;
    f.open(my_data->filename);

    res.clear();
    if (my_data->stream_length > 0) {
        res.reserve(static_cast<std::string::size_type>(my_data->stream_length));
    }

    unsigned int start_index = my_data->start_pos + ((my_data->stream_length - (my_data->start_pos / 4)) * my_data->thread_id); // calculate start index
    unsigned int end_index = start_index + my_data->stream_length;                                 // calculate end index


    f.seekg(start_index);                      // set the buffer to the position we want to start reading from
    while (getline(f, line)) {
        if (f.tellg() == end_index) { break; } // read until we reach end_index
        (res += line) += "\n";
    }
    f.close();

    unsigned int j = 0; // to start reading from pos 0 in res

    for (auto i { my_data->size * my_data->thread_id }; i < (my_data->size * (1 + my_data->thread_id)); i++) {  // for-loop from start index
                                                                                                   // to end index for each thread
        my_data->R[i] = res[j];
        my_data->G[i] = res[j + 1];
        my_data->B[i] = res[j + 2];
        j += 3;


        if (&(my_data->R[i]) == nullptr || &(my_data->G[i]) == nullptr || &(my_data->B[i]) == nullptr) {
            delete[] my_data->R;
            delete[] my_data->G;
            delete[] my_data->B;
            // return { nullptr, nullptr, nullptr };
        }
    }
    pthread_exit(NULL);
}

std::tuple<unsigned char*, unsigned char*, unsigned char*> Reader::get_data_par(unsigned x_size, unsigned y_size, std::string filename, int MAX_THREADS)
{
    auto size { x_size * y_size };
    auto R { new char[size] }, G { new char[size] }, B { new char[size] };


    int original_pos = stream.tellg();
    stream.seekg(0, std::istream::end);
    std::streamoff len = stream.tellg(); // get length of stream
    stream.seekg(0);

    thread_data thread_data_array[MAX_THREADS];
    pthread_t p_threads[MAX_THREADS];


    for (int i = 0; i < MAX_THREADS; i++){
            thread_data_array[i].thread_id = i;
            thread_data_array[i].R = R;
            thread_data_array[i].G = G;
            thread_data_array[i].B = B;
            thread_data_array[i].size = size / MAX_THREADS;
            thread_data_array[i].start_pos = original_pos;
            thread_data_array[i].stream_length = len / MAX_THREADS;
            thread_data_array[i].filename = filename;

            pthread_create(&p_threads[i], NULL, read_rgb_par, (void*) &thread_data_array[i]);
        }

        for (int i = 0; i < MAX_THREADS; i++) {
            pthread_join(p_threads[i], NULL);
        }

    return { reinterpret_cast<unsigned char*>(R), reinterpret_cast<unsigned char*>(G), reinterpret_cast<unsigned char*>(B) };
}

Matrix Reader::operator()(std::string filename, int MAX_THREADS)
{
    try {
        fill(filename);

        if (stream.fail()) {
            throw std::runtime_error { "couldn't open file " + filename };
        }

        auto magic { get_magic_number() };

        if (magic != magic_number) {
            throw std::runtime_error { "incorrect magic number: " + magic };
        }

        auto [x_size, y_size] { get_dimensions() };

        if (x_size == 0 || y_size == 0) {
            throw std::runtime_error { "couldn't read dimensions" };
        }

        auto total_size { x_size * y_size };

        if (total_size > max_pixels) {
            throw std::runtime_error { "image size is too big: " + std::to_string(total_size) };
        }

        auto color_max { get_color_max() };

        if (color_max == 0) {
            throw std::runtime_error { "couldn't read color max" };
        }

        auto [R, G, B] { get_data_par(x_size, y_size, filename, MAX_THREADS) };

        if (!R || !G || !B) {
            throw std::runtime_error { "couldn't read image data" };
        }

        stream.clear();
        return Matrix { R, G, B, x_size, y_size, color_max };
    } catch (std::runtime_error e) {
        error("reading", e.what());
        stream.clear();
        return Matrix {};
    }
}

void error(std::string op, std::string what)
{
    std::cerr << "Encountered PPM error during " << op << ": " << what << std::endl;
}

void Writer::operator()(Matrix m, std::string filename)
{
    try {
        std::ofstream f {};

        f.open(filename);

        if (!f) {
            throw std::runtime_error { "failed to open " + filename };
        }

        f << magic_number << std::endl;

        f << m.get_x_size() << " " << m.get_y_size() << std::endl;
        f << m.get_color_max() << std::endl;

        auto size { m.get_x_size() * m.get_y_size() };
        auto R { m.get_R() }, G { m.get_G() }, B { m.get_B() };
        auto it_R { R }, it_G { G }, it_B { B };

        while (it_R < R + size && it_G < G + size && it_B < B + size) {
            f << *it_R++
              << *it_G++
              << *it_B++;
        }

        f.close();
    } catch (std::runtime_error e) {
        error("writing", e.what());
    }
}

}