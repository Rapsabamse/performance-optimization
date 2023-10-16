/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "ppm.hpp"
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>

namespace PPM {

void Reader::fill(std::string filename)
{
    std::ifstream f {};

    f.open(filename);

    if (!f) {
        stream.setstate(std::ios::failbit);
        return;
    }

    stream << f.rdbuf(); // put file data in stream for get_magic, get_dimensions and get_color_max to use

    f.seekg(0, std::istream::end);
    std::streamoff len = f.tellg(); // get length of file
    f.seekg(0);


    std::string line;
    res.clear();
    if (len > 0) {
        res.reserve(static_cast<std::string::size_type>(len)); // reserve space to prevent expensive reallocations
    }
    while (getline(f, line)) {
        (res += line) += "\n"; // read file into string for get_data to read pixel color values from
    }

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

std::tuple<unsigned char*, unsigned char*, unsigned char*> Reader::get_data(unsigned x_size, unsigned y_size)
{
    auto size { x_size * y_size };
    auto R { new char[size] }, G { new char[size] }, B { new char[size] };

    unsigned int j = stream.tellg(); // variable to keep reading the next element, starts at current position in stream


    for (auto i { 0 }, read { 0 }; i < size; i++) { // iterate through the R, G and B arrays

        R[i] = res[j];      // store the red pixel value in R[i]
        G[i] = res[j + 1];  // store the green pixel value in G[i]
        B[i] = res[j + 2];  // store the blue pixel value in B[i]

        j+=3;  // increment index in res by 3 so next value that is red is the next green pixel value
        

        // delete memory if data was not present
        if (&R[i] == nullptr || &G[i] == nullptr || &B[i] == nullptr) {
            delete[] R;
            delete[] G;
            delete[] B;
            return { nullptr, nullptr, nullptr };
        }
    }

    return { reinterpret_cast<unsigned char*>(R), reinterpret_cast<unsigned char*>(G), reinterpret_cast<unsigned char*>(B) };
}

Matrix Reader::operator()(std::string filename)
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

        auto [R, G, B] { get_data(x_size, y_size) };

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
