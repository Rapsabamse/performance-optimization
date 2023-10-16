/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"
#include "ppm.hpp"
#include "filters.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char const* argv[])
{
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " [radius] [infile] [outfile] [num_threads]" << std::endl;
        std::exit(1);
    }

    PPM::Reader reader {};
    PPM::Writer writer {};

    auto m { reader(argv[2]) };
    auto radius { static_cast<unsigned>(std::stoul(argv[1])) };

    Filter::blur_par(m, radius, std::stoul(argv[4]));
    writer(m, argv[3]);

    return 0;
}