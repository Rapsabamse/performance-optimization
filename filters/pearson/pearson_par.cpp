/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis_par.hpp"
#include "dataset.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char const* argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [num_threads]" << std::endl;
        std::exit(1);
    }

    auto datasets { Dataset::read(argv[1]) };
    auto corrs { Analysis::correlation_coefficients(datasets, std::stoul(argv[3])) };
    Dataset::write(corrs, argv[2]);

    return 0;
}
