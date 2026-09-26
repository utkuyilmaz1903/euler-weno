#include "euler/version.hpp"

#include <cstdlib>
#include <iostream>

int main() {
    std::cout << "euler-weno " << euler::version() << '\n';
    return EXIT_SUCCESS;
}
