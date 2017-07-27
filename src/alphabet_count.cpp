#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <array>

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage: prog FILE\n";
        return 1;
    }

    auto hist = std::array<uint64_t, 256>();

    for(size_t i = 0; i < 256; i++) {
        hist[i] = 0;
    }

    std::ifstream myfile (argv[1]);
    if (myfile.is_open())
    {
        char c;

        while(myfile.get(c)) {
            hist[uint8_t(c)]++;
        }

        myfile.close();
    }

    size_t alpha = 0;
    size_t size = 0;

    for(size_t i = 0; i < 256; i++) {
        alpha += (hist[i] != 0);
        size += hist[i];
    }

    for(size_t i = 0; i < 256; i++) {
        if (hist[i] != 0) {
            std::cout
                << "char ["
                << std::setw(3) << int(i)
                << "]: " << hist[i]
                << ", " << (double(hist[i] * 10000 / size) / 100.0) << "%"
                << "\n";
        }
    }

    std::cout << "size: " << size << "\n";
    std::cout << "alphabet size: " << alpha << "\n";

    return 0;
}
