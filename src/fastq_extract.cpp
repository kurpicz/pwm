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

    std::ifstream myfile (argv[1]);

        char c;
    size_t pos = 0;

    auto check = [&](auto c_should) {
        if (c != c_should) {
            std::cerr << "File err @" << pos << ": " << c << "should be " << c_should << "\n";
            std::exit(1);
        }
    };

    if (myfile.is_open())
    {

        while(myfile) {
            while(myfile.get(c)) {
                pos++;
                if (c == '\n') {
                    break;
                }
            }

            while(myfile.get(c)) {
                pos++;
                if (c == '\n') {
                    break;
                }
                std::cout << c;
            }

            check('+');

            while(myfile.get(c)) {
                pos++;
                if (c == '\n') {
                    break;
                }
            }

            while(myfile.get(c)) {
                pos++;
                if (c == '\n') {
                    break;
                }
            }
        }

        myfile.close();
    }


    return 0;
}
