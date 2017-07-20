#include "algorithm.hpp"

int main(int argc, char const *argv[]) {
  auto& algo_list = algorithm_list::get_algorithm_list();
  for (const auto& a : algo_list) {
    a->print_info();
  }
  // std::cout << "HIER" << std::endl;
  return 0;
}