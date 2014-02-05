
#pragma once

#include <vector>

#include <complex>

namespace cctbx {
  class Aspherical{
  public:
    Aspherical();

    void assign_atom_types(std::vector<size_t> &atomic_numbers,
                           std::vector<double> &xyz,
                           std::vector<std::string> &bank_atom_types,
                           std::vector<std::string> &local_coordinates);
    void set_atom_types(std::vector<std::string> &atom_types,
                        std::vector<std::string> &local_coordinates);
    void calculate_f(std::vector<double> &xyz,
                     std::vector<double> &occupancies,
                     std::vector<double> &adps,
                     std::vector<int> &miller_indices,
                     std::vector<std::complex<double> > &sf);
  }; //class Aspherical
} //namespace cctbx

