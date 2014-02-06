
#pragma once

#include <vector>

#include <complex>

#include <scitbx/constants.h>
#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/vec3.h>
#include <boost/python.hpp>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/vec3.h>
#include <boost/python.hpp>

/*
struct Disorder
{
    std::vector<std::vector<int> > disorder_group;
    std:;vector<std::vector<int> > alternative_groups;
};
*/


namespace cctbx {
namespace multipolar{

using namespace scitbx;

void multipolar_test();

af::shared<std::string> assign_atom_types( af::shared<int> atomic_numbers,
                           af::shared<scitbx::vec3<double> > coordinates );

 /* class Aspherical{
  public:
    Aspherical();




    void assign_atom_types_1(const std::vector<size_t> &atomic_numbers,
                           const std::vector<double> &xyz,
                           const std::vector<std::vector<int> > connectivity,
                           const Disorder &disorder,
                           std::vector<std::string> &bank_atom_types,
                           std::vector<std::string> &local_coordinates);

    void set_atom_types(std::vector<std::string> &atom_types,
                        std::vector<std::string> &local_coordinates);
    void calculate_f(std::vector<double> &xyz,
                     std::vector<double> &occupancies,
                     std::vector<double> &adps,
                     std::vector<int> &miller_indices,
                     std::vector<std::complex<double> > &sf);
  }; //class Aspherical*/
  } //namespace multipolar
} //namespace cctbx
