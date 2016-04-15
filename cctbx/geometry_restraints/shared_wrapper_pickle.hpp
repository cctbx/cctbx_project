
#ifndef CCTBX_GEOMETRY_RESTRAINTS_SHARED_WRAPPER_PICKLE
#define CCTBX_GEOMETRY_RESTRAINTS_SHARED_WRAPPER_PICKLE


#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/class.hpp>

/*
A workaround for enabling pickling in some classes used by the restraints manager:
pair_sym_table
shared_dihedral_proxy
shared_parallelity_proxy
shared_planarity_proxy
shared_chirality_proxy
shared_angle_proxy
bond_params_table

These all are templated from scitbx::af::boost_python::shared_wrapper
Alternatively pickling could have been implemented in
scitbx/array_family/boost_python/shared_wrapper.h.
But this would conflict with the pickling functions in
scitbx/array_family/shared.py which would need to be
commented out which in turn might potentially break backwards
compatibility when loading old projects.

*/

template<typename Array>
struct shared_wrapper_pickle_suite : boost::python::pickle_suite
{
  static
    boost::python::tuple
    getinitargs(Array const& stat)
  {
    return boost::python::make_tuple(boost::python::list(stat));
  }
};


#endif
