/*
 * nexus_ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/slice.hpp>
#include <scitbx/array_family/flex_types.h>
#include <dxtbx/error.h>
#include <vector>
#include <hdf5.h>

namespace dxtbx { namespace format { namespace boost_python {

  using namespace boost::python;


  /**
   * A function to extract data from hdf5 dataset directly into a flex array
   */
  inline
  scitbx::af::versa<int, scitbx::af::flex_grid<> > dataset_as_flex_int(
      hid_t dataset_id,
      boost::python::tuple selection) {

    // The number of dimensions
    std::size_t ndims = boost::python::len(selection);

    // Get the file space
    hid_t file_space_id = H5Dget_space(dataset_id);
    std::size_t rank = H5Sget_simple_extent_ndims(file_space_id);
    DXTBX_ASSERT(rank == ndims);
    std::vector<hsize_t> dataset_dims(rank);
    H5Sget_simple_extent_dims(file_space_id, &dataset_dims[0], NULL);

    // Create the grid
    scitbx::af::flex_grid<>::index_type dims(ndims);
    std::vector<hsize_t> start(ndims);
    std::vector<hsize_t> count(ndims);
    for (std::size_t i = 0; i < ndims; ++i) {
      boost::python::slice slice = boost::python::extract<boost::python::slice>(selection[i]);
      int slice_start = boost::python::extract<int>(slice.start());
      int slice_stop  = boost::python::extract<int>(slice.stop());
      int slice_step  = boost::python::extract<int>(slice.step());
      DXTBX_ASSERT(slice_step == 1);
      DXTBX_ASSERT(slice_stop > slice_start);
      dims[i] = slice_stop - slice_start;
      start[i] = slice_start;
      count[i] = dims[i];
      DXTBX_ASSERT(start[i] + count[i] <= dataset_dims[i]);
    }

    // Create the data array
    scitbx::af::flex_grid<> grid(dims);
    scitbx::af::versa<int, scitbx::af::flex_grid<> > data(grid, scitbx::af::init_functor_null<int>());

    // Create the dataspace id
    herr_t status1 = H5Sselect_hyperslab(
        file_space_id,
        H5S_SELECT_SET,
        &start[0],
        NULL,
        &count[0],
        NULL);
    DXTBX_ASSERT(status1 >= 0);

    // Create the memory space size
    hid_t mem_space_id = H5Screate_simple(ndims, &count[0], NULL);

    // Copy the data
    herr_t status2 = H5Dread(
        dataset_id,
        H5T_NATIVE_INT,
        mem_space_id,
        file_space_id,
        H5P_DEFAULT,
        &data[0]);
    DXTBX_ASSERT(status2 >= 0);

    // Close some stuff
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);

    // Return the data
    return data;
  }

  BOOST_PYTHON_MODULE(dxtbx_format_nexus_ext)
  {
    def("dataset_as_flex_int", &dataset_as_flex_int);
  }

}}} // namespace = dxtbx::format::boost_python
