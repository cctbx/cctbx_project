/*
 * image_ext.cc
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
#include <dxtbx/format/image.h>
#include <dxtbx/format/smv.h>
#include <vector>
#include <hdf5.h>

namespace dxtbx { namespace format { namespace boost_python {

  using namespace boost::python;



  BOOST_PYTHON_MODULE(dxtbx_format_image_ext)
  {
    class_<ImageReader>("ImageReader", no_init)
      .def("filename", &ImageReader::filename)
      .def("size", &ImageReader::size)
      .def("type", &ImageReader::type)
      .def("byte_order", &ImageReader::byte_order)
      .def("is_int", &ImageReader::is_int)
      .def("is_float", &ImageReader::is_float)
      .def("as_int", &ImageReader::as_int)
      .def("as_double", &ImageReader::as_double)
      ;

    class_<SMVReader, bases<ImageReader> >("SMVReader", no_init)
      .def(init<const char*>((
              arg("filename"))))
      ;
  }

}}} // namespace = dxtbx::format::boost_python
