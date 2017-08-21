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
#include <dxtbx/format/image_reader.h>
#include <dxtbx/format/smv_reader.h>
#include <dxtbx/format/tiff_reader.h>
#include <dxtbx/format/cbf_reader.h>
#include <dxtbx/format/hdf5_reader.h>
#include <vector>
#include <hdf5.h>

namespace dxtbx { namespace format { namespace boost_python {

  using namespace boost::python;


  template <typename ImageReaderType>
  void image_list_reader_suite(const char *name) {

    typedef ImageListReader<ImageReaderType> image_list_reader;

    class_<image_list_reader>(name, no_init)
      .def(
          init<const scitbx::af::const_ref<std::string>&>((
            arg("filenames"))))
      .def("filenames",
          &image_list_reader::filenames)
      .def("image",
          &image_list_reader::image, (
            arg("index")))
      .def("__len__",
          &image_list_reader::size)
      ;
  }

  template <typename T>
  void image_tile_wrapper(const char *name) {

    typedef ImageTile<T> image_tile_type;

    class_<image_tile_type>(name, no_init)
      .def("name", &image_tile_type::name)
      .def("data", &image_tile_type::data)
      .def("accessor", &image_tile_type::accessor)
      .def("empty", &image_tile_type::empty)
      ;
  }

  template <typename T>
  void image_wrapper(const char *name) {

    typedef Image<T> image_type;

    class_<image_type>(name, no_init)
      .def("tile", &image_type::tile)
      .def("tile_names", &image_type::tile_names)
      .def("n_tiles", &image_type::n_tiles)
      .def("push_back", &image_type::push_back)
      ;

  }

  BOOST_PYTHON_MODULE(dxtbx_format_image_ext)
  {
    image_tile_wrapper<int>("ImageTileInt");
    image_tile_wrapper<double>("ImageTileDouble");
    image_wrapper<int>("ImageInt");
    image_wrapper<double>("ImageDouble");

    class_<ImageReader>("ImageReader", no_init)
      .def("filename", &ImageReader::filename)
      .def("image", &ImageReader::image)
      ;

    class_<SMVReader, bases<ImageReader> >("SMVReader", no_init)
      .def(init<const char*>((
              arg("filename"))))
      .def("byte_order", &SMVReader::byte_order)
      ;

    class_<TIFFReader, bases<ImageReader> >("TIFFReader", no_init)
      .def(init<const char*>((
              arg("filename"))))
      ;

    class_<CBFFastReader, bases<ImageReader> >("CBFFastReader", no_init)
      .def(init<const char*>((
              arg("filename"))))
      ;

    class_<CBFReader, bases<ImageReader> >("CBFReader", no_init)
      .def(init<const char*>((
              arg("filename"))))
      ;

    class_<HDF5Reader>("HDF5Reader", no_init)
      .def(init<hid_t,
                const scitbx::af::const_ref<std::string>&>((
                    arg("handle"),
                    arg("datasets"))))
      .def(init<hid_t,
                const scitbx::af::const_ref<std::string>&,
                std::size_t,
                std::size_t>((
                    arg("handle"),
                    arg("datasets"),
                    arg("first"),
                    arg("last"))))
      .def("filename", &HDF5Reader::filename)
      .def("handle", &HDF5Reader::handle)
      .def("datasets", &HDF5Reader::datasets)
      .def("first", &HDF5Reader::first)
      .def("last", &HDF5Reader::last)
      .def("image", &HDF5Reader::image)
      .def("__len__", &HDF5Reader::size)
      ;

    image_list_reader_suite<SMVReader>("SMVImageListReader");
    image_list_reader_suite<TIFFReader>("TIFFImageListReader");
    image_list_reader_suite<CBFFastReader>("CBFFastImageListReader");
    image_list_reader_suite<CBFReader>("CBFImageListReader");

  }

}}} // namespace = dxtbx::format::boost_python
