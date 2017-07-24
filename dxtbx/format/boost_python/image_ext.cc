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
#include <dxtbx/format/tiff.h>
#include <dxtbx/format/cbf.h>
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


  BOOST_PYTHON_MODULE(dxtbx_format_image_ext)
  {

    class_<ImageTile>("ImageTile", no_init)
      .def("as_int", &ImageTile::as_int)
      .def("as_double", &ImageTile::as_double)
      .def("is_int", &ImageTile::is_int)
      .def("is_double", &ImageTile::is_double)
      .def("name", &ImageTile::name)
      ;

    class_<Image>("Image", no_init)
      .def("tile", &Image::tile)
      .def("tile_names", &Image::tile_names)
      .def("n_tiles", &Image::n_tiles)
      ;

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
