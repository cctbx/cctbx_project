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
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
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
  boost::shared_ptr< ImageTile<T> >
  make_image_tile(typename scitbx::af::flex<T>::type data) {
    DXTBX_ASSERT(data.accessor().all().size() == 2);
    return boost::make_shared< ImageTile<T> >(
        scitbx::af::versa<T, scitbx::af::c_grid<2> >(
          data.handle(),
          scitbx::af::c_grid<2>(data.accessor())));
  }

  template <typename T>
  boost::shared_ptr< ImageTile<T> >
  make_image_tile_with_name(
      typename scitbx::af::flex<T>::type data,
      const char *name) {
    std::cout << "H" << std::endl;
    DXTBX_ASSERT(data.accessor().all().size() == 2);
    return boost::make_shared< ImageTile<T> >(
        scitbx::af::versa<T, scitbx::af::c_grid<2> >(
          data.handle(),
          scitbx::af::c_grid<2>(data.accessor())),
        name);
  }

  template <typename T>
  boost::shared_ptr< Image<T> >
  make_image_from_tuple(boost::python::tuple data) {
    typedef typename scitbx::af::flex<T>::type flex_type;
    boost::shared_ptr< Image<T> > result(new Image<T>());
    for (std::size_t i = 0; i < boost::python::len(data); ++i) {
      flex_type a = boost::python::extract<flex_type>(data[i])();
      DXTBX_ASSERT(a.accessor().all().size() == 2);
      result->push_back(ImageTile<T>(
        scitbx::af::versa<T, scitbx::af::c_grid<2> >(
          a.handle(),
          scitbx::af::c_grid<2>(a.accessor()))));
    }
    return result;
  }

  template <typename T>
  boost::shared_ptr< Image<T> >
  make_image_from_object(boost::python::object data) {

    if (data != boost::python::object()) {
      DXTBX_ERROR("No conversion to Image");
    }
    return boost::make_shared< Image<T> >();
  }

  template <typename T>
  boost::shared_ptr< Image<T> >
  make_image_from_flex(typename scitbx::af::flex<T>::type data) {
    DXTBX_ASSERT(data.accessor().all().size() == 2);
    return boost::make_shared< Image<T> >(ImageTile<T>(
        scitbx::af::versa<T, scitbx::af::c_grid<2> >(
          data.handle(),
          scitbx::af::c_grid<2>(data.accessor()))));
  }

  template <typename T>
  struct ImageTilePickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getinitargs(ImageTile<T> obj) {
      return boost::python::make_tuple(obj.data(), obj.name());
    }

  };

  template <typename T>
  struct ImagePickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getstate(const Image<T> &obj) {
      boost::python::list tile_list;
      for (std::size_t i = 0; i < obj.n_tiles(); ++i) {
        tile_list.append(obj.tile(i));
      }
      return boost::python::make_tuple(tile_list);
    }

    static
    void setstate(Image<T> &obj, boost::python::tuple state) {
      DXTBX_ASSERT(boost::python::len(state) == 1);
      boost::python::list tile_list = boost::python::extract<boost::python::list>(state[0])();
      for (std::size_t i = 0; i < boost::python::len(tile_list); ++i) {
        obj.push_back(boost::python::extract< ImageTile<T> >(tile_list[i])());
      }
    }

  };

  template <typename T>
  void image_tile_wrapper(const char *name) {

    typedef ImageTile<T> image_tile_type;
    typedef typename image_tile_type::array_type array_type;

    class_<image_tile_type, boost::shared_ptr< ImageTile<T> > >(name, no_init)
      .def("__init__", make_constructor(&make_image_tile<T>))
      .def("__init__", make_constructor(&make_image_tile_with_name<T>))
      .def("name", &image_tile_type::name)
      .def("data", &image_tile_type::data)
      .def("empty", &image_tile_type::empty)
      .def_pickle(ImageTilePickleSuite<T>())
      ;
  }

  template <typename T>
  void image_wrapper(const char *name) {

    typedef Image<T> image_type;
    typedef typename image_type::tile_type tile_type;

    class_<image_type>(name)
      .def(init< tile_type >())
      .def("__init__", make_constructor(&make_image_from_flex<T>))
      .def("__init__", make_constructor(&make_image_from_tuple<T>))
      .def("__getitem__", &image_type::tile)
      .def("tile", &image_type::tile)
      .def("tile_names", &image_type::tile_names)
      .def("n_tiles", &image_type::n_tiles)
      .def("empty", &image_type::empty)
      .def("append", &image_type::push_back)
      .def("__len__", &image_type::n_tiles)
      .def("__iter__", range(&image_type::begin, &image_type::end))
      .def_pickle(ImagePickleSuite<T>())
      ;

  }


  BOOST_PYTHON_MODULE(dxtbx_format_image_ext)
  {
    image_tile_wrapper<bool>("ImageTileBool");
    image_tile_wrapper<int>("ImageTileInt");
    image_tile_wrapper<double>("ImageTileDouble");
    image_wrapper<bool>("ImageBool");
    image_wrapper<int>("ImageInt");
    image_wrapper<double>("ImageDouble");

    class_<ImageBuffer>("ImageBuffer")
      .def(init< Image<int> >())
      .def(init< Image<double> >())
      .def("is_empty", &ImageBuffer::is_empty)
      .def("is_int", &ImageBuffer::is_int)
      .def("is_double", &ImageBuffer::is_double)
      .def("as_int", &ImageBuffer::as_int)
      .def("as_double", &ImageBuffer::as_double)
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
