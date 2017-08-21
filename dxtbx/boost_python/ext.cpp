#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <boost_adaptbx/python_streambuf.h>
#include <boost/cstdint.hpp>
#include <cctype>
#include <fstream>
#include <vector>
#include <limits>
#include <dxtbx/imageset.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace boost_python {

  using boost::uint32_t;

  bool
  is_big_endian()
  {
    uint32_t val = 1;
    char * buff = (char *) & val;

    return (buff[0] == 0);
  }

  scitbx::af::shared<int>
  read_uint8(boost_adaptbx::python::streambuf & input,
             size_t count)
  {
    scitbx::af::shared<int> result;
    boost_adaptbx::python::streambuf::istream is(input);
    std::vector<unsigned char> data;
    data.resize(count);

    is.read((char *) &data[0], count * sizeof(unsigned char));

    for (size_t j = 0; j < count; j++) {
      result.push_back((int) data[j]);
    }

    return result;
  }

  scitbx::af::shared<int>
  read_uint16(boost_adaptbx::python::streambuf & input,
              size_t count)
  {
    scitbx::af::shared<int> result;
    boost_adaptbx::python::streambuf::istream is(input);
    std::vector<unsigned short> data;
    data.resize(count);

    is.read((char *) &data[0], count * sizeof(unsigned short));

    for (size_t j = 0; j < count; j++) {
      result.push_back((int) data[j]);
    }

    return result;
  }

  scitbx::af::shared<int>
  read_uint32(boost_adaptbx::python::streambuf & input,
              size_t count)
  {
    scitbx::af::shared<int> result;
    boost_adaptbx::python::streambuf::istream is(input);
    std::vector<unsigned int> data;
    data.resize(count);

    is.read((char *) &data[0], count * sizeof(unsigned int));

    for (size_t j = 0; j < count; j++) {
      DXTBX_ASSERT (data[j] <= std::numeric_limits<int>::max());
      result.push_back((int) data[j]);
    }

    return result;
  }

  scitbx::af::shared<int>
  read_uint16_bs(boost_adaptbx::python::streambuf & input,
                 size_t count)
  {
    scitbx::af::shared<int> result;
    boost_adaptbx::python::streambuf::istream is(input);
    std::vector<unsigned short> data;
    data.resize(count);

    is.read((char *) &data[0], count * sizeof(unsigned short));

    /* swap bytes */

    for (size_t j = 0; j < count; j++) {
      unsigned short x = data[j];
      data[j] = x >> 8 | x << 8;
    }

    for (size_t j = 0; j < count; j++) {
      result.push_back((int) data[j]);
    }

    return result;
  }

  scitbx::af::shared<int>
  read_uint32_bs(boost_adaptbx::python::streambuf & input,
                 size_t count)
  {
    scitbx::af::shared<int> result;
    boost_adaptbx::python::streambuf::istream is(input);
    std::vector<unsigned int> data;
    data.resize(count);

    is.read((char *) &data[0], count * sizeof(unsigned int));

    /* swap bytes */

    for (size_t j = 0; j < count; j++) {
      unsigned int x = data[j];
      data[j] = (x << 24) | (x << 8 & 0xff0000) |
        (x >> 8 & 0xff00) | (x >> 24);
    }

    for (size_t j = 0; j < count; j++) {
      DXTBX_ASSERT (data[j] <= std::numeric_limits<int>::max());
      result.push_back((int) data[j]);
    }

    return result;
  }

  scitbx::af::shared<int>
  read_int16(boost_adaptbx::python::streambuf & input,
             size_t count)
  {
    scitbx::af::shared<int> result;
    boost_adaptbx::python::streambuf::istream is(input);
    std::vector<short> data;
    data.resize(count);

    is.read((char *) &data[0], count * sizeof(short));

    for (size_t j = 0; j < count; j++) {
      result.push_back((int) data[j]);
    }

    return result;
  }

  scitbx::af::shared<int>
  read_int32(boost_adaptbx::python::streambuf & input,
             size_t count)
  {
    scitbx::af::shared<int> result;
    boost_adaptbx::python::streambuf::istream is(input);
    std::vector<int> data;
    data.resize(count);

    is.read((char *) &data[0], count * sizeof(int));

    for (size_t j = 0; j < count; j++) {
      result.push_back((int) data[j]);
    }

    return result;
  }

  scitbx::af::shared<double>
  read_float32(boost_adaptbx::python::streambuf & input,
             size_t count)
  {
    scitbx::af::shared<double> result;
    boost_adaptbx::python::streambuf::istream is(input);
    std::vector<float> data;
    data.resize(count);

    is.read((char *) &data[0], count * sizeof(float));

    for (size_t j = 0; j < count; j++) {
      result.push_back((double) data[j]);
    }

    return result;
  }

  void init_module()
  {
    using namespace boost::python;
    def("read_uint8", read_uint8, (arg("file"), arg("count")));
    def("read_uint16", read_uint16, (arg("file"), arg("count")));
    def("read_uint32", read_uint32, (arg("file"), arg("count")));
    def("read_uint16_bs", read_uint16_bs, (arg("file"), arg("count")));
    def("read_uint32_bs", read_uint32_bs, (arg("file"), arg("count")));
    def("read_int16", read_int16, (arg("file"), arg("count")));
    def("read_int32", read_int32, (arg("file"), arg("count")));
    def("read_float32", read_float32, (arg("file"), arg("count")));
    def("is_big_endian", is_big_endian);
  }


  void export_to_ewald_sphere_helpers();

  void export_imageset() {
    using namespace boost::python;


    class_< ExternalLookupItem<double> >("ExternalLookupItemDouble")
      .add_property("filename",
          &ExternalLookupItem<double>::get_filename,
          &ExternalLookupItem<double>::set_filename)
      .add_property("data",
          &ExternalLookupItem<double>::get_data,
          &ExternalLookupItem<double>::set_data)
      ;

    class_< ExternalLookupItem<bool> >("ExternalLookupItemBool")
      .add_property("filename",
          &ExternalLookupItem<bool>::get_filename,
          &ExternalLookupItem<bool>::set_filename)
      .add_property("data",
          &ExternalLookupItem<bool>::get_data,
          &ExternalLookupItem<bool>::set_data)
      ;

    class_<ExternalLookup>("ExternalLookup")
      .add_property("mask",
          make_function(
            &ExternalLookup::mask,
            return_internal_reference<>()))
      .add_property("gain",
          make_function(
            &ExternalLookup::gain,
            return_internal_reference<>()))
      .add_property("pedestal",
          make_function(
            &ExternalLookup::pedestal,
            return_internal_reference<>()))
      ;

    class_<ImageSetData>("ImageSetDataWrapper", no_init)
      .def(init<boost::python::object>())
      .def("get_data", &ImageSetData::get_data)
      .def("get_mask", &ImageSetData::get_mask)
      .def("has_single_file_reader", &ImageSetData::has_single_file_reader)
      .def("get_path", &ImageSetData::get_path)
      .def("get_master_path", &ImageSetData::get_master_path)
      .def("get_image_identifier", &ImageSetData::get_image_identifier)
      .def("get_property", &ImageSetData::get_property)
      .def("set_property", &ImageSetData::set_property)
      .def("get_beam", &ImageSetData::get_beam)
      .def("get_detector", &ImageSetData::get_detector)
      .def("get_goniometer", &ImageSetData::get_goniometer)
      .def("get_scan", &ImageSetData::get_scan)
      .def("set_beam", &ImageSetData::set_beam)
      .def("set_detector", &ImageSetData::set_detector)
      .def("set_goniometer", &ImageSetData::set_goniometer)
      .def("set_scan", &ImageSetData::set_scan)
      .add_property("external_lookup",
          make_function(
            &ImageSet::external_lookup,
            return_internal_reference<>()))
      ;

    class_<ImageSet>("ImageSet", no_init)
      .def(init<
          const ImageSetData &
          >())
      .def(init<
          const ImageSetData &,
          const scitbx::af::const_ref<std::size_t> &
          >())
      .def("data", &ImageSet::data)
      .def("indices", &ImageSet::indices)
      .def("size", &ImageSet::size)
      .def("get_raw_data", &ImageSet::get_raw_data)
      .def("get_corrected_data", &ImageSet::get_corrected_data)
      .def("get_gain", &ImageSet::get_gain)
      .def("get_pedestal", &ImageSet::get_pedestal)
      .def("get_mask", &ImageSet::get_mask)
      .def("get_property", &ImageSet::get_property)
      .def("set_property", &ImageSet::set_property)
      .def("get_beam", &ImageSet::get_beam)
      .def("get_detector", &ImageSet::get_detector)
      .def("get_goniometer", &ImageSet::get_goniometer)
      .def("get_scan", &ImageSet::get_scan)
      .def("set_beam", &ImageSet::set_beam)
      .def("set_detector", &ImageSet::set_detector)
      .def("set_goniometer", &ImageSet::set_goniometer)
      .def("set_scan", &ImageSet::set_scan)
      .def("get_path", &ImageSet::get_path)
      .def("get_image_identifier", &ImageSet::get_image_identifier)
      .def("as_imageset", &ImageSet::as_imageset)
      .def("complete_set", &ImageSet::complete_set)
      .def("partial_set", &ImageSet::partial_set)
      .def("__eq__", &ImageSet::operator==)
      .def("__ne__", &ImageSet::operator!=)
      .add_property("external_lookup",
          make_function(
            &ImageSet::external_lookup,
            return_internal_reference<>()))
      ;

    class_<ImageGrid, bases<ImageSet> >("ImageGrid", no_init)
      .def(init<
          const ImageSetData &,
          int2
          >())
      .def(init<
          const ImageSetData &,
          const scitbx::af::const_ref<std::size_t> &,
          int2
          >())
      .def("get_grid_size", &ImageGrid::get_grid_size)
      .def("from_imageset", &ImageGrid::from_imageset)
      ;

    class_<ImageSweep, bases<ImageSet> >("ImageSweep", no_init)
      .def(init<
          const ImageSetData &,
          const Beam &,
          const Detector &,
          const Goniometer &,
          const Scan &
          >())
      .def(init<
          const ImageSetData &,
          const scitbx::af::const_ref<std::size_t> &,
          const Beam &,
          const Detector &,
          const Goniometer &,
          const Scan &
          >())
      .def("get_array_range", &ImageSweep::get_array_range)
      .def("complete_set", &ImageSweep::complete_sweep)
      .def("partial_set", &ImageSweep::partial_sweep)
      ;
  }


  BOOST_PYTHON_MODULE(dxtbx_ext)
  {
    init_module();
    export_to_ewald_sphere_helpers();
    export_imageset();
  }

}} //namespace dxtbx::boost_python
