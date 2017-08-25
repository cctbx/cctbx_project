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

  BOOST_PYTHON_MODULE(dxtbx_ext)
  {
    init_module();
    export_to_ewald_sphere_helpers();
  }

}} //namespace dxtbx::boost_python
