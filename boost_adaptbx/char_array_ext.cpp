#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

namespace {

  template <std::size_t N>
  struct char_n_holder
  {
    typedef char char_n[N];
    char_n value;

    char_n_holder()
    {
      const char* v = "barte";
      for(std::size_t i=0;i<N;i++) {
        value[i] = v[i];
      }
    }
  };

  template <unsigned N>
  std::string
  use_char_n(char const(&s)[N])
  {
    std::string result;
    for(unsigned i=N;i>0;) {
      i--;
      result += std::string(&s[i], 1);
    }
    return result;
  }

  template <std::size_t N>
  struct char_n_to_python_str
  {
    static PyObject*
    convert(char const(&s)[N])
    {
      return boost::python::incref(
        boost::python::object(
          std::string(s, N)).ptr());
    }
  };

  template <std::size_t N>
  struct char_n_from_python_str
  {
    typedef char char_n[N];

    char_n_from_python_str()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<char_n>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
      if (!PyString_Check(obj_ptr)) return 0;
      if (PyString_GET_SIZE(obj_ptr) != N) return 0;
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      char* py_str = PyString_AsString(obj_ptr);
      if (py_str == 0) boost::python::throw_error_already_set();
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<char_n>*)
          data)->storage.bytes;
      new (storage) char_n();
      data->convertible = storage;
      char_n& result = *((char_n*)storage);
      for(std::size_t i=0;i<N;i++) {
        result[i] = py_str[i];
      }
    }
  };

  void init_module()
  {
    using namespace boost::python;

    typedef char char_3[3];
    to_python_converter<char_3, char_n_to_python_str<3> >();
    char_n_from_python_str<3>();
    class_<char_n_holder<3> >("char_3_holder")
      .def_readonly("value", &char_n_holder<3>::value)
    ;
    def("use_char_n", (std::string(*)(char_3 const&)) use_char_n);

    typedef char char_5[5];
    to_python_converter<char_5, char_n_to_python_str<5> >();
    char_n_from_python_str<5>();
    class_<char_n_holder<5> >("char_5_holder")
      .def_readonly("value", &char_n_holder<5>::value)
    ;
    def("use_char_n", (std::string(*)(char_5 const&)) use_char_n);
  }

} // namespace <anonymous>

BOOST_PYTHON_MODULE(boost_adaptbx_char_array_ext)
{
  init_module();
}
