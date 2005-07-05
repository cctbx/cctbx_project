#ifndef BOOST_ADAPTBX_SWIG_ARG_HPP
# define BOOST_ADAPTBX_SWIG_ARG_HPP

# include <boost/python/detail/prefix.hpp>

# include <boost/python/type_id.hpp>
# include <boost/python/converter/registry.hpp>

# include <cstring>

extern "C" {

// definition lifted from SWIG-1.3.24/Lib/python/pyrun.swg
typedef struct {
  PyObject_HEAD
  void *ptr;
  const char *desc;
} PySwigObject;

}

namespace boost { namespace python {

template <typename T>
struct swig_arg
{
    static std::string desc;

    swig_arg(const char* name, const char* prefix="_p_")
    {
        desc = std::string(prefix) + name;
        converter::registry::insert(&extract, type_id<T>());
    }
 private:
    static void* extract(PyObject* op)
    {
        if (std::strcmp(op->ob_type->tp_name, "PySwigObject") != 0) return 0;

        PySwigObject* swig_obj_ptr = reinterpret_cast<PySwigObject*>(op);
        if (std::strcmp(swig_obj_ptr->desc, desc.c_str()) != 0) return 0;
        return swig_obj_ptr->ptr;
    }
};

template <typename T> std::string swig_arg<T>::desc;

}} // namespace boost::python

#define BOOST_PYTHON_SWIG_ARG(T) boost::python::swig_arg<T>(#T);

#endif // BOOST_ADAPTBX_SWIG_ARG_HPP
