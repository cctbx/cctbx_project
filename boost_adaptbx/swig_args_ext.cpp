#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/lvalue_from_pytype.hpp>

#include <iostream>
#include <exception>
#include <cstring>

//############################################################################
// copy of SWIG-1.3.24/Examples/python/class/example.h
//############################################################################
/* File : example.h */

class Shape {
public:
  Shape() {
    nshapes++;
  }
  virtual ~Shape() {
    nshapes--;
  };
  double  x, y;
  void    move(double dx, double dy);
  virtual double area(void) = 0;
  virtual double perimeter(void) = 0;
  static  int nshapes;
};

class Circle : public Shape {
private:
  double radius;
public:
  Circle(double r) : radius(r) { };
  virtual double area(void);
  virtual double perimeter(void);
};

class Square : public Shape {
private:
  double width;
public:
  Square(double w) : width(w) { };
  virtual double area(void);
  virtual double perimeter(void);
};


//############################################################################
// this should go somewhere into the boost/boost/python include directory
//############################################################################
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

//############################################################################
// user code
//############################################################################

void
show_circle(Circle* circle)
{
  std::cout << "x: " << circle->x << std::endl;
  std::cout << "y: " << circle->y << std::endl;
  std::cout << "area: " << circle->area() << std::endl;
}

void
show_square(Square* square)
{
  std::cout << "x: " << square->x << std::endl;
  std::cout << "y: " << square->y << std::endl;
  std::cout << "area: " << square->area() << std::endl;
  std::cout << "area: " << square->perimeter() << std::endl;
}

BOOST_PYTHON_MODULE(boost_python_swig_args_ext)
{
  using namespace boost::python;
  BOOST_PYTHON_SWIG_ARG(Circle);
  BOOST_PYTHON_SWIG_ARG(Square);
  def("show", show_circle);
  def("show", show_square);
}
