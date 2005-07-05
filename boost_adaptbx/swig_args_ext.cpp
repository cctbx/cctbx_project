#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include "swig_arg.hpp"

#include <iostream>

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
