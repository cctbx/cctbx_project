// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Dec 21: Using iterator-based fftbx interface (rwgk)
     2001 Nov 03: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/fftbx/gridding.h>
#include <cctbx/fftbx/complex_to_complex_3d.h>
#include <cctbx/fftbx/real_to_complex_3d.h>
#include <cctbx/std_vector_bpl.h>
#include <cctbx/carray_bpl.h>

using namespace cctbx;

namespace {

  int adjust_gridding_2(int min_grid,
                                int max_prime) {
    return fftbx::adjust_gridding(min_grid, max_prime);
  }
  int adjust_gridding_3(int min_grid,
                                int max_prime,
                                int mandatory_factor) {
    return fftbx::adjust_gridding(min_grid, max_prime, mandatory_factor);
  }

  int3 adjust_gridding_triple_2(
    const int3& min_grid,
    int max_prime) {
    return fftbx::adjust_gridding_array(min_grid, max_prime);
  }

  int3 adjust_gridding_triple_3(
    const int3& min_grid,
    int max_prime,
    const int3& mandatory_factors) {
    return fftbx::adjust_gridding_array(min_grid, max_prime,
                                        mandatory_factors);
  }

  void throw_size_error() {
    PyErr_SetString(PyExc_RuntimeError, "Vector is too small.");
    throw boost::python::error_already_set();
  }

  void throw_index_error() {
    PyErr_SetString(PyExc_IndexError, "Index is out of range.");
    throw boost::python::error_already_set();
  }

  template <typename ElementType, typename CastElementType>
  struct v3d_accessor : public vecrefnd<ElementType, dimension<3> >
  {
    v3d_accessor() {}
    v3d_accessor(const int3& dim,
                 std::vector<ElementType>& vec,
                 bool resize_vector)
      : vecrefnd<ElementType, dimension<3> >((void*) 0, dim)
    {
      if (resize_vector) vec.resize(this->dim().size1d());
      m_begin = &(*(vec.begin()));
    }
    v3d_accessor(const int3& dim,
                 std::vector<CastElementType>& vec,
                 bool resize_vector
#if (defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
                 , bool MSVC_DUMMY = false
#endif
                )
      : vecrefnd<ElementType, dimension<3> >((void*) 0, dim)
    {
      std::size_t se = sizeof(ElementType);
      std::size_t sc = sizeof(CastElementType);
      if (se >= sc) {
        if (resize_vector) vec.resize(this->dim().size1d() * se / sc);
      }
      else {
        if (dim[2] % (sc / se)) {
          PyErr_SetString(PyExc_RuntimeError,
            "Vector not properly dimensioned:"
            " number of elements in third dimension must be even.");
          throw boost::python::error_already_set();
        }
        if (resize_vector) vec.resize(this->dim().size1d() * se / sc);
      }
      m_begin = reinterpret_cast<ElementType*>(&(*(vec.begin())));
    }
    ElementType
    getitem(const int3& I) const {
      if (!dim().is_valid_index(I)) throw_index_error();
      return operator()(I);
    }
    void setitem(const int3& I,
                 ElementType value) {
      if (!dim().is_valid_index(I)) throw_index_error();
      operator()(I) = value;
    }
  };

  typedef v3d_accessor<double, std::complex<double> > vd3d_accessor;
  typedef v3d_accessor<std::complex<double>, double>  vc3d_accessor;

  void cc_forward_complex(fftbx::complex_to_complex<double>& fft,
                          std::vector<std::complex<double> >& vec) {
    if (vec.size() < fft.N()) throw_size_error();
    fft.forward(vec.begin());
  }
  void cc_backward_complex(fftbx::complex_to_complex<double>& fft,
                           std::vector<std::complex<double> >& vec) {
    if (vec.size() < fft.N()) throw_size_error();
    fft.backward(vec.begin());
  }
  void cc_forward_real(fftbx::complex_to_complex<double>& fft,
                       std::vector<double>& vec) {
    if (vec.size() < 2 * fft.N()) throw_size_error();
    fft.forward(vec.begin());
  }
  void cc_backward_real(fftbx::complex_to_complex<double>& fft,
                        std::vector<double>& vec) {
    if (vec.size() < 2 * fft.N()) throw_size_error();
    fft.backward(vec.begin());
  }

  void rc_forward_complex(fftbx::real_to_complex<double>& fft,
                          std::vector<std::complex<double> >& vec) {
    if (vec.size() < fft.Ncomplex()) throw_size_error();
    fft.forward(vec.begin());
  }
  void rc_backward_complex(fftbx::real_to_complex<double>& fft,
                           std::vector<std::complex<double> >& vec) {
    if (vec.size() < fft.Ncomplex()) throw_size_error();
    fft.backward(vec.begin());
  }
  void rc_forward_real(fftbx::real_to_complex<double>& fft,
                       std::vector<double>& vec) {
    if (vec.size() < fft.Mreal()) throw_size_error();
    fft.forward(vec.begin());
  }
  void rc_backward_real(fftbx::real_to_complex<double>& fft,
                        std::vector<double>& vec) {
    if (vec.size() < fft.Mreal()) throw_size_error();
    fft.backward(vec.begin());
  }

  void cc_3d_forward_complex(fftbx::complex_to_complex_3d<double>& fft,
                     vc3d_accessor map) {
    fft.forward(map);
  }
  void cc_3d_backward_complex(fftbx::complex_to_complex_3d<double>& fft,
                      vc3d_accessor map) {
    fft.backward(map);
  }
  void cc_3d_forward_real(fftbx::complex_to_complex_3d<double>& fft,
                          vd3d_accessor map) {
    fft.forward(map);
  }
  void cc_3d_backward_real(fftbx::complex_to_complex_3d<double>& fft,
                          vd3d_accessor map) {
    fft.backward(map);
  }

  void rc_3d_forward_complex(fftbx::real_to_complex_3d<double>& fft,
                             vc3d_accessor map) {
    fft.forward(map);
  }
  void rc_3d_backward_complex(fftbx::real_to_complex_3d<double>& fft,
                              vc3d_accessor map) {
    fft.backward(map);
  }
  void rc_3d_forward_real(fftbx::real_to_complex_3d<double>& fft,
                          vd3d_accessor map) {
    fft.forward(map);
  }
  void rc_3d_backward_real(fftbx::real_to_complex_3d<double>& fft,
                           vd3d_accessor map) {
    fft.backward(map);
  }

#   include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    (void) python::wrap_std_vector(this_module,
      "vector_of_int", int());
    python::class_builder<
      std::vector<double>,
      python::std_vector_wrapper<double> >
    py_vector_of_double =
    python::wrap_std_vector(this_module,
      "vector_of_double",
      double());
    python::export_converters(py_vector_of_double);
    python::class_builder<
      std::vector<std::complex<double> >,
      python::std_vector_wrapper<std::complex<double> > >
    py_vector_of_complex =
    python::wrap_std_vector(this_module,
      "vector_of_complex",
      std::complex<double>());
    python::export_converters(py_vector_of_complex);

    class_builder<vd3d_accessor>
    py_vd3d_accessor(this_module, "vd3d_accessor");
    class_builder<vc3d_accessor>
    py_vc3d_accessor(this_module, "vc3d_accessor");

    class_builder<fftbx::factorization>
    py_factorization(this_module, "factorization");
    class_builder<fftbx::complex_to_complex<double> >
    py_complex_to_complex(this_module, "complex_to_complex");
    class_builder<fftbx::real_to_complex<double> >
    py_real_to_complex(this_module, "real_to_complex");
    class_builder<fftbx::complex_to_complex_3d<double> >
    py_complex_to_complex_3d(this_module, "complex_to_complex_3d");
    class_builder<fftbx::real_to_complex_3d<double> >
    py_real_to_complex_3d(this_module, "real_to_complex_3d");

    py_complex_to_complex.declare_base(
      py_factorization, boost::python::without_downcast);
    py_real_to_complex.declare_base(
      py_factorization, boost::python::without_downcast);

    this_module.def(adjust_gridding_2, "adjust_gridding");
    this_module.def(adjust_gridding_3, "adjust_gridding");
    this_module.def(adjust_gridding_triple_2, "adjust_gridding_triple");
    this_module.def(adjust_gridding_triple_3, "adjust_gridding_triple");

    py_vd3d_accessor.def(constructor<>());
    py_vd3d_accessor.def(constructor<
      const int3&,
      std::vector<double>&,
      bool>());
    py_vd3d_accessor.def(constructor<
      const int3&,
      std::vector<std::complex<double> >&,
      bool>());
    py_vd3d_accessor.def(&vd3d_accessor::getitem, "__getitem__");
    py_vd3d_accessor.def(&vd3d_accessor::setitem, "__setitem__");

    py_vc3d_accessor.def(constructor<>());
    py_vc3d_accessor.def(constructor<
      const int3&,
      std::vector<std::complex<double> >&,
      bool>());
    py_vc3d_accessor.def(constructor<
      const int3&,
      std::vector<double>&,
      bool>());
    py_vc3d_accessor.def(&vc3d_accessor::getitem, "__getitem__");
    py_vc3d_accessor.def(&vc3d_accessor::setitem, "__setitem__");

    py_factorization.def(constructor<>());
    py_factorization.def(constructor<std::size_t, bool>());
    py_factorization.def(&fftbx::factorization::N, "N");
    py_factorization.def(&fftbx::factorization::Factors, "Factors");

    py_complex_to_complex.def(constructor<>());
    py_complex_to_complex.def(constructor<std::size_t>());
    py_complex_to_complex.def(
      &fftbx::complex_to_complex<double>::WA, "WA");
    py_complex_to_complex.def(cc_forward_complex, "forward");
    py_complex_to_complex.def(cc_backward_complex, "backward");
    py_complex_to_complex.def(cc_forward_real, "forward");
    py_complex_to_complex.def(cc_backward_real, "backward");

    py_real_to_complex.def(constructor<>());
    py_real_to_complex.def(constructor<std::size_t>());
    py_real_to_complex.def(&fftbx::real_to_complex<double>::Nreal, "Nreal");
    py_real_to_complex.def(&fftbx::real_to_complex<double>::Mreal, "Mreal");
    py_real_to_complex.def(
      &fftbx::real_to_complex<double>::Ncomplex, "Ncomplex");
    py_real_to_complex.def(&fftbx::real_to_complex<double>::WA, "WA");
    py_real_to_complex.def(rc_forward_complex, "forward");
    py_real_to_complex.def(rc_backward_complex, "backward");
    py_real_to_complex.def(rc_forward_real, "forward");
    py_real_to_complex.def(rc_backward_real, "backward");

    py_complex_to_complex_3d.def(constructor<>());
    py_complex_to_complex_3d.def(
      constructor<std::size_t, std::size_t, std::size_t>());
    py_complex_to_complex_3d.def(
      constructor<const int3&>());
    py_complex_to_complex_3d.def(
      &fftbx::complex_to_complex_3d<double>::N, "N");
    py_complex_to_complex_3d.def(cc_3d_forward_complex, "forward");
    py_complex_to_complex_3d.def(cc_3d_backward_complex, "backward");
    py_complex_to_complex_3d.def(cc_3d_forward_real, "forward");
    py_complex_to_complex_3d.def(cc_3d_backward_real, "backward");

    py_real_to_complex_3d.def(constructor<>());
    py_real_to_complex_3d.def(
      constructor<std::size_t, std::size_t, std::size_t>());
    py_real_to_complex_3d.def(
      constructor<const int3&>());
    py_real_to_complex_3d.def(
      &fftbx::real_to_complex_3d<double>::Nreal, "Nreal");
    py_real_to_complex_3d.def(
      &fftbx::real_to_complex_3d<double>::Mreal, "Mreal");
    py_real_to_complex_3d.def(
      &fftbx::real_to_complex_3d<double>::Ncomplex, "Ncomplex");
    py_real_to_complex_3d.def(rc_3d_forward_complex, "forward");
    py_real_to_complex_3d.def(rc_3d_backward_complex, "backward");
    py_real_to_complex_3d.def(rc_3d_forward_real, "forward");
    py_real_to_complex_3d.def(rc_3d_backward_real, "backward");
  }

}

BOOST_PYTHON_MODULE_INIT(fftbx)
{
  boost::python::module_builder this_module("fftbx");
  init_module(this_module);
}
