// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/class_builder.hpp>
#include <cctbx/error.h>
#include <cctbx/miller.h>
#include <cctbx/hendrickson_lattman.h>
#include <cctbx/sftbx/xray_scatterer.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace af {

  namespace {

    struct getstate_manager
    {
      getstate_manager(std::size_t a_size, std::size_t size_per_element)
      {
        str_capacity = a_size * size_per_element + 50; // extra space for a_size
        str_obj = PyString_FromStringAndSize(
          0, static_cast<int>(str_capacity + 100)); // extra space for safety
        str_begin = PyString_AS_STRING(str_obj);
        str_end = str_begin;
        sprintf(str_end, "%lu", static_cast<unsigned long>(a_size));
        while (*str_end) str_end++;
        *str_end++ = ',';
      };

      void advance()
      {
        while (*str_end) str_end++;
        *str_end++ = ',';
        cctbx_assert(str_end - str_begin <= str_capacity);
      }

      PyObject* finalize()
      {
        cctbx_assert(_PyString_Resize(&str_obj,
          static_cast<int>(str_end - str_begin)) == 0);
        return str_obj;
      }

      std::size_t str_capacity;
      PyObject* str_obj;
      char* str_begin;
      char* str_end;
    };

    struct setstate_manager
    {
      setstate_manager(std::size_t a_size, PyObject* state)
      {
        cctbx_assert(a_size == 0);
        str_ptr = PyString_AsString(state);
        cctbx_assert(str_ptr != 0);
        cctbx_assert(sscanf(str_ptr, "%lu", &a_capacity) == 1);
        while (*str_ptr != ',') str_ptr++;
        str_ptr++;
      };

      void advance()
      {
        while (*str_ptr != ',') str_ptr++;
        str_ptr++;
      }

      void finalize()
      {
        cctbx_assert(*str_ptr == 0);
      }

      char* str_ptr;
      std::size_t a_capacity;
    };

    struct bool_picklers
    {
      static
      boost::python::ref
      getstate(shared<bool> const& a)
      {
        getstate_manager mgr(a.size(), 1);
        for(std::size_t i=0;i<a.size();i++) {
          if (a[i]) *mgr.str_end++ = '1';
          else      *mgr.str_end++ = '0';
        }
        return boost::python::ref(mgr.finalize());
      }

      static
      void
      setstate(shared<bool>& a, boost::python::ref state)
      {
        setstate_manager mgr(a.size(), state.get());
        a.reserve(mgr.a_capacity);
        for(std::size_t i=0;i<mgr.a_capacity;i++) {
          if (*mgr.str_ptr++ == '1') a.push_back(true);
          else                       a.push_back(false);
        }
        mgr.finalize();
      }
    };

    template <typename ElementType>
    struct num_picklers
    {
      static
      boost::python::ref
      getstate(
        shared<ElementType> const& a,
        std::size_t size_per_element,
        const char* fmt)
      {
        getstate_manager mgr(a.size(), size_per_element);
        for(std::size_t i=0;i<a.size();i++) {
          sprintf(mgr.str_end, fmt, a[i]);
          mgr.advance();
        }
        return boost::python::ref(mgr.finalize());
      }

      static
      void
      setstate(
        shared<ElementType>& a,
        boost::python::ref state,
        const char* fmt)
      {
        setstate_manager mgr(a.size(), state.get());
        a.reserve(mgr.a_capacity);
        for(std::size_t i=0;i<mgr.a_capacity;i++) {
          ElementType val;
          cctbx_assert(sscanf(mgr.str_ptr, fmt, &val) == 1);
          mgr.advance();
          a.push_back(val);
        }
        mgr.finalize();
      }
    };

    template <typename ElementType>
    struct complex_picklers
    {
      static
      boost::python::ref
      getstate(
        shared<std::complex<ElementType> > const& a,
        std::size_t size_per_element,
        const char* fmt)
      {
        getstate_manager mgr(a.size(), 2 * size_per_element);
        for(std::size_t i=0;i<a.size();i++) {
          sprintf(mgr.str_end, fmt, a[i].real());
          mgr.advance();
          sprintf(mgr.str_end, fmt, a[i].imag());
          mgr.advance();
        }
        return boost::python::ref(mgr.finalize());
      }

      static
      void
      setstate(
        shared<std::complex<ElementType> >& a,
        boost::python::ref state,
        const char* fmt)
      {
        setstate_manager mgr(a.size(), state.get());
        a.reserve(mgr.a_capacity);
        for(std::size_t i=0;i<mgr.a_capacity;i++) {
          ElementType val[2];
          for(std::size_t j=0;j<2;j++) {
            cctbx_assert(sscanf(mgr.str_ptr, fmt, val+j) == 1);
            mgr.advance();
          }
          a.push_back(std::complex<ElementType>(val[0], val[1]));
        }
        mgr.finalize();
      }
    };

    struct miller_index_picklers
    {
      static
      boost::python::ref
      getstate(shared<miller::Index> const& a)
      {
        getstate_manager mgr(a.size(), 3 * 6);
        for(std::size_t i=0;i<a.size();i++) {
          miller::Index const& h = a[i];
          for(std::size_t j=0;j<3;j++) {
            sprintf(mgr.str_end, "%d", h[j]);
            mgr.advance();
          }
        }
        return boost::python::ref(mgr.finalize());
      }

      static
      void
      setstate(shared<miller::Index>& a, boost::python::ref state)
      {
        setstate_manager mgr(a.size(), state.get());
        a.reserve(mgr.a_capacity);
        for(std::size_t i=0;i<mgr.a_capacity;i++) {
          int h[3];
          for(std::size_t j=0;j<3;j++) {
            cctbx_assert(sscanf(mgr.str_ptr, "%d", h+j) == 1);
            mgr.advance();
          }
          a.push_back(miller::Index(h));
        }
        mgr.finalize();
      }
    };

    template <typename FloatType>
    struct hendrickson_lattman_picklers
    {
      typedef hendrickson_lattman<FloatType> hl_type;

      static
      boost::python::ref
      getstate(
        shared<hl_type> const& a,
        std::size_t size_per_element,
        const char* fmt)
      {
        getstate_manager mgr(a.size(), 4 * size_per_element);
        for(std::size_t i=0;i<a.size();i++) {
          af::tiny<FloatType, 4> const& c = a[i].array();
          for(std::size_t j=0;j<4;j++) {
            sprintf(mgr.str_end, fmt, c[j]);
            mgr.advance();
          }
        }
        return boost::python::ref(mgr.finalize());
      }

      static
      void
      setstate(
        shared<hl_type>& a,
        boost::python::ref state,
        const char* fmt)
      {
        setstate_manager mgr(a.size(), state.get());
        a.reserve(mgr.a_capacity);
        for(std::size_t i=0;i<mgr.a_capacity;i++) {
          FloatType c[4];
          for(std::size_t j=0;j<4;j++) {
            cctbx_assert(sscanf(mgr.str_ptr, fmt, c+j) == 1);
            mgr.advance();
          }
          a.push_back(hl_type(c));
        }
        mgr.finalize();
      }
    };

    struct make_pickle_string
    {
      std::string buffer;

      make_pickle_string& operator<<(std::string const& val)
      {
        buffer += val + '\0' + ',';
        return *this;
      }

      make_pickle_string& operator<<(const char* val)
      {
        return *this << std::string(val);
      }

      make_pickle_string& operator<<(bool const& val)
      {
        if (val) buffer += "1,";
        else     buffer += "0,";
        return *this;
      }

      make_pickle_string& operator<<(int const& val)
      {
        char buf[64];
        sprintf(buf, "%d,", val);
        buffer += buf;
        return *this;
      }

      make_pickle_string& operator<<(unsigned int const& val)
      {
        char buf[64];
        sprintf(buf, "%u,", val);
        buffer += buf;
        return *this;
      }

      make_pickle_string& operator<<(long const& val)
      {
        char buf[64];
        sprintf(buf, "%ld,", val);
        buffer += buf;
        return *this;
      }

      make_pickle_string& operator<<(unsigned long const& val)
      {
        char buf[64];
        sprintf(buf, "%lu,", val);
        buffer += buf;
        return *this;
      }

      make_pickle_string& operator<<(float const& val)
      {
        char buf[64];
        sprintf(buf, "%.6g,", val);
        buffer += buf;
        return *this;
      }

      make_pickle_string& operator<<(double const& val)
      {
        char buf[64];
        sprintf(buf, "%.12g,", val);
        buffer += buf;
        return *this;
      }

      template <typename FloatType>
      make_pickle_string& operator<<(std::complex<FloatType> const& val)
      {
        return *this << val.real() << val.imag();
      }

      template <typename ElementType, typename AccessorType>
      make_pickle_string& operator<<(
        af::const_ref<ElementType, AccessorType> const& vals)
      {
        for(std::size_t i=0;i<vals.size();i++) {
          *this << vals[i];
        }
        return *this;
      }
    };

    struct read_from_pickle_string
    {
      const char* str_ptr;

      read_from_pickle_string(PyObject* str_obj)
      : str_ptr(PyString_AsString(str_obj))
      {
        cctbx_assert(str_ptr != 0);
      }

      void assert_end() const
      {
        cctbx_assert(*str_ptr == 0);
      }

      read_from_pickle_string& advance()
      {
        while (*str_ptr != ',') str_ptr++;
        str_ptr++;
        return *this;
      }

      read_from_pickle_string& operator>>(std::string& val)
      {
        val = std::string(str_ptr);
        while (*str_ptr) str_ptr++;
        return advance();
      }

      read_from_pickle_string& operator>>(bool& val)
      {
        *str_ptr == '1' ? val = true : val = false;
        return advance();
      }

      read_from_pickle_string& operator>>(int& val)
      {
        cctbx_assert(sscanf(str_ptr, "%d", &val) == 1);
        return advance();
      }

      read_from_pickle_string& operator>>(unsigned int& val)
      {
        cctbx_assert(sscanf(str_ptr, "%u", &val) == 1);
        return advance();
      }

      read_from_pickle_string& operator>>(long& val)
      {
        cctbx_assert(sscanf(str_ptr, "%ld", &val) == 1);
        return advance();
      }

      read_from_pickle_string& operator>>(unsigned long& val)
      {
        cctbx_assert(sscanf(str_ptr, "%lu", &val) == 1);
        return advance();
      }

      read_from_pickle_string& operator>>(float& val)
      {
        cctbx_assert(sscanf(str_ptr, "%g", &val) == 1);
        return advance();
      }

      read_from_pickle_string& operator>>(double& val)
      {
        cctbx_assert(sscanf(str_ptr, "%lg", &val) == 1);
        return advance();
      }

      template <typename FloatType>
      read_from_pickle_string& operator>>(std::complex<FloatType>& val)
      {
        FloatType r, i;
        *this >> r >> i;
        val = std::complex<FloatType>(r, i);
        return *this;
      }

      template <typename ElementType, typename AccessorType>
      read_from_pickle_string& operator>>(
        af::ref<ElementType, AccessorType> vals)
      {
        for(std::size_t i=0;i<vals.size();i++) {
          *this >> vals[i];
        }
        return *this;
      }
    };

    template <typename FloatType, typename CAASF_Type>
    struct xray_scatterer_picklers
    {
      typedef sftbx::XrayScatterer<FloatType, CAASF_Type> xs_type;

      static
      boost::python::ref
      getstate(shared<xs_type> const& a)
      {
        make_pickle_string result;
        result << a.size();
        for(std::size_t i=0;i<a.size();i++) {
          xs_type const& site = a[i];
          result << site.Label()
                 << site.CAASF().Label()
                 << site.fpfdp()
                 << site.Coordinates().const_ref()
                 << site.Occ()
                 << site.isAnisotropic();
          if (site.isAnisotropic()) result << site.Uaniso().const_ref();
          else                      result << site.Uiso();
        }
        return boost::python::make_ref(result.buffer);
      }

      static
      void
      setstate(shared<xs_type>& a, boost::python::ref state)
      {
        cctbx_assert(a.size() == 0);
        read_from_pickle_string inp(state.get());
        std::size_t a_capacity;
        inp >> a_capacity;
        a.reserve(a_capacity);
        std::string lbl, sf_lbl;
        std::complex<FloatType> fpfdp;
        fractional<FloatType> coor;
        FloatType occ;
        bool is_aniso;
        for(std::size_t i=0;i<a_capacity;i++) {
          inp >> lbl >> sf_lbl >> fpfdp >> coor.ref() >> occ >> is_aniso;
          if (is_aniso) {
            af::tiny<FloatType, 6> u_aniso;
            inp >> u_aniso.ref();
            a.push_back(
              xs_type(lbl, CAASF_Type(sf_lbl), fpfdp, coor, occ, u_aniso));
          }
          else {
            FloatType u_iso;
            inp >> u_iso;
            a.push_back(
              xs_type(lbl, CAASF_Type(sf_lbl), fpfdp, coor, occ, u_iso));
          }
        }
        inp.assert_end();
      }
    };

  } // namespace <anonymous>

  boost::python::ref shared_bool_getstate(shared<bool> const& a)
  {
    return bool_picklers::getstate(a);
  }

  void shared_bool_setstate(shared<bool>& a, boost::python::ref state)
  {
    bool_picklers::setstate(a, state);
  }

  boost::python::ref shared_int_getstate(shared<int> const& a)
  {
    return num_picklers<int>::getstate(a, 12, "%d");
  }

  void shared_int_setstate(shared<int>& a, boost::python::ref state)
  {
    num_picklers<int>::setstate(a, state, "%d");
  }

  boost::python::ref shared_long_getstate(shared<long> const& a)
  {
    return num_picklers<long>::getstate(a, 21, "%ld");
  }

  void shared_long_setstate(shared<long>& a, boost::python::ref state)
  {
    num_picklers<long>::setstate(a, state, "%ld");
  }

  boost::python::ref shared_float_getstate(shared<float> const& a)
  {
    return num_picklers<float>::getstate(a, 13, "%.6g");
  }

  void shared_float_setstate(shared<float>& a, boost::python::ref state)
  {
    num_picklers<float>::setstate(a, state, "%g");
  }

  boost::python::ref shared_double_getstate(shared<double> const& a)
  {
    return num_picklers<double>::getstate(a, 19, "%.12g");
  }

  void shared_double_setstate(shared<double>& a, boost::python::ref state)
  {
    num_picklers<double>::setstate(a, state, "%lg");
  }

  boost::python::ref shared_complex_double_getstate(
    shared<std::complex<double> > const& a)
  {
    return complex_picklers<double>::getstate(a, 19, "%.12g");
  }

  void shared_complex_double_setstate(
    shared<std::complex<double> >& a,
    boost::python::ref state)
  {
    complex_picklers<double>::setstate(a, state, "%lg");
  }

  boost::python::ref shared_miller_index_getstate(
    shared<miller::Index> const& a)
  {
    return miller_index_picklers::getstate(a);
  }

  void shared_miller_index_setstate(
    shared<miller::Index>& a,
    boost::python::ref state)
  {
    miller_index_picklers::setstate(a, state);
  }

  boost::python::ref shared_hendrickson_lattman_double_getstate(
    shared<hendrickson_lattman<double> > const& a)
  {
    return hendrickson_lattman_picklers<double>::getstate(a, 19, "%.12g");
  }

  void shared_hendrickson_lattman_double_setstate(
    shared<hendrickson_lattman<double> >& a,
    boost::python::ref state)
  {
    hendrickson_lattman_picklers<double>::setstate(a, state, "%lg");
  }

  boost::python::ref shared_xray_scatterer_double_wk1995_getstate(
    shared<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> > const& a)
  {
    return xray_scatterer_picklers<double, eltbx::CAASF_WK1995>::getstate(a);
  }

  void shared_xray_scatterer_double_wk1995_setstate(
    shared<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> >& a,
    boost::python::ref state)
  {
    xray_scatterer_picklers<double, eltbx::CAASF_WK1995>::setstate(a, state);
  }

}} // namespace cctbx::af
