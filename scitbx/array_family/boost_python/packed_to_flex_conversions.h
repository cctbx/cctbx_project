/// Conversion of versa<T, scitbx::matrix::packed_[ul]_accessor> to flex array.

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/array_family/boost_python/ref_flex_conversions.h>

#include <boost/python/object.hpp>
#include <boost/python/to_python_converter.hpp>

namespace scitbx { namespace af { namespace boost_python {

  struct packed_u_size_functor
  {
    static std::size_t get(std::size_t sz) {
      return dimension_from_packed_size(sz);
    }
  };

  template <typename ElementType, typename PackedAccessorType>
  struct versa_packed_to_flex
  {
    typedef versa<ElementType, PackedAccessorType> packed_versa_t;
    typedef versa<ElementType, flex_grid<> > flex_versa_t;

    static PyObject *convert(packed_versa_t const &a) {
      flex_grid<> acc(a.accessor().size_1d());
      flex_versa_t result(a, acc);
      return boost::python::incref(boost::python::object(result).ptr());
    }

    static const PyTypeObject *get_pytype() {
      return boost::python::converter::registered<
        flex_versa_t>::converters.to_python_target_type();
    }

    versa_packed_to_flex() {
      boost::python::to_python_converter<packed_versa_t,
                                         versa_packed_to_flex,
                                         true>();
    }
  };

  template <typename ElementType>
  struct default_packed_flex_conversions
  {
    default_packed_flex_conversions() {
      // Upper diagonal packed
      versa_packed_to_flex<ElementType, packed_u_accessor>();
      ref_from_flex<const_ref<ElementType, packed_u_accessor>,
                    packed_u_size_functor>();
      ref_from_flex<ref<ElementType, packed_u_accessor>,
                    packed_u_size_functor>();

      // Lower diagonal packed
      versa_packed_to_flex<ElementType, packed_l_accessor>();
      ref_from_flex<const_ref<ElementType, packed_l_accessor>,
                    packed_u_size_functor>();
      ref_from_flex<ref<ElementType, packed_l_accessor>,
                    packed_u_size_functor>();
    }

  };

}}}
