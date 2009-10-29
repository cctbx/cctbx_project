#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_PASSING_FLEX_AS_ARGUMENT_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_PASSING_FLEX_AS_ARGUMENT_H

#include <scitbx/array_family/boost_python/utils.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>

namespace scitbx { namespace af { namespace boost_python {

template<class ElementType>
class flex_1d
  : public af::versa<ElementType, af::flex_grid<> >::base_array_type
{
  public:
    typedef af::versa<ElementType, af::flex_grid<> > flex_array_type;
    typedef typename flex_array_type::base_array_type base_type;

    flex_1d(flex_array_type &array) :
      base_type(array.as_base_array()),
      a(array)
    {
      SCITBX_ASSERT(array.accessor().nd() == 1 && array.accessor().is_0_based())
                   (array.accessor().nd());
    }

    ~flex_1d() { a.resize(af::flex_grid<>(this->size())); }

  private:
    flex_array_type &a;
};


template <typename ElementType>
struct flex_1d_from_flex
{
  typedef flex_1d<ElementType> target_type;
  typedef typename target_type::flex_array_type flex_type;

  flex_1d_from_flex() {
    using namespace boost::python;
    converter::registry::push_back(
      &convertible,
      &construct,
      type_id<target_type>()
      #ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
      , &converter::expected_pytype_for_arg<flex_type>::get_pytype
      #endif
    );
  }

  static void* convertible(PyObject *obj_ptr) {
    using namespace boost::python;
    object obj(borrowed(obj_ptr));
    extract<flex_type&> flex_proxy(obj);
    return flex_proxy.check() ? obj_ptr : 0;
  }

  static void construct(
    PyObject *obj_ptr,
    boost::python::converter::rvalue_from_python_stage1_data *data)
  {
    using namespace boost::python;
    object obj(borrowed(obj_ptr));
    flex_type &a = extract<flex_type&>(obj)();
    if (!a.check_shared_size()) raise_shared_size_mismatch();
    typedef converter::rvalue_from_python_storage<target_type> rvalue_t;
    void *storage = reinterpret_cast<rvalue_t*>(data)->storage.bytes;
    new (storage) target_type(a);
    data->convertible = storage;
  }
};

}}} // namespace scitbx::af::boost_python

#endif // GUARD
