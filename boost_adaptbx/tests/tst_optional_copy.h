#include <tbxx/optional_copy.hpp>
#include <tbxx/error_utils.hpp>

namespace boost_adaptbx { namespace tst_optional_copy {

  template <typename ContainerType>
  tbxx::optional_container<ContainerType> const&
  oc_const_ref(tbxx::optional_container<ContainerType> const& o)
  {
    return o;
  }

  template <typename ContainerType>
  void
  exercise(
    ContainerType const& v1,
    ContainerType const& v4,
    bool value_is_shared)
  {
    typedef tbxx::optional_container<ContainerType> o;
    {
      o o1;
      TBXX_ASSERT(!o1);
      TBXX_ASSERT(o1.get() == 0);
      o o2(o1);
      TBXX_ASSERT(o2.get() == 0);
      o2 = o1;
      TBXX_ASSERT(o2.get() == 0);
      o2 = v1;
      TBXX_ASSERT(o2);
      TBXX_ASSERT((*o2.get())[0] == 1);
      TBXX_ASSERT(o2->begin() == o2.get()->begin());
      (*o2)[0] = 2;
      TBXX_ASSERT((*oc_const_ref(o2))[0] == 2);
      (*o2)[0] = 1;
      TBXX_ASSERT(o2[0] == 1);
      o1 = o2;
      TBXX_ASSERT(o1[0] == 1);
      o2.release();
      TBXX_ASSERT(o2.get() == 0);
      TBXX_ASSERT(o1[0] == 1);
      TBXX_ASSERT(oc_const_ref(o2).get() == 0);
      TBXX_ASSERT(o1->begin() == oc_const_ref(o1).get()->begin());
      TBXX_ASSERT(oc_const_ref(o1)[0] == 1);
    }
    {
      o o1(v1);
      TBXX_ASSERT(o1.get() != 0);
      TBXX_ASSERT(o1[0] == 1);
      o o2(o1);
      TBXX_ASSERT(o2.get() != 0);
      TBXX_ASSERT(o2.get() != o1.get());
      TBXX_ASSERT(o2[0] == 1);
      o2[0] = 2;
      TBXX_ASSERT(o1[0] == (value_is_shared ? 2 : 1));
      TBXX_ASSERT(o2[0] == 2);
      o1 = o2;
      TBXX_ASSERT(o1[0] == 2);
      TBXX_ASSERT(o2[0] == 2);
      o1[0] = 3;
      TBXX_ASSERT(o1[0] == 3);
      TBXX_ASSERT(o2[0] == (value_is_shared ? 3 : 2));
      o2 = v4;
      TBXX_ASSERT(o1[0] == 3);
      TBXX_ASSERT(o2[0] == 4);
    }
  };

}} // namespace boost_adaptbx::tst_optional_copy
