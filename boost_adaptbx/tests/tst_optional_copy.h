#include <boost_adaptbx/optional_copy.h>
#include <boost_adaptbx/error_utils.h>

namespace boost_adaptbx { namespace tst_optional_copy {

  template <typename ValueType>
  optional_copy<ValueType> const&
  oc_const_ref(optional_copy<ValueType> const& o)
  {
    return o;
  }

  template <typename ValueType>
  void
  exercise(
    ValueType const& v1,
    ValueType const& v4,
    bool value_is_shared)
  {
    typedef optional_copy<ValueType> o;
    {
      o o1;
      ASSERTBX(!o1);
      ASSERTBX(o1.get() == 0);
      o o2(o1);
      ASSERTBX(o2.get() == 0);
      o2 = o1;
      ASSERTBX(o2.get() == 0);
      o2 = v1;
      ASSERTBX(o2);
      ASSERTBX((*o2.get())[0] == 1);
      ASSERTBX(o2->begin() == o2.get()->begin());
      (*o2)[0] = 2;
      ASSERTBX((*oc_const_ref(o2))[0] == 2);
      (*o2)[0] = 1;
      ASSERTBX(o2[0] == 1);
      o1 = o2;
      ASSERTBX(o1[0] == 1);
      o2.release();
      ASSERTBX(o2.get() == 0);
      ASSERTBX(o1[0] == 1);
      ASSERTBX(oc_const_ref(o2).get() == 0);
      ASSERTBX(o1->begin() == oc_const_ref(o1).get()->begin());
      ASSERTBX(oc_const_ref(o1)[0] == 1);
    }
    {
      o o1(v1);
      ASSERTBX(o1.get() != 0);
      ASSERTBX(o1[0] == 1);
      o o2(o1);
      ASSERTBX(o2.get() != 0);
      ASSERTBX(o2.get() != o1.get());
      ASSERTBX(o2[0] == 1);
      o2[0] = 2;
      ASSERTBX(o1[0] == (value_is_shared ? 2 : 1));
      ASSERTBX(o2[0] == 2);
      o1 = o2;
      ASSERTBX(o1[0] == 2);
      ASSERTBX(o2[0] == 2);
      o1[0] = 3;
      ASSERTBX(o1[0] == 3);
      ASSERTBX(o2[0] == (value_is_shared ? 3 : 2));
      o2 = v4;
      ASSERTBX(o1[0] == 3);
      ASSERTBX(o2[0] == 4);
    }
  };

}} // namespace boost_adaptbx::tst_optional_copy
