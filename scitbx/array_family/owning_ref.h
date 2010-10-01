#ifndef SCITBX_ARRAY_FAMILY_OWNING_REF_H
#define SCITBX_ARRAY_FAMILY_OWNING_REF_H

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>

namespace scitbx { namespace af {

  /// Like af::ref but it owns the versa array it refers to
  /** This array shall never be resized after its creation. Otherwise the
      ref would get stale.
   */
  template <typename T, class AccessorType=typename af::versa<T>::accessor_type>
  class ref_owning_versa : public af::ref<T, AccessorType>
  {
  public:
    typedef typename AccessorType::index_value_type index_value_type;

    /// Uninitialised
    ref_owning_versa() {}

    /// Takes ownership of a and construct a reference to it
    ref_owning_versa(versa<T, AccessorType> const &array)
      : a(array)
    {
      init();
    }

    //@{
    /** Construct the owned array and a reference to it */
    ref_owning_versa(AccessorType const &ac)
      : a(ac)
    {
      init();
    }

    ref_owning_versa(index_value_type n0)
      : a(n0)
    {
      init();
    }

    ref_owning_versa(AccessorType const &ac, T const &x)
      : a(ac, x)
    {
      init();
    }

    ref_owning_versa(index_value_type n0, T const &x)
      : a(n0, x)
    {
      init();
    }

    template <typename FunctorType>
    ref_owning_versa(AccessorType const& ac,
                     init_functor<FunctorType> const& ftor)
      : a(ac, ftor)
    {
      init();
    }

    template <typename FunctorType>
    ref_owning_versa(index_value_type n0,
                     init_functor<FunctorType> const& ftor)
      : a(n0, ftor)
    {
      init();
    }
    //@}

    /// Release ownership of its array and take ownership of other,
    /// building a reference to it
    ref_owning_versa &operator=(versa<T, AccessorType> const &other) {
      a = other;
      af::ref<T, AccessorType>::operator=(a.ref());
      return *this;
    }

    /// The owned array, returned by value for transfer of ownership
    versa<T, AccessorType> array() const { return a; }

    /// The mere ref equal to this
    af::ref<T, AccessorType> ref() const { return *this; }

    /// A const ref to the owned array
    af::const_ref<T, AccessorType> const_ref() const { return *this; }

  private:
    versa<T, AccessorType> a;

    void init() {
      af::ref<T, AccessorType>::operator=(a.ref());
    }
  };


  /// Like af::ref but it owns the shared array it refers to
  /** This array shall never be resized after its creation. Otherwise the
   ref would get stale.
   */
  template <typename T>
  class ref_owning_shared : public af::ref<T>
  {
  public:
    typedef typename shared<T>::size_type size_type;

    /// Uninitialised
    ref_owning_shared() {}

    /// Takes ownership of a and construct a reference to it
    ref_owning_shared(shared<T> const &array)
      : a(array)
    {
      init();
    }

    //@{
    /** Construct the owned array and a reference to it */
    ref_owning_shared(size_type sz)
      : a(sz)
    {
      init();
    }

    ref_owning_shared(reserve const& sz)
      : a(sz)
    {
      init();
    }

    ref_owning_shared(size_type sz, T const &x)
      : a(sz, x)
    {
      init();
    }

    template <typename FunctorType>
    ref_owning_shared(size_type const& sz,
                      init_functor<FunctorType> const& ftor)
      : a(sz, ftor)
    {
      init();
    }
    //@}

    /// Release ownership of its array and take ownership of other,
    /// building a reference to it
    ref_owning_shared &operator=(shared<T> const &other) {
      a = other;
      af::ref<T>::operator=(a.ref());
      return *this;
    }

    /// The owned array, returned by value for transfer of ownership
    shared<T> array() const { return a; }

    /// The mere ref equal to this
    af::ref<T> ref() const { return *this; }

    /// A const ref to the owned array
    af::const_ref<T> const_ref() const { return *this; }

  private:
    shared<T> a;

    void init() {
      af::ref<T>::operator=(a.ref());
    }
  };



}}



#endif // GUARD
