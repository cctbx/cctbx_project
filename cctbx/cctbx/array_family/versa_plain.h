// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_VERSA_PLAIN_H
#define CCTBX_ARRAY_FAMILY_VERSA_PLAIN_H

#include <cctbx/array_family/shared_base.h>
#include <cctbx/array_family/versa_helpers.h>

namespace cctbx { namespace af {

  template <typename ElementType,
            typename AccessorType = grid_accessor<1> >
  class versa_plain : public shared_base<ElementType>
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef detail::char_block handle_type;

      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;

      CCTBX_ARRAY_FAMILY_VERSA_CONSTRUCTORS(versa_plain)

      versa_plain(const handle_type& handle, const accessor_type& ac)
        : shared_base<ElementType>(handle) {
        this->resize(ac);
      }

      versa_plain(const handle_type& handle, const size_type& sz)
        : shared_base<ElementType>(handle) {
        this->resize(accessor_type(sz));
      }

      void resize(const accessor_type& ac) {
        m_accessor = ac;
        shared_base<ElementType>(*this).resize(m_accessor.size1d());
      }

      const accessor_type& accessor() const { return m_accessor; }
      size_type size() const { return m_accessor.size1d(); }

      versa_plain<ElementType> as_1d() {
        return versa_plain<ElementType>(this->handle(), this->size());
      }

      CCTBX_ARRAY_FAMILY_TAKE_VERSA_REF(this->begin(), this->accessor())

            value_type& operator()(const index_type& i)       {
        return this->begin()[m_accessor(i)];
      }
      const value_type& operator()(const index_type& i) const {
        return this->begin()[m_accessor(i)];
      }

      // Convenience operator()

      const value_type& operator()(int i0) const {
        return operator()(index_type(i0));
      }
            value_type& operator()(int i0)       {
        return operator()(index_type(i0));
      }
      const value_type& operator()(int i0,
                                   int i1) const {
        return operator()(index_type(i0, i1));
      }
            value_type& operator()(int i0,
                                   int i1)       {
        return operator()(index_type(i0, i1));
      }
      const value_type& operator()(int i0,
                                   int i1,
                                   int i2) const {
        return operator()(index_type(i0, i1, i2));
      }
            value_type& operator()(int i0,
                                   int i1,
                                   int i2)       {
        return operator()(index_type(i0, i1, i2));
      }

    protected:
      accessor_type m_accessor;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_VERSA_PLAIN_H
