/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 Jul: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

/*! \file
    Toolbox for Miller indices.
 */

#ifndef CCTBX_MILLER_H
#define CCTBX_MILLER_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/misc_functions.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx {
  //! Miller index namespace.
  namespace miller {

  //! Miller index class.
  template <typename NumType = int>
  class index : public scitbx::vec3<NumType>
  {
    public:
      typedef scitbx::vec3<NumType> base_type;

      //! Default constructor: (h,k,l) = (0,0,0)
      index() : base_type(0,0,0) {}

      //! Construction with an instance of the base type.
      index(base_type const& h) : base_type(h) {}

      //! Construction from the highest array type in the inheritance tree.
      template <typename OtherNumType>
      index(af::tiny_plain<OtherNumType, 3> const& v)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = v[i];
      }

      //! Construction with raw pointer to index elements.
      template <typename OtherNumType>
      explicit
      index(const OtherNumType* hkl)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = hkl[i];
      }

      //! Construction from individual h,k,l.
      index(NumType const& h, NumType const& k, NumType const& l)
      : base_type(h, k, l)
      {}

      //! Definition of sort order for human-readable listings.
      /*! This comparison is computationally more expensive than
          hash_compare below.
       */
      bool operator<(index const& other) const
      {
        using scitbx::fn::absolute;
        const int P[3] = {2, 0, 1};
        for(std::size_t i=0;i<3;i++) {
          if (this->elems[P[i]] >= 0 && other[P[i]] <  0) return true;
          if (this->elems[P[i]] <  0 && other[P[i]] >= 0) return false;
        }
        for(std::size_t i=0;i<3;i++) {
          if (  absolute(this->elems[P[i]])
              < absolute(other[P[i]])) return true;
          if (  absolute(this->elems[P[i]])
              > absolute(other[P[i]])) return false;
        }
        return false;
      }

      //! Test this > other, based on sort order defined by operator<().
      bool operator>(index const& other) const
      {
        if (*this < other) return false;
        for(std::size_t i=0;i<3;i++) {
          if (this->elems[i] != other[i]) return true;
        }
        return false;
      }

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300 // VC++ 7.0
      // work around compiler bug
      index operator-() const
      {
        return index(-static_cast<base_type>(*this));
      }
#endif
  };

  /*! \brief Definition of fast comparison for use in,
      e.g., std::map<miller::index<> >.
   */
  template <typename NumType = int>
  class hash_compare
  {
    public:
      //! This fast comparison function is implemented as operator().
      bool operator()(index<NumType> const& h1, index<NumType> const& h2) const
      {
        for(std::size_t i=0;i<3;i++) {
          if (h1[i] < h2[i]) return true;
          if (h1[i] > h2[i]) return false;
        }
        return false;
      }
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_H
