/*! \file
    Handling of Miller indices.
 */

#ifndef CCTBX_MILLER_H
#define CCTBX_MILLER_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/misc_functions.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cstdio>

namespace cctbx {
  //! Miller index namespace.
  namespace miller {

  //! Miller index class.
  template <typename NumType=int>
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
          fast_less_than below.
       */
      bool operator<(index const& other) const
      {
        const int P[3] = {2, 0, 1};
        for(std::size_t i=0;i<3;i++) {
          if (this->elems[P[i]] >= 0 && other[P[i]] <  0) return true;
          if (this->elems[P[i]] <  0 && other[P[i]] >= 0) return false;
        }
        for(std::size_t i=0;i<3;i++) {
          if (  scitbx::fn::absolute(this->elems[P[i]])
              < scitbx::fn::absolute(other[P[i]])) return true;
          if (  scitbx::fn::absolute(this->elems[P[i]])
              > scitbx::fn::absolute(other[P[i]])) return false;
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

#if defined(__GNUC__) \
 && __GNUC__ == 3 \
 && __GNUC_MINOR__ == 2 \
 && __GNUC_PATCHLEVEL__ != 0 // 3 known to fail, 0 OK, others unknown
      // work around compiler bug
      bool operator!=(index const& other) const
      {
        return static_cast<base_type>(*this) != static_cast<base_type>(other);
      }
#endif

      std::string
      as_string() const
      {
        char buf[128];
        buf[127] = '\0';
        std::snprintf(buf, sizeof(buf), "(%ld,%ld,%ld)",
          static_cast<long>(this->elems[0]),
          static_cast<long>(this->elems[1]),
          static_cast<long>(this->elems[2]));
        CCTBX_ASSERT(buf[127] == '\0');
        return std::string(buf);
      }
  };

  /*! \brief Definition of fast comparison for use in,
      e.g., std::map<miller::index<> >.
   */
  template <typename NumType=int>
  class fast_less_than
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
