#ifndef SCITBX_ARRAY_FAMILY_COUNTS_H
#define SCITBX_ARRAY_FAMILY_COUNTS_H

#include <scitbx/array_family/ref.h>
#include <boost/shared_ptr.hpp>
#include <exception>

namespace scitbx { namespace af {

  template <typename IntType, typename MapType>
  struct counts
  {
    static
    boost::shared_ptr<MapType>
    unlimited(af::const_ref<IntType> const& self)
    {
      boost::shared_ptr<MapType> result(new MapType);
      MapType& m = *result;
      for(std::size_t i=0;i<self.size();i++) {
        m[static_cast<typename MapType::key_type>(self[i])]++;
      }
      return result;
    }

    static
    boost::shared_ptr<MapType>
    limited(af::const_ref<IntType> const& self, std::size_t max_keys)
    {
      boost::shared_ptr<MapType> result(new MapType);
      MapType& m = *result;
      for(std::size_t i=0;i<self.size();i++) {
        m[static_cast<typename MapType::key_type>(self[i])]++;
        if (m.size() > max_keys) {
          throw std::runtime_error(
            "scitbx::af::counts::limited: max_keys exceeded.");
        }
      }
      return result;
    }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_COUNTS_H
