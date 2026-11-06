#ifndef CCTBX_GEOMETRY_RESTRAINTS_N_PROXIES_EXCLUDING_H
#define CCTBX_GEOMETRY_RESTRAINTS_N_PROXIES_EXCLUDING_H

#include <cctbx/geometry_restraints/angle.h>
#include <cctbx/geometry_restraints/dihedral.h>
#include <cctbx/geometry_restraints/chirality.h>
#include <cctbx/geometry_restraints/planarity.h>
#include <cctbx/geometry_restraints/parallelity.h>
#include <cctbx/geometry_restraints/bond_similarity.h>

namespace cctbx { namespace geometry_restraints {

  //! True if none of i_seqs is flagged in exclude.
  template <typename IndexArrayType>
  inline bool
  proxy_is_not_excluded(
    IndexArrayType const& i_seqs,
    af::const_ref<bool> const& exclude)
  {
    for(std::size_t i=0;i<i_seqs.size();i++) {
      CCTBX_ASSERT(i_seqs[i] < exclude.size());
      if (exclude[i_seqs[i]]) return false;
    }
    return true;
  }

  /*! \brief Number of proxies none of whose atoms is flagged in exclude.

      Used for the restraint count that normalizes the geometry target.
   */
  inline std::size_t
  n_proxies_excluding(
    af::const_ref<angle_proxy> const& proxies,
    af::const_ref<bool> const& exclude)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      if (proxy_is_not_excluded(proxies[i].i_seqs, exclude)) result++;
    }
    return result;
  }

  inline std::size_t
  n_proxies_excluding(
    af::const_ref<dihedral_proxy> const& proxies,
    af::const_ref<bool> const& exclude)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      if (proxy_is_not_excluded(proxies[i].i_seqs, exclude)) result++;
    }
    return result;
  }

  inline std::size_t
  n_proxies_excluding(
    af::const_ref<chirality_proxy> const& proxies,
    af::const_ref<bool> const& exclude)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      if (proxy_is_not_excluded(proxies[i].i_seqs, exclude)) result++;
    }
    return result;
  }

  inline std::size_t
  n_proxies_excluding(
    af::const_ref<planarity_proxy> const& proxies,
    af::const_ref<bool> const& exclude)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      if (proxy_is_not_excluded(proxies[i].i_seqs.const_ref(), exclude)) {
        result++;
      }
    }
    return result;
  }

  inline std::size_t
  n_proxies_excluding(
    af::const_ref<parallelity_proxy> const& proxies,
    af::const_ref<bool> const& exclude)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      if (   proxy_is_not_excluded(proxies[i].i_seqs.const_ref(), exclude)
          && proxy_is_not_excluded(proxies[i].j_seqs.const_ref(), exclude)) {
        result++;
      }
    }
    return result;
  }

  inline std::size_t
  n_proxies_excluding(
    af::const_ref<bond_similarity_proxy> const& proxies,
    af::const_ref<bool> const& exclude)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      af::const_ref<af::tiny<std::size_t, 2> > pairs
        = proxies[i].i_seqs.const_ref();
      bool keep = true;
      for(std::size_t j=0;j<pairs.size();j++) {
        if (!proxy_is_not_excluded(pairs[j], exclude)) {
          keep = false;
          break;
        }
      }
      if (keep) result++;
    }
    return result;
  }

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_N_PROXIES_EXCLUDING_H
