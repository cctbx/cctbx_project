#ifndef CCTBX_GEOMETRY_RESTRAINTS_SELECT_H
#define CCTBX_GEOMETRY_RESTRAINTS_SELECT_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/selections.h>

namespace cctbx { namespace geometry_restraints {

  namespace detail {

    template <typename ArrayType>
    struct new_i_seqs;

    // For angle, dihedral or chirality proxies.
    template <std::size_t N>
    struct new_i_seqs<af::tiny<unsigned, N> > : af::tiny<unsigned, N>
    {
      bool is_valid;

      new_i_seqs(
        af::tiny<unsigned, N> const& i_seqs,
        af::const_ref<std::size_t> const& reindexing_array)
      :
        is_valid(false)
      {
        std::size_t n_seq = reindexing_array.size();
        for(unsigned i=0;i<N;i++) {
          unsigned i_seq = i_seqs[i];
          CCTBX_ASSERT(i_seq < n_seq);
          this->elems[i] = static_cast<unsigned>(reindexing_array[i_seq]);
          if (this->elems[i] == n_seq) return;
        }
        is_valid = true;
      }
    };

    // For planarity proxies.
    template <>
    struct new_i_seqs<af::shared<std::size_t> > : af::shared<std::size_t>
    {
      bool is_valid;

      new_i_seqs(
        af::shared<std::size_t> const& i_seqs_,
        af::const_ref<std::size_t> const& reindexing_array)
      :
        is_valid(false)
      {
        af::const_ref<std::size_t> i_seqs = i_seqs_.const_ref();
        std::size_t n_seq = reindexing_array.size();
        for(std::size_t i=0;i<i_seqs.size();i++) {
          std::size_t i_seq = i_seqs[i];
          CCTBX_ASSERT(i_seq < n_seq);
          std::size_t new_i_seq = reindexing_array[i_seq];
          if (new_i_seq != n_seq) {
            this->push_back(new_i_seq);
          }
        }
        is_valid = (this->size() > 3);
      }
    };

  } // namespace detail

  /*! \brief Applies selection to array of angle, dihedral, chirality
      or planarity proxies.
   */
  template <typename ProxyType>
  af::shared<ProxyType>
  shared_proxy_select(
    af::const_ref<ProxyType> const& self,
    std::size_t n_seq,
    af::const_ref<std::size_t> const& iselection)
  {
    af::shared<ProxyType> result;
    af::shared<std::size_t>
      reindexing_array = scitbx::af::reindexing_array(n_seq, iselection);
    af::const_ref<std::size_t>
      reindexing_array_ref = reindexing_array.const_ref();
    for(std::size_t i_proxy=0;i_proxy<self.size();i_proxy++) {
      ProxyType const& p = self[i_proxy];
      detail::new_i_seqs<typename ProxyType::i_seqs_type>
        new_i_seqs(p.i_seqs, reindexing_array_ref);
      if (new_i_seqs.is_valid) {
        result.push_back(ProxyType(new_i_seqs, p));
      }
    }
    return result;
  }

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_SELECT_H
