#ifndef CCTBX_GEOMETRY_RESTRAINTS_PROXY_SELECT_H
#define CCTBX_GEOMETRY_RESTRAINTS_PROXY_SELECT_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/selections.h>

namespace cctbx { namespace geometry_restraints {

  /*! \brief Applies selection on array of atoms to array of
      angle, dihedral or chirality proxies.
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
      typename ProxyType::i_seqs_type new_i_seqs;
      bool is_valid = true;
      for(unsigned i=0;i<p.i_seqs.size();i++) {
        unsigned i_seq = p.i_seqs[i];
        CCTBX_ASSERT(i_seq < n_seq);
        new_i_seqs[i] = static_cast<unsigned>(reindexing_array_ref[i_seq]);
        if (new_i_seqs[i] == n_seq) {
          is_valid = false;
          break;
        }
      }
      if (is_valid) {
        result.push_back(ProxyType(new_i_seqs, p));
      }
    }
    return result;
  }

  template <typename ProxyType>
  af::shared<ProxyType>
  shared_proxy_select_origin(
    af::const_ref<ProxyType> const& self,
    unsigned char origin_id)
  {
    af::shared<ProxyType> result;
    for(std::size_t i_proxy=0; i_proxy<self.size(); i_proxy++) {
      ProxyType const& p = self[i_proxy];
      if (p.origin_id == origin_id) {
        result.push_back(ProxyType(p));
      }
    }
    return result;
  }

  /*! \brief Applies selection on array of atoms to array of
      planarity proxies.
   */
  template <typename ProxyType>
  af::shared<ProxyType>
  shared_planarity_proxy_select(
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
      af::const_ref<typename ProxyType::i_seqs_type::value_type>
        i_seqs = p.i_seqs.const_ref();
      af::const_ref<double>
        weights = p.weights.const_ref();
      af::shared<typename ProxyType::i_seqs_type::value_type> new_i_seqs;
      af::shared<double> new_weights;
      for(std::size_t i=0;i<i_seqs.size();i++) {
        std::size_t i_seq = i_seqs[i];
        CCTBX_ASSERT(i_seq < n_seq);
        std::size_t new_i_seq = reindexing_array_ref[i_seq];
        if (new_i_seq != n_seq) {
          new_i_seqs.push_back(new_i_seq);
          new_weights.push_back(weights[i]);
        }
      }
      if (new_i_seqs.size() > 3) {
        result.push_back(ProxyType(new_i_seqs, new_weights, p.origin_id));
      }
    }
    return result;
  }

  template <typename ProxyType>
  af::shared<ProxyType>
  shared_parallelity_proxy_select(
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
      af::const_ref<typename ProxyType::i_seqs_type::value_type>
          i_seqs = p.i_seqs.const_ref();
      af::const_ref<typename ProxyType::i_seqs_type::value_type>
          j_seqs = p.j_seqs.const_ref();
      af::shared<typename ProxyType::i_seqs_type::value_type> new_i_seqs;
      af::shared<typename ProxyType::i_seqs_type::value_type> new_j_seqs;
      for(std::size_t i=0;i<i_seqs.size();i++) {
        std::size_t i_seq = i_seqs[i];
        CCTBX_ASSERT(i_seq < n_seq);
        std::size_t new_i_seq = reindexing_array_ref[i_seq];
        if (new_i_seq != n_seq) {
          new_i_seqs.push_back(new_i_seq);
        }
      }
      for(std::size_t j=0;j<j_seqs.size();j++) {
        std::size_t j_seq = j_seqs[j];
        CCTBX_ASSERT(j_seq < n_seq);
        std::size_t new_j_seq = reindexing_array_ref[j_seq];
        if (new_j_seq != n_seq) {
          new_j_seqs.push_back(new_j_seq);
        }
      }
      if (new_i_seqs.size() > 2 && new_j_seqs.size() > 2) {
        result.push_back(ProxyType(new_i_seqs, new_j_seqs, p));
      }
    }
    return result;
  }



  template <typename ProxyType>
  af::shared<ProxyType>
  shared_proxy_remove(
    af::const_ref<ProxyType> const& self,
    af::const_ref<bool> const& selection)
  {
    af::shared<ProxyType> result;
    for(std::size_t i_proxy=0;i_proxy<self.size();i_proxy++) {
      ProxyType const& p = self[i_proxy];
      af::const_ref<typename ProxyType::i_seqs_type::value_type>
        i_seqs = p.i_seqs.const_ref();
      for(unsigned i=0;i<i_seqs.size();i++) {
        unsigned i_seq = i_seqs[i];
        CCTBX_ASSERT(i_seq < selection.size());
        if (!selection[i_seq]) {
          result.push_back(p);
          break;
        }
      }
    }
    return result;
  }

  template <typename ProxyType>
  af::shared<ProxyType>
  shared_proxy_remove(
    af::const_ref<ProxyType> const& self,
    af::const_ref<size_t> const& selection)
  {
    af::shared<ProxyType> result;
    for(std::size_t i_proxy=0;i_proxy<self.size();i_proxy++) {
      ProxyType const& p = self[i_proxy];
      af::const_ref<typename ProxyType::i_seqs_type::value_type>
        i_seqs = p.i_seqs.const_ref();
      int match_ctr = 0;
      for(unsigned i=0;i<i_seqs.size();i++) {
        unsigned i_seq = i_seqs[i];
        for(unsigned j=0;j<selection.size();j++) {
          if (i_seq == selection[j]) {
            match_ctr++;
          }
        }
      }
      if (match_ctr != i_seqs.size()) {
        result.push_back(p);
      }
    }
    return result;
  }

  template <typename ProxyType>
  af::shared<ProxyType>
  shared_proxy_remove(
    af::const_ref<ProxyType> const& self,
    unsigned char origin_id)
  {
    af::shared<ProxyType> result;
    for(std::size_t i_proxy=0; i_proxy<self.size(); i_proxy++) {
      ProxyType const& p = self[i_proxy];
      if (p.origin_id != origin_id) {
        result.push_back(p);
      }
    }
    return result;
  }

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_PROXY_SELECT_H
