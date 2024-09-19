#pragma once
#include <cctbx/xray/scatterer.h>

namespace cctbx { 
namespace xray {

  template <typename FloatType, class mask_info, uint64_t cell_m>
  struct scatterer_lookup {
    typedef scatterer<> scatterer_t;

    std::map<uint64_t, const scatterer_t*> map;
    scatterer_lookup() {}

    scatterer_lookup(const af::shared<scatterer_t>& scatterers,
      FloatType multiplier = 1)
    {
      init(scatterers, multiplier);
    }
    
    scatterer_lookup(const af::shared<scatterer_t>& scatterers,
      const af::shared<int>& data, FloatType multiplier = 1)
    {
      init_with_data(scatterers, data, multiplier);
    }

    void init(const af::shared<scatterer_t>& scatterers,
        FloatType multiplier = 1)
      {
        for (size_t i = 0; i < scatterers.size(); i++) {
          map.insert(
            std::make_pair(
              scatterers[i].get_id<mask_info, cell_m>(0, multiplier).id,
              &scatterers[i]));
        }
      }

      void init_with_data(const af::shared<scatterer_t>& scatterers,
        const af::shared<int>& data, FloatType multiplier=1)
      {
        CCTBX_ASSERT(data.size() == 0 || scatterers.size() == data.size());
        for (size_t i = 0; i < scatterers.size(); i++) {
          map.insert(
            std::make_pair(
              scatterers[i].get_id<mask_info, cell_m>(
                data.size() == 0 ? 0 : data[i], multiplier).id,
              &scatterers[i]));
        }
      }

      const scatterer_t& find(uint64_t id) const {
        typename std::map<uint64_t, const scatterer_t*>::const_iterator si =
          map.find(id);
        CCTBX_ASSERT(si != map.end());
        return *si->second;
      }

      uint64_t get_id(int z, const fractional<>& site,
        short data = 0, FloatType multiplier = 1) const
      {
        return scatterer_id_base<FloatType, mask_info, cell_m>(z, site, data).id;
      }
    };

    template <typename FloatType>
    struct scatterer_cart_lookup {
      typedef scatterer<FloatType> scatterer_t;
      const uctbx::unit_cell& u_cell;
      typedef scitbx::vec3<FloatType> cart_t;
      typedef std::pair<FloatType, size_t> pair_t;
      af::shared<scatterer_t> scatterers;
      af::shared<int> data;
      std::vector<cart_t> crds;

      std::vector<pair_t> map;
      scatterer_cart_lookup(const uctbx::unit_cell& u_cell)
        : u_cell(u_cell)
      {}

        scatterer_cart_lookup(const uctbx::unit_cell& u_cell,
          const af::shared<scatterer_t>& scatterers)
        : u_cell(u_cell)
      {
        init(scatterers);
      }
      
      scatterer_cart_lookup(const uctbx::unit_cell& u_cell,
        const af::shared<scatterer_t>& scatterers,
        const af::shared<int>& data)
        : u_cell(u_cell)
      {
        init_with_data(scatterers, data);
      }
      
      void init(const af::shared<scatterer_t>& scatterers) {
        this->scatterers = scatterers;
        crds.reserve(scatterers.size());
        map.reserve(scatterers.size());
        for (size_t i = 0; i < scatterers.size(); i++) {
          crds.push_back(u_cell.orthogonalize(scatterers[i].site));
          map.push_back(std::make_pair(crds[i].length(), i));
        }
        std::sort(map.begin(), map.end(), &sort_qd);
      }

      void init_with_data(const af::shared<scatterer_t>& scatterers,
        const af::shared<int>& data)
      {
        this->scatterers = scatterers;
        this->data = data;
        crds.reserve(scatterers.size());
        map.reserve(scatterers.size());
        for (size_t i = 0; i < scatterers.size(); i++) {
          crds.push_back(u_cell.orthogonalize(scatterers[i].site));
          map.push_back(std::make_pair(crds[i].length(), i));
        }
        std::sort(map.begin(), map.end(), &sort_qd);
      }

      const scatterer_t& find_fractional(const fractional<FloatType>& site,
        int Z, int sdata = 0, FloatType eps = 1e-3) const
      {
        return find_cartesian(u_cell.orthogonalize(site), Z, sdata, eps);
      }

      const scatterer_t& find_cartesian(const cart_t& crd, int Z, int sdata = 0,
        FloatType eps = 1e-3) const
      {
        size_t idx = index_of_cartesian(crd, Z, sdata, eps);
        if (idx == ~0) {
          throw CCTBX_ERROR("Could not locate scatterer at the given coordinate");
        }
        return scatterers[idx];
      }

      size_t index_of_fractional(const fractional<FloatType>& site,
        int Z, int sdata = 0, FloatType eps = 1e-3) const
      {
        return index_of_cartesian(u_cell.orthogonalize(site), Z, sdata, eps);
      }

      size_t index_of_cartesian(const cart_t& crd, int Z, int sdata = 0,
        FloatType eps = 1e-3) const
      {
        typedef typename std::vector<pair_t>::const_iterator itr_t;
        FloatType sql = crd.length();
        FloatType eps_qd = sql * eps;
        itr_t itr = std::upper_bound(
          map.begin(), map.end(),
          std::make_pair(sql, size_t(~0)), less_qd());
        CCTBX_ASSERT(itr != map.end());
        itr_t itr1 = itr;
        FloatType diff = std::abs((*itr).first - sql);
        while (diff < eps_qd && itr != map.end()) {
          size_t idx = (*itr).second;
          const scatterer_t& s = scatterers[idx];
          FloatType qd = (crd - crds[idx]).length();
          if (s.get_atomic_number() == Z && qd < eps) {
            if (data.size() == 0 || data[idx] == sdata) {
              return idx;
            }
          }
          itr++;
          if (itr != map.end()) {
            diff = std::abs((*itr).first - sql);
          }
        }
        itr = itr1;
        CCTBX_ASSERT(itr != map.begin());
        itr--;
        diff = std::abs((*itr).first - sql);
        while (diff < eps_qd) {
          size_t idx = (*itr).second;
          const scatterer_t& s = scatterers[idx];
          FloatType qd = (crd - crds[idx]).length();
          if (s.get_atomic_number() == Z && qd < eps) {
            if (data.size() == 0 || data[idx] == sdata) {
              return idx;
            }
          }
          if (itr == map.begin()) {
            break;
          }
          itr--;
          diff = std::abs((*itr).first - sql);
        }
        return ~0;
      }

      static bool sort_qd(const pair_t& a, const pair_t& b) {
        return a.first < b.first;
      }

      struct less_qd {
        bool operator()(const pair_t& a, const pair_t& b) const {
          return sort_qd(a, b);
        }
      };
    };
}} // namespace cctbx::xray