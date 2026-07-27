// (c) O.V.D., OlexSys Ltd, 2025
#pragma once
#include <cctbx/xray/scatterer.h>

namespace cctbx {
namespace xray {

  template <typename FloatType, class crd_t, class mask_info, uint64_t cell_m>
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

      uint64_t get_id(int z, const fractional<FloatType>& site,
        short data = 0, FloatType multiplier = 1) const
      {
        return scatterer_id_base<FloatType, crd_t, mask_info, cell_m>(z, site, data).id;
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

      size_t index_of_fractional(const fractional<FloatType>& site, int Z, int sdata = 0, FloatType eps = 1e-3) const
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

    template <typename FloatType>
    struct scatterer_ID_lookup {
      typedef scatterer<FloatType> scatterer_t;
      const uctbx::unit_cell& u_cell;
      typedef scitbx::vec3<FloatType> cart_t;
      typedef std::pair<uint64_t,uint64_t> scatterer_id_t;
      af::shared<scatterer_t> scatterers;
      std::vector<cart_t> crds;
      std::vector<scatterer_id_t> scatterer_as_ids;
      af::shared<int> data;
      

      scatterer_ID_lookup(const uctbx::unit_cell& u_cell)
        : u_cell(u_cell)
      {}

      scatterer_ID_lookup(const uctbx::unit_cell& u_cell,
          const af::shared<scatterer_t>& scatterers,
          const af::shared<int>& data = af::shared<int>())
        : u_cell(u_cell)
      {
        this->data = data;
        this->scatterers = scatterers;
        crds.reserve(scatterers.size());
        scatterer_as_ids.reserve(scatterers.size());
        for (size_t i = 0; i < scatterers.size(); i++) {
          crds.push_back(u_cell.orthogonalize(scatterers[i].site));
          scatterer_as_ids.push_back(scatterers[i].get_id_big(data[i]).as_uint64());
        }
      }

      size_t index_from_id(const scatterer_id_big<FloatType, fractional<FloatType> >& sc_id) const
      {
        scatterer_id_t sc_id_pair = sc_id.as_uint64();
        tmp_logs << "TEST ID: " << sc_id_pair.first << ", " << sc_id_pair.second << std::endl;
        for (size_t i = 0; i < scatterers.size(); i++) {
          tmp_logs << "CHECK ID: " << scatterer_as_ids[i].first << ", " << scatterer_as_ids[i].second << std::endl;
          if (scatterer_as_ids[i] == sc_id_pair) {
            return i;
          }
        }
        return ~0;
      }

      size_t index_cartesian(const scatterer_id_big<FloatType, fractional<FloatType>>& sc_test, FloatType eps = 1e-3) const
      {
        cart_t crd_test = u_cell.orthogonalize(sc_test.get_crd());
        uint8_t Z_test = sc_test.get_z();
        int16_t data_test = sc_test.get_data();

        for (size_t i = 0; i < scatterers.size(); i++) {
          if (scatterers[i].get_atomic_number() == Z_test && data[i] == data_test) {
            FloatType qd = (crds[i] - crd_test).length();
            if (qd < eps) {
              return i;
            }
          }
        }
        return ~0;
      }

      size_t get_index(const scatterer_id_big<FloatType, fractional<FloatType> >& sc_id, FloatType eps = 1e-3) const
      {
        size_t idx = index_from_id(sc_id);
        if (idx != ~0) {
          return idx;
        }
        
        idx = index_cartesian(sc_id, eps);
        if (idx != ~0) {
          return idx;
        }
        throw CCTBX_ERROR("Could not locate scatterer");
      }
    };
}} // namespace cctbx::xray
