#pragma once
#include <cctbx/coordinates.h>
namespace cctbx {
  namespace xray {
    /* Stores coordinates with float type precision, adopted from Olex2. */
    
    template <typename FloatType, class mask_info, uint64_t cell_m> 
    struct scatterer_id_base : public mask_info {
      uint64_t id;

      scatterer_id_base(const scatterer_id_base& sid)
        : id(sid.id)
      {}

      scatterer_id_base(uint64_t id = ~0)
        : id(id)
      {}
      // use multiplier like 1e3 to 'crop' values
      scatterer_id_base(int8_t z, const fractional<FloatType>& crd, short data,
        FloatType multiplier = 1)
      {
        id = ((uint64_t)z) & mask_info::z_mask;
        static const int64_t k = mask_info::mask_m / cell_m;
        int64_t x = multiplier == 1 ? (int64_t)(crd[0] * k)
          : ((int64_t)(crd[0] * multiplier)) / multiplier * k;
        if (x < 0) {
          id |= mask_info::a_sig;
          id |= (((-x) << 8) & mask_info::a_mask);
        }
        else {
          id |= ((std::abs(x) << 8) & mask_info::a_mask);
        }
        x = multiplier == 1 ? (int64_t)(crd[1] * k)
          : ((int64_t)(crd[1] * multiplier)) / multiplier * k;
        if (x < 0) {
          id |= mask_info::b_sig;
          id |= (((-x) << (8 + mask_info::a_shift)) & mask_info::b_mask);
        }
        else {
          id |= ((x << (8 + mask_info::a_shift)) & mask_info::b_mask);
        }
        x = multiplier == 1 ? (int64_t)(crd[2] * k)
          : ((int64_t)(crd[2] * multiplier)) / multiplier * k;
        if (x < 0) {
          id |= mask_info::c_sig;
          id |= (((-x) << (8 + mask_info::a_shift * 2)) & mask_info::c_mask);
        }
        else {
          id |= ((x << (8 + mask_info::a_shift * 2)) & mask_info::c_mask);
        }
        x = data;
        id |= (((int64_t)data << (8 + mask_info::a_shift * 3) + 3) & mask_info::d_mask);
      }

      bool test(const fractional<FloatType>& crd) const {
        return std::abs(crd[0]) <= cell_m &&
          std::abs(crd[1]) <= cell_m &&
          std::abs(crd[2]) <= cell_m;
      }

      scatterer_id_base& operator = (const scatterer_id_base& i) {
        this->id = i.id;
        return *this;
      }

      scatterer_id_base& operator = (const uint64_t& i) {
        this->id = i;
        return *this;
      }

      bool operator == (const scatterer_id_base& i) const {
        return id == i.id;
      }

      bool operator != (const scatterer_id_base& i) const {
        return id != i.id;
      }

      int Compare(const scatterer_id_base& i) const {
        return id < i.id ? -1 : (id > i.id ? 1 : 0);
      }

      fractional<FloatType> get_crd() const {
        static const double k = (double)cell_m / mask_info::mask_m;
        fractional<FloatType> r(
          static_cast<FloatType>((id & mask_info::a_mask) >> 8) * k,
          static_cast<FloatType>((id & mask_info::b_mask) >> (8 + mask_info::a_shift)) * k,
          static_cast<FloatType>((id & mask_info::c_mask) >> (8 + mask_info::a_shift * 2)) * k);
        if ((id & mask_info::a_sig) != 0) {
          r[0] = -r[0];
        }
        if ((id & mask_info::b_sig) != 0) {
          r[1] = -r[1];
        }
        if ((id & mask_info::c_sig) != 0) {
          r[2] = -r[2];
        }
        return r;
      }

      int8_t get_z() const {
        return (int8_t)(id & mask_info::z_mask);
      }

      // this wil always be returned as unsigned value
      short get_data() const {
        return (short)((id & mask_info::d_mask) >> (8 + mask_info::a_shift * 3 + 3));
      }
    };

    /*  0-8 - z, 8-25, 25-42, 42-59 - a, b c, 59-61 - signs, 62-64 - data
    precision 16: ~0.0025, 1: ~0.0000077
    */
    struct scatterer_id_masks_d2 {
      const static uint64_t
        z_mask = 0x00000000000000FF,
        a_mask = 0x0000000001FFFF00,
        b_mask = 0x000003FFFE000000,
        c_mask = 0x07FFFC0000000000,
        a_sig  = 0x0800000000000000,
        b_sig  = 0x1000000000000000,
        c_sig  = 0x2000000000000000,
        d_mask = 0xC000000000000000,
        mask_m = 0x000000000001FFFF; // max crd value
      const static int a_shift = 17;
    };

    /* 2 extra bits for data */
    template <typename FloatType, uint64_t cell_m>
    struct scatterer_id_2 : public scatterer_id_base<FloatType, scatterer_id_masks_d2, cell_m> {
      typedef scatterer_id_base<FloatType, scatterer_id_masks_d2, cell_m> parent_t;
      scatterer_id_2(const scatterer_id_2& sid)
        : parent_t(sid)
      {}

      scatterer_id_2(uint64_t id = ~0)
        : parent_t(id)
      {}

      scatterer_id_2(int8_t z, const fractional<FloatType>& crd, short data,
        FloatType multiplier = 1)
        : parent_t(z, crd, data, multiplier)
      {}
    };

    /*  0-8 - z, 8-24, 24-40, 40-56 - a, b c, 56-59 - signs, 59-64 - data
    precision 16: ~0.004, 1: ~0.00002
    */
    struct scatterer_id_masks_d5 {
      const static uint64_t
        z_mask = 0x00000000000000FF,
        a_mask = 0x0000000000FFFF00,
        b_mask = 0x000000FFFF000000,
        c_mask = 0x00FFFF0000000000,
        a_sig =  0x0100000000000000,
        b_sig =  0x0200000000000000,
        c_sig =  0x0400000000000000,
        d_mask = 0xF800000000000000,
        mask_m = 0x000000000000FFFF; // max crd value
      const static int a_shift = 16;
    };

    /* 5 extra bits for data */
    template <typename FloatType, uint64_t cell_m>
    struct scatterer_id_5 : public scatterer_id_base<FloatType, scatterer_id_masks_d5, cell_m> {
      typedef scatterer_id_base<FloatType, scatterer_id_masks_d5, cell_m> parent_t;
      scatterer_id_5(const scatterer_id_5& sid)
        : parent_t(sid)
      {}

      scatterer_id_5(uint64_t id = ~0)
        : parent_t(id)
      {}
      // use multiplier like 1e3 to 'crop' values
      scatterer_id_5(int8_t z, const fractional<FloatType>& crd, short data,
        FloatType multiplier = 1)
        : parent_t(z, crd, data, multiplier)
      {}
    };

}} // namespace cctbx::xray