#pragma once
#include <cctbx/coordinates.h>
namespace cctbx {
  namespace xray {
    /* Stores coordinates with float type precision, adopted from Olex2 */
    template <typename FloatType>
    struct scatterer_id {
    const static uint64_t
      z_mask = 0x00000000000000FF,
      a_mask = 0x0000000001FFFF00,
      b_mask = 0x000003FFFE000000,
      c_mask = 0x07FFFC0000000000,
      a_sig = 0x0800000000000000,
      b_sig = 0x1000000000000000,
      c_sig = 0x2000000000000000,
      mask_m = 0x000000000001FFFF, // max crd value
      cell_m = 16 // max unit cells in each direction
      ;
    /*  0-8 - z
    9-25, 26-42, 43-59 - a, b c, 60-62 - signs
    */
    uint64_t id;

    scatterer_id(const scatterer_id& sid)
      : id(sid.id)
    {}

    scatterer_id(uint64_t id = ~0)
      : id(id)
    {}

    scatterer_id(int8_t z, const fractional<FloatType>& crd, FloatType multiplier = 1) {
      id = ((uint64_t)z) & z_mask;
      static const int64_t k = mask_m / cell_m;
      int64_t x = multiplier == 1 ? (int64_t)(crd[0] * k)
        : ((int64_t)(crd[0] * multiplier)) / multiplier * k;
      if (x < 0) {
        id |= a_sig;
        id |= (((-x) << 8) & a_mask);
      }
      else {
        id |= ((std::abs(x) << 8) & a_mask);
      }
      x = multiplier == 1 ? (int64_t)(crd[1] * k)
        : ((int64_t)(crd[1] * multiplier)) / multiplier * k;
      if (x < 0) {
        id |= b_sig;
        id |= (((-x) << 25) & b_mask);
      }
      else {
        id |= ((x << 25) & b_mask);
      }
      x = multiplier == 1 ? (int64_t)(crd[2] * k)
        : ((int64_t)(crd[2] * multiplier)) / multiplier * k;
      if (x < 0) {
        id |= c_sig;
        id |= (((-x) << 42) & c_mask);
      }
      else {
        id |= ((x << 42) & c_mask);
      }
    }

    scatterer_id& operator = (const scatterer_id& i) {
      this->id = i.id;
      return *this;
    }

    scatterer_id& operator = (const uint64_t& i) {
      this->id = i;
      return *this;
    }

    bool operator == (const scatterer_id& i) const {
      return id == i.id;
    }

    int Compare(const scatterer_id& i) const {
      return id < i.id ? -1 : (id > i.id ? 1 : 0);
    }

    fractional<FloatType> get_crd() const {
      static const double k = (double)cell_m / mask_m;
      fractional<FloatType> r(
        static_cast<FloatType>((id & a_mask) >> 8) * k,
        static_cast<FloatType>((id & b_mask) >> 25) * k,
        static_cast<FloatType>((id & c_mask) >> 42) * k);
      if ((id & a_sig) != 0) {
        r[0] = -r[0];
      }
      if ((id & b_sig) != 0) {
        r[1] = -r[1];
      }
      if ((id & c_sig) != 0) {
        r[2] = -r[2];
      }
      return r;
    }

    int8_t get_z() const {
      return (int8_t)(id & z_mask);
    }
  };

}} // namespace cctbx::xray