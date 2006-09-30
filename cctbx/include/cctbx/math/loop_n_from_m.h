#ifndef CCTBX_MATH_LOOP_N_FROM_M_H
#define CCTBX_MATH_LOOP_N_FROM_M_H

#include <algorithm>
#include <cstddef>
#include <cctbx/error.h>

namespace cctbx { namespace math {

  template <std::size_t MaxN>
  class loop_n_from_m
  {
    public:
      loop_n_from_m()
      : m_(0), n_(0), over_(1)
      {}

      loop_n_from_m(std::size_t m, std::size_t n)
      : m_(m), n_(n), over_(0)
      {
        CCTBX_ASSERT(m_ >= n_);
        CCTBX_ASSERT(MaxN >= n_);
        for(std::size_t i=0;i<n_;i++) current_[i] = i;
      }

      bool incr()
      {
        if (n_ > 0) {
          std::size_t p, l;
          p = l = n_ - 1;
          for (;;) {
                current_[p]++;
            if (current_[p] + l == m_ + p) {
              if (p == 0) break;
              p--;
            }
            else if (p < l) {
              current_[p + 1] = current_[p];
                        p++;
            }
            else {
              return true;
            }
          }
        }
        over_++;
        return false;
      }

      std::size_t m() const { return m_; }

      std::size_t n() const { return n_; }

      std::size_t operator[](std::size_t i) const { return current_[i]; }

      std::size_t over() const { return over_; }

    protected:
      std::size_t m_;
      std::size_t n_;
      std::size_t current_[MaxN];
      std::size_t over_;
  };

}} // namespace cctbx::math

#endif // CCTBX_MATH_LOOP_N_FROM_M_H
