#ifndef SCITBX_MATH_UNIMODULAR_GENERATOR_H
#define SCITBX_MATH_UNIMODULAR_GENERATOR_H

#include <scitbx/mat3.h>

namespace scitbx { namespace math {

  template <typename IntType>
  class unimodular_generator
  {
    private:
      IntType range_;
      bool at_end_;
      unsigned return_target_id_;
      IntType i0, i1, i2, i3, i4, i5, i6, i7, i8;
      IntType i37, i38, i48, i4857, i3746, i3856, i04857;
      IntType i04857_one, one_i04857_13856;

    public:
      unimodular_generator() {}

      unimodular_generator(IntType const& range)
      :
        range_(range),
        at_end_(false),
        return_target_id_(0)
      {
        incr();
      }

      bool
      at_end() const { return at_end_; }

      scitbx::mat3<IntType>
      next()
      {
        SCITBX_ASSERT(!at_end_);
        scitbx::mat3<IntType> result(i0,i1,i2,i3,i4,i5,i6,i7,i8);
        incr();
        return result;
      }

      std::size_t
      count()
      {
        std::size_t result = 0;
        while (!at_end_) {
          result++;
          incr();
        }
        return result;
      }

    private:
      void
      incr()
      {
        if (return_target_id_ == 1) goto return_target_1;
        if (return_target_id_ == 2) goto return_target_2;
        if (return_target_id_ == 3) goto return_target_3;
        // 0 1 2
        // 3 4 5
        // 6 7 8
        // Manually optimized search for matrices with determinant 1.
        for(i4=-range_;i4<=range_;i4++) {
        for(i8=-range_;i8<=range_;i8++) { i48 = i4*i8;
        for(i5=-range_;i5<=range_;i5++) {
        for(i7=-range_;i7<=range_;i7++) { i4857 = i48-i5*i7;
        for(i3=-range_;i3<=range_;i3++) { i38 = i3*i8;
                                          i37 = i3*i7;
        for(i6=-range_;i6<=range_;i6++) { i3746 = i37-i4*i6;
                                          i3856 = i38-i5*i6;
          if (i3746 == 0) {
            if (i3856 == 0) {
              if (i4857 == -1 || i4857 == 1) {
                return_target_id_ = 1;
                i0 = i4857;
                for(i1=-range_;i1<=range_;i1++) {
                  for(i2=-range_;i2<=range_;i2++) {
                    return;
                    return_target_1:;
                  }
                }
              }
            }
            else {
              return_target_id_ = 2;
              for(i0=-range_;i0<=range_;i0++) {
                i04857_one = i0*i4857 - 1;
                i1 = i04857_one / i3856;
                if (-range_ <= i1 && i1 <= range_ && i1*i3856 == i04857_one) {
                  for(i2=-range_;i2<=range_;i2++) {
                    return;
                    return_target_2:;
                  }
                }
              }
            }
          }
          else {
            return_target_id_ = 3;
            for(i0=-range_;i0<=range_;i0++) {
              i04857 = i0*i4857;
              for(i1=-range_;i1<=range_;i1++) {
                one_i04857_13856 = 1 - i04857 + i1*i3856;
                i2 = one_i04857_13856 / i3746;
                if (-range_ <= i2 && i2 <= range_
                    && i2*i3746 == one_i04857_13856) {
                  return;
                  return_target_3:;
                }
              }
            }
          }
        }}}}}}
        at_end_ = true;
      }
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_UNIMODULAR_GENERATOR_H
