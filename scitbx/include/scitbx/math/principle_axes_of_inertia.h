#include <scitbx/math/eigensystem.h>

namespace scitbx { namespace math {

  template <typename FloatType=double>
  class principle_axes_of_inertia
  {
    private:
      vec3<FloatType> center_of_mass_;
      sym_mat3<FloatType> inertia_tensor_;
      eigensystem::real_symmetric<FloatType> eigensystem_;

    public:
      principle_axes_of_inertia() {}

      principle_axes_of_inertia(
        af::const_ref<vec3<FloatType> > const& points)
      :
        center_of_mass_(0,0,0),
        inertia_tensor_(0,0,0,0,0,0)
      {
        if (points.size() != 0) {
          for(std::size_t i_p=0;i_p<points.size();i_p++) {
            center_of_mass_ += points[i_p];
          }
          center_of_mass_ /= static_cast<FloatType>(points.size());
          for(std::size_t i_p=0;i_p<points.size();i_p++) {
            vec3<FloatType> p = points[i_p] - center_of_mass_;
            vec3<FloatType> pp;
            for(int i=0;i<3;i++) {
              pp[i] = p[i] * p[i];
            }
            FloatType pp_sum = pp.sum();
            for(int i=0;i<3;i++) {
              inertia_tensor_[i] += pp_sum - pp[i];
            }
            inertia_tensor_(0,1) -= p[0] * p[1];
            inertia_tensor_(0,2) -= p[0] * p[2];
            inertia_tensor_(1,2) -= p[1] * p[2];
          }
        }
        eigensystem_ = math::eigensystem::real_symmetric<FloatType>(
          inertia_tensor_);
      }

      principle_axes_of_inertia(
        af::const_ref<vec3<FloatType> > const& points,
        af::const_ref<FloatType> const& weights)
      :
        center_of_mass_(0,0,0),
        inertia_tensor_(0,0,0,0,0,0)
      {
        SCITBX_ASSERT(weights.size() == points.size());
        FloatType sum_weights = 0;
        for(std::size_t i_p=0;i_p<points.size();i_p++) {
          FloatType w = weights[i_p];
          center_of_mass_ += w * points[i_p];
          sum_weights += w;
        }
        if (sum_weights != 0) {
          center_of_mass_ /= sum_weights;
          for(std::size_t i_p=0;i_p<points.size();i_p++) {
            vec3<FloatType> p = points[i_p] - center_of_mass_;
            vec3<FloatType> pp;
            for(int i=0;i<3;i++) {
              pp[i] = p[i] * p[i];
            }
            FloatType pp_sum = pp.sum();
            FloatType w = weights[i_p];
            for(int i=0;i<3;i++) {
              inertia_tensor_[i] += w * (pp_sum - pp[i]);
            }
            inertia_tensor_(0,1) -= w * p[0] * p[1];
            inertia_tensor_(0,2) -= w * p[0] * p[2];
            inertia_tensor_(1,2) -= w * p[1] * p[2];
          }
        }
        eigensystem_ = math::eigensystem::real_symmetric<FloatType>(
          inertia_tensor_);
      }

      vec3<FloatType> const&
      center_of_mass() const { return center_of_mass_; }

      sym_mat3<FloatType> const&
      inertia_tensor() const { return inertia_tensor_; }

      math::eigensystem::real_symmetric<FloatType> const&
      eigensystem() const { return eigensystem_; }
  };

}} // namespace scitbx::math
