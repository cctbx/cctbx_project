#include <scitbx/math/eigensystem.h>

namespace scitbx { namespace math {

  /*! Principal axes of inertia given discrete points in 3-dimensional
      space and optionally weights.
   */
  template <typename FloatType=double>
  class principal_axes_of_inertia
  {
    private:
      vec3<FloatType> center_of_mass_;
      sym_mat3<FloatType> inertia_tensor_;
      eigensystem::real_symmetric<FloatType> eigensystem_;

    public:
      //! Default constructor. Some data members are not initialized!
      principal_axes_of_inertia() {}

      //! Intitialization given discrete points with unit weights.
      principal_axes_of_inertia(
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
            vec3<FloatType> pp(p[0]*p[0], p[1]*p[1], p[2]*p[2]);
            inertia_tensor_(0,0) += pp[1] + pp[2];
            inertia_tensor_(1,1) += pp[0] + pp[2];
            inertia_tensor_(2,2) += pp[0] + pp[1];
            inertia_tensor_(0,1) -= p[0] * p[1];
            inertia_tensor_(0,2) -= p[0] * p[2];
            inertia_tensor_(1,2) -= p[1] * p[2];
          }
        }
        eigensystem_ = math::eigensystem::real_symmetric<FloatType>(
          inertia_tensor_);
      }

      //! Intitialization given discrete points and weights.
      principal_axes_of_inertia(
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
            vec3<FloatType> pp(p[0]*p[0], p[1]*p[1], p[2]*p[2]);
            FloatType w = weights[i_p];
            inertia_tensor_(0,0) += w * (pp[1] + pp[2]);
            inertia_tensor_(1,1) += w * (pp[0] + pp[2]);
            inertia_tensor_(2,2) += w * (pp[0] + pp[1]);
            inertia_tensor_(0,1) -= w * p[0] * p[1];
            inertia_tensor_(0,2) -= w * p[0] * p[2];
            inertia_tensor_(1,2) -= w * p[1] * p[2];
          }
        }
        eigensystem_ = math::eigensystem::real_symmetric<FloatType>(
          inertia_tensor_);
      }

      //! Center of mass.
      /*! The weighted average of the coordinates of the points as passed
          to the constructor.
       */
      vec3<FloatType> const&
      center_of_mass() const { return center_of_mass_; }

      //! Real-symmetric 3x3 inertia tensor.
      /*! See e.g. http://kwon3d.com/theory/moi/iten.html or
          search for "inertia tensor" at http://www.google.com/
          or another search engine.
       */
      sym_mat3<FloatType> const&
      inertia_tensor() const { return inertia_tensor_; }

      //! Eigenvectors and eigenvalues of inertia_tensor().
      math::eigensystem::real_symmetric<FloatType> const&
      eigensystem() const { return eigensystem_; }
  };

}} // namespace scitbx::math
