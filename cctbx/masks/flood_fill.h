#ifndef CCTBX_MASKS_FLOOD_FILL_H
#define CCTBX_MASKS_FLOOD_FILL_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/uctbx.h>
#include <scitbx/array_family/accessors/c_grid_periodic.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/math/accumulators.h>

#include <stack>

namespace cctbx { namespace masks {

  namespace af = scitbx::af;

  //! Separates a mask into individual voids.
  /*!
     Takes a mask of 0's and 1's where 0 indicates outside a void and 1 is
     inside the void.  Individual voids will be labelled as 2, 3, etc...

     This is not the most efficient implementation, but the time spent in
     this function is small relative to that taken computing the mask in the
     first place.
   */
  template <typename DataType, typename FloatType>
  class flood_fill
  {
    //! Convenience typedefs
    typedef scitbx::sym_mat3<FloatType> sym_mat3;
    typedef scitbx::vec3<FloatType> vec3;

  public:
    flood_fill(
      af::ref<DataType, af::c_grid_periodic<3> > const & data,
      uctbx::unit_cell const & unit_cell)
      :
    n_voids_(0),
    gridding_n_real(data.accessor()),
    unit_cell_(unit_cell)
    {
      std::stack<index_t> stack;
      DataType target = 1;
      DataType replacement = 2;
      for (std::size_t i=0; i<gridding_n_real[0]; i++) {
        for (std::size_t j=0; j<gridding_n_real[1]; j++) {
          for (std::size_t k=0; k<gridding_n_real[2]; k++) {
            if (data(i,j,k) == target) {
              stack.push(index_t(i,j,k));
              data(i,j,k) = replacement;
              accumulators.push_back(accumulator_t());
              n_voids_ += 1;
              grid_points_per_void.push_back(0);
              while (!stack.empty()) {
                index_t index = stack.top();
                stack.pop();
                accumulators[accumulators.size() -1](vec3(index));
                grid_points_per_void[grid_points_per_void.size() - 1]++;
                for (std::size_t i=0; i<3; i++) {
                  index_t index_1 = index;
                  index_1[i]++;
                  if (data(index_1) == target) {
                    data(index_1) = replacement;
                    stack.push(index_1);
                  }
                  index_1[i] = index[i]-1;
                  if (data(index_1) == target) {
                    data(index_1) = replacement;
                    stack.push(index_1);
                  }
                }
              }
              replacement++;
            }
          }
        }
      }
    }

    unsigned n_voids() { return n_voids_; }

    //! Provides the centre of mass or average indices for each void.
    /*! Move coordinates (if necessary) so the fractional coordinates will be
        between -1 and 1. Where there is a channel, an index along the
        direction of the channel may drift many unit cells along. In the case
        of a channel this index (along the channel) is not really applicable
        anyway.
     */
    af::shared<vec3> centres_of_mass()
    {
      af::shared<vec3> result((af::reserve(n_voids_)));
      for (std::size_t i=0; i<n_voids_; i++) {
        result.push_back(accumulators[i].center_of_mass());
        for (std::size_t j=0; j<3; j++) {
          while (result[i][j] > gridding_n_real[j]) {
            result[i][j] -= gridding_n_real[j];
          }
          while (result[i][j] < -gridding_n_real[j]) {
            result[i][j] += gridding_n_real[j];
          }
        }
      }
      return result;
    }

    //! Provides the centre of mass in fractional coordinates for each void.
    af::shared<vec3> centres_of_mass_frac()
    {
      af::shared<vec3> result = centres_of_mass();
      for (std::size_t i=0; i<n_voids_; i++) {
        result[i] = result[i] / vec3(gridding_n_real);
      }
      return result;
    }

    //! Provides the centre of mass in Cartesian coordinates for each void.
    af::shared<vec3> centres_of_mass_cart()
    {
      return unit_cell_.orthogonalize(centres_of_mass_frac().const_ref());
    }

    /*! The covariance matrix has been accumulated for the grid indices,
        and this must now be converted to fractional coordinates.
     */
    af::shared<sym_mat3> covariance_matrices_frac()
    {
      af::shared<sym_mat3> result((af::reserve(n_voids_)));
      for (std::size_t i=0; i<n_voids_; i++) {
        sym_mat3 cov = accumulators[i].covariance_matrix();
        cov[0] /= (gridding_n_real[0] * gridding_n_real[0]);
        cov[1] /= (gridding_n_real[1] * gridding_n_real[1]);
        cov[2] /= (gridding_n_real[2] * gridding_n_real[2]);
        cov[3] /= (gridding_n_real[0] * gridding_n_real[1]);
        cov[4] /= (gridding_n_real[0] * gridding_n_real[2]);
        cov[5] /= (gridding_n_real[1] * gridding_n_real[2]);
        sym_mat3 G = unit_cell_.metrical_matrix();
        cov = sym_mat3(G * cov * G, 1e-6);
        result.push_back(cov);
      }
      return result;
    }

    af::shared<sym_mat3> covariance_matrices_cart()
    {
      af::shared<sym_mat3> result
        = covariance_matrices_frac();
      scitbx::mat3<FloatType> F = unit_cell_.fractionalization_matrix();
      for (std::size_t i=0; i<n_voids_; i++) {
        result[i] = sym_mat3(F.transpose() * result[i] * F, 1e-6);
      }
      return result;
    }

    //! The inertia tensor in fractional coordinates for each void.
    af::shared<sym_mat3> inertia_tensors_frac()
    {
       return inertia_tensors_impl(covariance_matrices_frac());
    }

    //! The inertia tensor in Cartesian coordinates for each void.
    af::shared<sym_mat3> inertia_tensors_cart()
    {
       return inertia_tensors_impl(covariance_matrices_cart());
    }

    uctbx::unit_cell const & unit_cell() const { return unit_cell_; }

    typedef af::c_grid_periodic<3>::index_type index_t;
    af::shared<int> grid_points_per_void;
    index_t const gridding_n_real;

  private:
    /*! We obtain the inertia tensor from the relationship:

          inertia_tensor = sum_weights * (identity * trace(covariance) - covariance)

        see also: http://en.wikipedia.org/wiki/Variance#Moment_of_inertia
     */
    af::shared<sym_mat3> inertia_tensors_impl(
      af::shared<sym_mat3> const & covariance_matrices)
    {
      af::shared<sym_mat3> result = covariance_matrices;
      for (std::size_t i=0; i<n_voids_; i++) {
        scitbx::sym_mat3<FloatType> cov = result[i];
        FloatType trace = cov.trace();
        result[i] = sym_mat3(trace,trace,trace,0,0,0);
        result[i] -= cov;
        result[i] *= static_cast<FloatType>(grid_points_per_void[i]);
      }
      return result;
    }

    typedef scitbx::math::accumulator::inertia_accumulator<FloatType>
      accumulator_t;
    unsigned n_voids_;
    af::shared<accumulator_t> accumulators;
    uctbx::unit_cell const unit_cell_;
  };

}} // namespace cctbx::masks

#endif // CCTBX_MASKS_FLOOD_FILL_H
