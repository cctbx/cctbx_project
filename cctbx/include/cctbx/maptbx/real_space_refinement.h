#ifndef CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
#define CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
// Done by Erik McKee && Jacob Smith

#include <cctbx/maptbx/eight_point_interpolation.h>
#if !(defined(__linux__) && defined(__GNUC__) \
 && __GNUC__ == 2 && __GNUC_MINOR__ == 96)
#include <limits>
#endif

namespace cctbx { namespace maptbx { namespace real_space_refinement {

static const bool no_assert = true;

template <typename FloatType>
FloatType
residual(
  af::const_ref<FloatType, af::flex_grid<> > const& map,
  scitbx::mat3<FloatType> const& gridding_matrix,
  af::const_ref<scitbx::vec3<FloatType> > const& sites_cart,
  af::const_ref<FloatType> const& weights)
{
  FloatType val = 0.0;
  std::size_t num = sites_cart.size();

  CCTBX_ASSERT(weights.size()==num);

  for(std::size_t i = 0; i < num; ++i) {
    val -= weights[i] * cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map, gridding_matrix, sites_cart[i], no_assert);
  }
  return val / af::sum(weights);
}

// computes the gradient
template < typename FloatType >
scitbx::vec3<FloatType>
gradient (
  af::const_ref<FloatType, af::flex_grid<> > const& map,
  scitbx::mat3<FloatType> const& gridding_matrix,
  scitbx::vec3<FloatType> const& pos,
  FloatType const& delta_h )
{
  scitbx::vec3<FloatType> result;
  for ( std::size_t i=0; i<3; ++i )
  {
    scitbx::vec3<FloatType> h_pos = pos;
    scitbx::vec3<FloatType> h_neg = pos;
    h_pos[i] += delta_h;
    h_neg[i] -= delta_h;
        FloatType e_h_pos = cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map,gridding_matrix,h_pos,no_assert);
    FloatType e_h_neg = cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map,gridding_matrix,h_neg,no_assert);
    result[i] = e_h_pos - e_h_neg;
    // "Numerical Analysis", Richard L. Burden, J. Douglas Faires
    // pg. 170:  f'(x) = (0.5*h)*(-f(x-h)+f(x+h))
    // error: -(h*h)*(1/6.0)*f'''(P)
    // The above book has an explanation of the "2", which is nonobvious (pgs.
    // 168-171)
    result[i] /= 2*delta_h;
  }
  // make sure the length of the gradient vector is no shorter than or
  // equal to delta-h
  FloatType len = result.length();
  if ( len > delta_h )
  {
    result *= delta_h / len;
  }
  return result;
}

template <typename FloatType>
af::shared<scitbx::vec3<FloatType> >
gradients(
  af::const_ref<FloatType, af::flex_grid<> > const& map,
  scitbx::mat3<FloatType> const& gridding_matrix,
  af::const_ref<scitbx::vec3<FloatType> > const& sites_cart,
  FloatType delta_h=1.0,
  std::size_t search_iter=0)
{
  if (search_iter == 0) {
#if defined(__linux__) && defined(__GNUC__) \
 && __GNUC__ == 2 && __GNUC_MINOR__ == 96
    search_iter = 53;
#else
    search_iter = std::numeric_limits<FloatType>::digits;
#endif
  }
  af::shared<scitbx::vec3<FloatType> > grad_vals;
  scitbx::vec3<FloatType> grad_val;
  std::size_t num = sites_cart.size();

  // compute the 3-point numerical derivative
  // and sign it in the proper direction
  for(std::size_t i = 0; i < num; ++i)
  {
    grad_val = gradient(map,gridding_matrix,sites_cart[i],delta_h);
    FloatType len = grad_val.length();

    // Full-Newton Predictor-Corrector With Early Abort
    // Get energy at predicted value
    // If predicted energy is less than current density
    //   halve the delta-h and get NEW slope value
    //   make sure new gradient-value is less than or equal
    //   to new delta-h
    // If the predicted energy is equal to or greater than
    // the current energy
    //    set the zero-out flag to false
    // If the zero-out flag is true, zero out the gradient
    // Add to the list of gradients
    FloatType current_density = cctbx::maptbx::
      non_crystallographic_eight_point_interpolation<FloatType>(
      map,gridding_matrix,sites_cart[i],no_assert);
    FloatType local_delta_h = delta_h;
    bool zero_out = true;
    for ( std::size_t iter=0; iter<search_iter; ++iter )
    {
      scitbx::vec3<FloatType> pred = sites_cart[i] + grad_val;
      FloatType pred_dens = cctbx::maptbx::
        non_crystallographic_eight_point_interpolation<FloatType>(
        map,gridding_matrix,pred,no_assert);
      if ( pred_dens < current_density )
      {
        local_delta_h *= 0.5;
        grad_val = gradient(map,gridding_matrix,sites_cart[i],local_delta_h);
      }
      else
      {
        zero_out = false;
        break;
      }
    }
    if ( zero_out )
    {
      grad_val[0] = grad_val[1] = grad_val[2] = 0;
    }
    grad_vals.push_back(-grad_val);
  }
  return grad_vals;
}

}}} // namespace cctbx::maptbx::real_space_refinement

#endif // CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
