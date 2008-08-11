#ifndef CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
#define CCTBX_MAPTBX_REAL_SPACE_REFINEMENT_H
// Done by Erik McKee && Jacob Smith

#include <cctbx/maptbx/basic_map.h>
#include <scitbx/array_family/versa_matrix.h>
#if !(defined(__linux__) && defined(__GNUC__) \
 && __GNUC__ == 2 && __GNUC_MINOR__ == 96)
#include <limits>
#endif

namespace cctbx { namespace maptbx { namespace real_space_refinement {

static const bool no_assert = true;

template < typename FloatType, typename IntType >
FloatType
residual(
  basic_map<FloatType,IntType> const& map,
  af::const_ref<scitbx::vec3<FloatType> > const& sites,
  af::const_ref<FloatType> const& weights)
{
  CCTBX_ASSERT(weights.size()==sites.size());

  return -af::matrix_multiply(weights,map.get_cart_values(sites).const_ref()) / af::sum(weights);
}

// computes the gradient
template < typename FloatType, typename IntType >
scitbx::vec3<FloatType>
gradient (
  basic_map<FloatType,IntType> const& map,
  scitbx::vec3<FloatType> const& site,
  FloatType const& delta_h )
{
  scitbx::vec3<FloatType> result;
  for ( std::size_t i=0; i<3; ++i )
  {
    cartesian<FloatType> h_pos(site);
    cartesian<FloatType> h_neg(site);
    h_pos[i] += delta_h;
    h_neg[i] -= delta_h;
    // "Numerical Analysis", Richard L. Burden, J. Douglas Faires
    // pg. 170:  f'(x) = (0.5*h)*(-f(x-h)+f(x+h))
    // error: -(h*h)*(1/6.0)*f'''(P)
    // The above book has an explanation of the "2", which is nonobvious (pgs.
    // 168-171)
    result[i] = (map[h_pos] - map[h_neg]) / 2*delta_h;
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

template < typename FloatType, typename IntType >
af::shared<scitbx::vec3<FloatType> >
gradients(
  basic_map<FloatType,IntType> const& map,
  af::const_ref<scitbx::vec3<FloatType> > const& sites,
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
  std::size_t num = sites.size();

  // compute the 3-point numerical derivative
  // and sign it in the proper direction
  for(std::size_t i = 0; i < num; ++i)
  {
    grad_val = gradient(map,sites[i],delta_h);

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
    FloatType current_density = map.get_cart_value(sites[i]);
    FloatType local_delta_h = delta_h;
    bool zero_out = true;
    for ( std::size_t iter=0; iter<search_iter; ++iter )
    {
      scitbx::vec3<FloatType> pred = sites[i] + grad_val;
      FloatType pred_dens = map.get_cart_value(pred);
      if ( pred_dens < current_density )
      {
        local_delta_h *= 0.5;
        grad_val = gradient(map,sites[i],local_delta_h);
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
