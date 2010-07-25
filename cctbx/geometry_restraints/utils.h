#ifndef CCTBX_GEOMETRY_RESTRAINTS_UTILS_H
#define CCTBX_GEOMETRY_RESTRAINTS_UTILS_H

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/rt_mx.h>
#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>
#include <tbxx/optional_copy.hpp>

namespace cctbx { namespace geometry_restraints {

  using tbxx::optional_container;

  /*! \brief Rotation part of rt_mx.inverse() in the cartesian system.
   */
  /*! Useful for mapping gradients in cartesian space
      from the asymmetric unit to the original sites.

      Inefficient implementation, the rotation matrices are NOT cached.
   */
  inline
  scitbx::mat3<double>
  r_inv_cart(
    uctbx::unit_cell const& unit_cell,
    sgtbx::rt_mx const& rt_mx)
  {
    typedef scitbx::type_holder<double> t_h;
    scitbx::mat3<double> o(unit_cell.orthogonalization_matrix());
    scitbx::mat3<double> f(unit_cell.fractionalization_matrix());
    scitbx::mat3<double> r_inv_cart_
      = o*rt_mx.r().inverse().as_floating_point(t_h())*f;
    return r_inv_cart_;
  }

  /*! \brief Difference between angle_1 and angle_2 (in degrees) taking the
      periodicity into account.
   */
  inline
  double
  angle_delta_deg(double angle_1, double angle_2, int periodicity=1)
  {
    double half_period = 180./std::max(1,std::abs(periodicity));
    double d = std::fmod(angle_2-angle_1, 2*half_period);
    if      (d < -half_period) d += 2*half_period;
    else if (d >  half_period) d -= 2*half_period;
    return d;
  }

  namespace detail {

    template <typename ProxyType, typename RestraintType>
    struct generic_deltas
    {
      static
      af::shared<double>
      get(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        af::const_ref<ProxyType> const& proxies)
      {
        af::shared<double> result((af::reserve(proxies.size())));
        for(std::size_t i=0;i<proxies.size();i++) {
          ProxyType const& proxy = proxies[i];
          RestraintType restraint(sites_cart, proxy);
          result.push_back(restraint.delta);
        }
        return result;
      }

      static
      af::shared<double>
      get(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        af::const_ref<ProxyType> const& proxies)
      {
        af::shared<double> result((af::reserve(proxies.size())));
        for(std::size_t i=0;i<proxies.size();i++) {
          ProxyType const& proxy = proxies[i];
          RestraintType restraint(unit_cell, sites_cart, proxy);
          result.push_back(restraint.delta);
        }
        return result;
      }
    };

    template <typename ProxyType, typename RestraintType>
    struct generic_residuals
    {
      static
      af::shared<double>
      get(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        af::const_ref<ProxyType> const& proxies)
      {
        af::shared<double> result((af::reserve(proxies.size())));
        for(std::size_t i=0;i<proxies.size();i++) {
          ProxyType const& proxy = proxies[i];
          RestraintType restraint(sites_cart, proxy);
          result.push_back(restraint.residual());
        }
        return result;
      }

      static
      af::shared<double>
      get(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        af::const_ref<ProxyType> const& proxies)
      {
        af::shared<double> result((af::reserve(proxies.size())));
        for(std::size_t i=0;i<proxies.size();i++) {
          ProxyType const& proxy = proxies[i];
          RestraintType restraint(unit_cell, sites_cart, proxy);
          result.push_back(restraint.residual());
        }
        return result;
      }
    };

    template <typename ProxyType, typename RestraintType>
    struct generic_residual_sum
    {
      static
      double
      get(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        af::const_ref<ProxyType> const& proxies,
        af::ref<scitbx::vec3<double> > const& gradient_array)
      {
        CCTBX_ASSERT(   gradient_array.size() == 0
                     || gradient_array.size() == sites_cart.size());
        double result = 0;
        for(std::size_t i=0;i<proxies.size();i++) {
          ProxyType const& proxy = proxies[i];
          RestraintType restraint(sites_cart, proxy);
          result += restraint.residual();
          if (gradient_array.size() != 0) {
            restraint.add_gradients(gradient_array, proxy.i_seqs);
          }
        }
        return result;
      }

      static
      double
      get(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        af::const_ref<ProxyType> const& proxies,
        af::ref<scitbx::vec3<double> > const& gradient_array)
      {
        CCTBX_ASSERT(   gradient_array.size() == 0
                     || gradient_array.size() == sites_cart.size());
        double result = 0;
        for(std::size_t i=0;i<proxies.size();i++) {
          ProxyType const& proxy = proxies[i];
          RestraintType restraint(unit_cell, sites_cart, proxy);
          result += restraint.residual();
          if (gradient_array.size() != 0) {
            restraint.add_gradients(unit_cell, gradient_array, proxy);
          }
        }
        return result;
      }
    };

  } // namespace detail

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_UTILS_H
