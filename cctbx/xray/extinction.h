#ifndef CCTBX_XRAY_EXTINCTION_H
#define CCTBX_XRAY_EXTINCTION_H
#include <cctbx/uctbx.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/sparse/matrix.h>
#include <cctbx/xray/fc_correction.h>

namespace cctbx { namespace xray {

  /* shelx-like extinction correction, shelx manual 7-7
  Fc = Fc*(1+0.001*x*Fc^2*lambda^3/sin(2theta))^0.25,
  where x is the extinction coefficient
  */
  template <typename FloatType>
  struct shelx_extinction_correction : public fc_correction<FloatType> {
    typedef fc_correction<FloatType> fc_correction_t;

    shelx_extinction_correction(uctbx::unit_cell const &u_cell,
      FloatType wavelength, FloatType value)
      : u_cell(u_cell),
      wavelength(wavelength),
      value(value),
      grad_value(0), grad_fc_mod(1),
      grad_index(-1)
    {}

    virtual FloatType compute(miller::index<> const &h,
      FloatType fc_sq,
      bool compute_grad) const
    {
      const FloatType
        p = calc_factor(h, fc_sq),
        k = 1 + p * value,
        k_sqrt = std::sqrt(k),
        rv = 1 / k_sqrt;
      if (compute_grad) {
        grad_value = -fc_sq * p / (2 * k * k_sqrt);
        // Laura Midgley, Oct 2020
        //grad_fc_mod = 0.5*(rv + rv*rv*rv);
      }
      return rv;
    }
    virtual int get_grad_index() const { return grad_index; }
    virtual af::const_ref<FloatType> get_gradients() const {
      return af::const_ref<FloatType>(&grad_value, 1);
    }
    virtual FloatType get_grad_Fc_multiplier() const { return grad_fc_mod; }
    virtual boost::shared_ptr<fc_correction_t> fork() const {
      shelx_extinction_correction *cr =
        new shelx_extinction_correction(u_cell, wavelength, value);
      cr->grad = fc_correction_t::grad;
      cr->grad_index = grad_index;
      return boost::shared_ptr<fc_correction_t>(cr);
    }
    virtual size_t n_param() const { return 1; }
  protected:
    FloatType calc_factor(const miller::index<>& h,
      FloatType fc_sq) const
    {
      const FloatType sin_2t = u_cell.sin_two_theta(h, wavelength);
      return fc_sq*std::pow(wavelength,3)/(sin_2t*1000);
    }

    uctbx::unit_cell const &u_cell;
  public:
    FloatType wavelength, value;
    mutable FloatType grad_value, grad_fc_mod;
    int grad_index;
  };

  /* shelx-like extinction correction, shelx manual
  Fc_sq = Fc_sq*(1 - g*exp(-8 * PI^2 * U * stol_sq))^2
  coefficient
  */
  template <typename FloatType>
  struct shelx_SWAT_correction : public fc_correction<FloatType> {
    typedef fc_correction<FloatType> fc_correction_t;

    shelx_SWAT_correction(uctbx::unit_cell const& u_cell,
      FloatType g, FloatType U)
      : u_cell(u_cell),
      grad_index(-1)
    {
      values[0] = g;
      values[1] = U;
    }

    // for testing
    FloatType calc(miller::index<> const& h, FloatType Fsq, FloatType g, FloatType U) const {
      FloatType stol_sq = u_cell.stol_sq(h);
      return Fsq * std::pow(1-g*std::exp(-scitbx::constants::eight_pi_sq * stol_sq *U), 2);
    }

    virtual FloatType compute(miller::index<> const& h,
      FloatType fc_sq,
      bool compute_grad) const
    {
      using namespace scitbx::constants;
      FloatType stol_sq = u_cell.stol_sq(h);
      FloatType f = eight_pi_sq * stol_sq;
      FloatType ef = exp(-f * values[1]);
      FloatType rv = 1 - values[0] * ef;
      if (compute_grad) {
        FloatType dm = 2 * fc_sq * rv * ef;
        grads[0] = -dm;
        grads[1] = dm * values[0] * f;
        /* Testing derivatives shows consistensy with analytical expressions above
        FloatType eps = 1e-6;
        FloatType v1 = calc(h, fc_sq, values[0] - eps, values[1]);
        FloatType v2 = calc(h, fc_sq, values[0] + eps, values[1]);
        FloatType diff = (v2 - v1) / (2*eps);

        v1 = calc(h, fc_sq, values[0], values[1] - eps);
        v2 = calc(h, fc_sq, values[0], values[1] + eps);
        diff = (v2 - v1) / (2 * eps);
        v1 = 0; // allow for a breakpoint here
        */
      }
      return rv*rv;
    }

    FloatType get_g() const { return values[0]; }
    void set_g(FloatType v) { values[0] = v; }

    FloatType get_U() const { return values[1]; }
    void set_U(FloatType v) { values[1] = v; }

    virtual int get_grad_index() const { return grad_index; }
    virtual af::const_ref<FloatType> get_gradients() const {
      return grads.const_ref();
    }
    virtual FloatType get_grad_Fc_multiplier() const { return 1; }
    virtual boost::shared_ptr<fc_correction_t> fork() const {
      shelx_SWAT_correction* cr =
        new shelx_SWAT_correction(u_cell, values[0], values[1]);
      cr->grad = fc_correction_t::grad;
      cr->grad_index = grad_index;
      return boost::shared_ptr<fc_correction_t>(cr);
    }
    virtual size_t n_param() const { return 2; }
  protected:
    uctbx::unit_cell const &u_cell;
  public:
    af::tiny<FloatType, 2> values; // g, U
    mutable af::tiny<FloatType, 2> grads; // d_Fc_d_g, d_Fc_d_U;
    int grad_index;
  };
}} // namespace cctbx::xray

#endif // GUARD
