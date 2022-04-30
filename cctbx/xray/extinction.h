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
        //grad_fc_mod = rv - fc_sq*value/(2 * k * k_sqrt);
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
  Fc = Fc*(1 - g*exp(-8 * M_PI^2 * U * stol_sq))
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

    virtual FloatType compute(miller::index<> const& h,
      FloatType fc_sq,
      bool compute_grad) const
    {
      using namespace scitbx::constants;
      FloatType stol_sq = u_cell.stol_sq(h);
      FloatType ef = exp(-eight_pi_sq * values[1] * stol_sq);
      if (compute_grad) {
        grads[0] = -ef;
        grads[1] = -eight_pi_sq * stol_sq * ef;
      }
      return 1 - values[0] * ef;
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
  protected:
    uctbx::unit_cell const &u_cell;
  public:
    af::tiny<FloatType, 2> values; // g, U
    mutable af::tiny<FloatType, 2> grads; // d_Fc_d_g, d_Fc_d_U;
    int grad_index;
  };
}} // namespace cctbx::xray

#endif // GUARD
