#ifndef CCTBX_XRAY_EXTINCTION_H
#define CCTBX_XRAY_EXTINCTION_H
#include <cctbx/uctbx.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/sparse/matrix.h>

namespace cctbx { namespace xray {

  /// extinction correction iterface
  template <typename FloatType>
  struct extinction_correction {
    virtual ~extinction_correction() {}
    // return multiplier for Fc_sq
    virtual FloatType compute(const miller::index<>& h,
      FloatType Fc_sq,
      af::shared<FloatType> &gradient,
      bool calc_grad) = 0;
    virtual FloatType& get_value() = 0;
    virtual bool grad_value() const = 0;
  };

  /// dummy extinction implementation
  template <typename FloatType>
  struct dummy_extinction_correction : public extinction_correction<FloatType> {
    virtual FloatType compute(const miller::index<>& h,
      FloatType Fc_sq,
      af::shared<FloatType> &gradient,
      bool calc_grad) { return 1; }

    virtual FloatType& get_value() {
      static FloatType value = 0;
      return value;
    }
    virtual bool grad_value() const { return false; }
  };

  /* shelx-like extinction correction, shelx manual 7-7
  Fc = Fc*(1+0.001*x*Fc^2/sin(2theta))^0.25, where x is the extinction
  coefficient
  */
  template <typename FloatType>
  struct shelx_extinction_correction : extinction_correction<FloatType> {
    shelx_extinction_correction(uctbx::unit_cell const &_u_cell,
      FloatType _lambda, FloatType _value)
      : u_cell(_u_cell),
        lambda(_lambda),
        value(_value),
        grad_index(-1),
        grad(false) {}

    FloatType compute(const miller::index<>& h,
      FloatType Fc_sq,
      af::shared<FloatType> &gradient,
      bool calc_grad)
    {
      const FloatType
        p = calc_factor(h, Fc_sq),
        k = 1+p*value,
        k_sqrt = std::sqrt(k);
      if( calc_grad && grad ) {
        gradient[grad_index] += Fc_sq*(1-p/(2*k*k_sqrt)); // pow(k,-3/2)
      }
      return 1/k_sqrt;
    }
    virtual FloatType& get_value() { return value; }
    virtual bool grad_value() const { return grad; }
  protected:
    FloatType calc_factor(const miller::index<>& h,
      FloatType Fc_sq) const
    {
      const FloatType sin_2t = std::abs(std::sin(
        uctbx::d_star_sq_as_two_theta(u_cell.d_star_sq(h), lambda)));
      return Fc_sq*std::pow(lambda,3)/(sin_2t*1000);
    }

    uctbx::unit_cell const &u_cell;
  public:
    FloatType lambda, value;
    int grad_index;
    bool grad;
  };

}} // namespace cctbx::xray

#endif // GUARD
