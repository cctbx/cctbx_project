#ifndef CCTBX_XRAY_EXTINCTION_H
#define CCTBX_XRAY_EXTINCTION_H
#include <cctbx/uctbx.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/sparse/matrix.h>

namespace cctbx { namespace xray {

  /// extinction correction iterface
  template <typename FloatType>
  struct extinction_correction {
    virtual ~extinction_correction() {}
    // return multiplier for Fc_sq
    virtual af::tiny<FloatType,2> compute(miller::index<> const &h,
      FloatType fc_sq,
      bool compute_gradient) = 0;
    virtual FloatType& get_value() = 0;
    virtual bool grad_value() const = 0;
    virtual int get_grad_index() const = 0;
    static af::tiny<FloatType,2> build_return_value(FloatType v, FloatType g=0) {
      af::tiny<FloatType,2> rv;
      rv[0] = v;
      rv[1] = g;
      return rv;
    }
  };

  /// dummy extinction implementation
  template <typename FloatType>
  struct dummy_extinction_correction : public extinction_correction<FloatType> {
    dummy_extinction_correction() {}
    virtual af::tiny<FloatType,2> compute(miller::index<> const &h,
      FloatType fc_sq,
      bool compute_gradient)
      {
        return extinction_correction<FloatType>::build_return_value(1);
      }

    virtual FloatType& get_value() {
      static FloatType value = 0;
      return value;
    }
    virtual bool grad_value() const { return false; }
    virtual int get_grad_index() const { return -1; }
  };

  /* shelx-like extinction correction, shelx manual 7-7
  Fc = Fc*(1+0.001*x*Fc^2/sin(2theta))^0.25, where x is the extinction
  coefficient
  */
  template <typename FloatType>
  struct shelx_extinction_correction : public extinction_correction<FloatType> {
    typedef extinction_correction<FloatType> exti_t;
    shelx_extinction_correction(uctbx::unit_cell const &u_cell_,
      FloatType wavelength_, FloatType value_)
      : u_cell(u_cell_),
        wavelength(wavelength_),
        value(value_),
        grad_index(-1),
        grad(false) {}

    virtual af::tiny<FloatType,2> compute(miller::index<> const &h,
      FloatType fc_sq,
      bool compute_grad)
    {
      const FloatType
        p = calc_factor(h, fc_sq),
        k = 1+p*value,
        k_sqrt = std::sqrt(k);
      if( grad && compute_grad ) {
        return exti_t::build_return_value(1/k_sqrt, -fc_sq*p/(2*k*k_sqrt));
      }
      return exti_t::build_return_value(1/k_sqrt);
    }
    virtual FloatType& get_value() { return value; }
    virtual bool grad_value() const { return grad; }
    virtual int get_grad_index() const { return grad_index; }
  protected:
    FloatType calc_factor(const miller::index<>& h,
      FloatType fc_sq) const
    {
      const FloatType sin_2t = u_cell.sin_two_theta(h, wavelength);
      return fc_sq*std::pow(wavelength,3)/(sin_2t*1000);
    }

    uctbx::unit_cell u_cell;
  public:
    FloatType wavelength, value;
    int grad_index;
    bool grad;
  };

}} // namespace cctbx::xray

#endif // GUARD
