#ifndef CCTBX_XRAY_TWIN_COMPONENT_H
#define CCTBX_XRAY_TWIN_COMPONENT_H

#include <cctbx/sgtbx/rot_mx.h>
#include <scitbx/array_family/shared.h>

namespace cctbx { namespace xray {

  template <typename FloatType>
  class twin_component
  {
  public:
    //! Default constructor. Data members are not initialized!
    twin_component() {}

    twin_component(
      sgtbx::rot_mx const &twin_law_,
      FloatType const &twin_fraction_,
      bool grad_twin_fraction_=false)
      :
    twin_law(twin_law_),
    twin_fraction(twin_fraction_),
    grad_twin_fraction(grad_twin_fraction_),
    grad_index(-1)
    {}

    void set_grad_twin_fraction(bool grad_twin_fraction_=true)
    {
      grad_twin_fraction = grad_twin_fraction_;
    }

    sgtbx::rot_mx twin_law;
    FloatType twin_fraction;
    int grad_index;
    bool grad_twin_fraction;
  };

  template <typename FloatType>
  FloatType sum_twin_fractions(
    af::shared<twin_component<FloatType> *> twin_components)
  {
    FloatType result = 0.;
    for (std::size_t i_twin=0; i_twin<twin_components.size(); i_twin++) {
      result += twin_components[i_twin]->twin_fraction;
    }
    return result;
  }

  template <typename FloatType>
  void set_grad_twin_fraction(
    af::shared<twin_component<FloatType> *> twin_components,
    bool grad_twin_fraction=true)
  {
    for (std::size_t i_twin=0; i_twin<twin_components.size(); i_twin++) {
      twin_components[i_twin]->set_grad_twin_fraction(grad_twin_fraction);
    }
  }

}} // namespace cctbx::xray

#endif // GUARD
