#include <cctbx/sgtbx/change_of_basis_op.h>

namespace cctbx { namespace sgtbx {

  tr_vec change_of_basis_op::operator()(
    tr_vec const& t,
    int sign_identity) const
  {
    // (C|V)( I|T)(C^-1|W)=( C|CT+V)(C^-1|W)=( I|CT+V+CW)=( I|C(T+W)+V)
    // (C|V)(-I|T)(C^-1|W)=(-C|CT+V)(C^-1|W)=(-I|CT+V-CW)=(-I|C(T-W)+V)
    tr_vec tf = t.new_denominator(c_inv_.t().den());
    tr_vec tw;
    if (sign_identity >= 0) tw = tf + c_inv_.t();
    else                    tw = tf - c_inv_.t();
    return (  c_.r() * tw
            + c_.t().scale(c_.r().den())).new_denominator(t.den());
  }

  rot_mx change_of_basis_op::operator()(rot_mx const& r) const
  {
    CCTBX_ASSERT(r.den() == 1);
    return (c_.r() * r * c_inv_.r()).new_denominator(1);
  }

  rt_mx change_of_basis_op::operator()(rt_mx const& s) const
  {
    CCTBX_ASSERT(s.r().den() == 1);
    CCTBX_ASSERT(c_.t().den() % s.t().den() == 0);
    return (c_ * (s.scale(1, c_.t().den() / s.t().den()) * c_inv_))
      .new_denominators(s);
  }

  rt_mx change_of_basis_op::apply(rt_mx const& s) const
  {
    return c_.multiply(s.multiply(c_inv_));
  }

  af::shared<miller::index<> >
  change_of_basis_op::apply(
    af::const_ref<miller::index<> > const& miller_indices) const
  {
    af::shared<miller::index<> > result((af::reserve(miller_indices.size())));
    for(std::size_t i=0;i<miller_indices.size();i++) {
      result.push_back(apply(miller_indices[i]));
    }
    return result;
  }

}} // namespace cctbx::sgtbx
