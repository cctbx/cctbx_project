#ifndef SCITBX_RIGID_BODY_JOINT_T_H
#define SCITBX_RIGID_BODY_JOINT_T_H

#include <scitbx/rotr3.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/accessors/mat_grid.h>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>

namespace scitbx { namespace rigid_body {

  using boost::shared_ptr;

  //! Change-of-basis ("cb") matrix for joint alignment and its inverse.
  template <typename FloatType>
  struct alignment_t
  {
    //! global frame -> body frame
    rotr3<FloatType> cb_0b;
    //! body frame -> global frame
    rotr3<FloatType> cb_b0;

    alignment_t() {}

    alignment_t(
      rotr3<FloatType> const& cb_0b_,
      rotr3<FloatType> const& cb_b0_)
    :
      cb_0b(cb_0b_),
      cb_b0(cb_b0_)
    {}
  };

  //! Abstract joint model.
  template <typename FloatType>
  struct joint_t
  {
    typedef FloatType ft;

    unsigned degrees_of_freedom;
    unsigned q_size;

    //! Xj = cb_as_spatial_transform(cb_ps)
    rotr3<ft> cb_ps;
    rotr3<ft> cb_sp;

    joint_t(
      unsigned degrees_of_freedom_,
      unsigned q_size_)
    :
      degrees_of_freedom(degrees_of_freedom_),
      q_size(q_size_)
    {}

    virtual
    ~joint_t() {}

    virtual
    af::const_ref<ft>
    qd_zero() const = 0;

    virtual
    af::const_ref<ft>
    qdd_zero() const = 0;

    //! S in RBDA Tab. 4.1, p. 79.
    virtual
    af::const_ref<ft, af::mat_grid>
    motion_subspace() const = 0;

    virtual
    boost::optional<vec3<ft> >
    get_linear_velocity(
      af::const_ref<ft> const& qd) const = 0;

    virtual
    af::small<ft, 6>
    new_linear_velocity(
      af::const_ref<ft> const& qd,
      vec3<ft> const& value) const = 0;

    virtual
    shared_ptr<joint_t>
    time_step_position(
      af::const_ref<ft> const& qd,
      ft const& delta_t) const = 0;

    virtual
    af::small<ft, 6>
    time_step_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& qdd,
      ft const& delta_t) const = 0;

    virtual
    af::small<ft, 7>
    tau_as_d_e_pot_d_q(
      af::small<ft, 6> const& tau) const = 0;

    virtual
    af::small<ft, 7>
    get_q() const = 0;

    virtual
    shared_ptr<joint_t>
    new_q(
      af::const_ref<ft> const& q) const = 0;
  };

}} // namespace scitbx::rigid_body

#endif // GUARD
