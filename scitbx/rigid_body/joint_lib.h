#ifndef SCITBX_RIGID_BODY_JOINT_LIB_H
#define SCITBX_RIGID_BODY_JOINT_LIB_H

#include <scitbx/rigid_body/matrix_helpers.h>
#include <scitbx/rotr3.h>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>

namespace scitbx { namespace rigid_body { namespace joint_lib {

  template <typename FloatType>
  mat3<FloatType>
  rbda_eq_4_12(
    af::tiny<FloatType, 4> const& q)
  {
    typedef FloatType ft;
    ft p0 = q[0];
    ft p1 = q[1];
    ft p2 = q[2];
    ft p3 = q[3];
    return mat3<ft>(
      p0*p0+p1*p1-0.5,   p1*p2+p0*p3,     p1*p3-p0*p2,
        p1*p2-p0*p3,   p0*p0+p2*p2-0.5,   p2*p3+p0*p1,
        p1*p3+p0*p2,     p2*p3-p0*p1,   p0*p0+p3*p3-0.5) * ft(2);
  }

  template <typename FloatType>
  af::tiny<FloatType, 4*3>
  rbda_eq_4_13(
    af::tiny<FloatType, 4> const& q)
  {
    typedef FloatType ft;
    ft p0 = q[0];
    ft p1 = q[1];
    ft p2 = q[2];
    ft p3 = q[3];
    ft coeffs[] = {
      -p1, -p2, -p3,
      p0, -p3, p2,
      p3, p0, -p1,
      -p2, p1, p0};
    return af::tiny<ft, 4*3>(coeffs, coeffs+4*3) * ft(0.5);
  }

  template <typename FloatType>
  af::tiny<FloatType, 4*4>
  d_unit_quaternion_d_qe_matrix(
    af::tiny<FloatType, 4> const& q)
  {
    typedef FloatType ft;
    ft p0 = q[0];
    ft p1 = q[1];
    ft p2 = q[2];
    ft p3 = q[3];
    ft p0s = p0*p0;
    ft p1s = p1*p1;
    ft p2s = p2*p2;
    ft p3s = p3*p3;
    ft n3 = std::sqrt(fn::pow3(p0s+p1s+p2s+p3s));
    ft c00 = p1s+p2s+p3s;
    ft c11 = p0s+p2s+p3s;
    ft c22 = p0s+p1s+p3s;
    ft c33 = p0s+p1s+p2s;
    ft c01 = -p0*p1;
    ft c02 = -p0*p2;
    ft c03 = -p0*p3;
    ft c12 = -p1*p2;
    ft c13 = -p1*p3;
    ft c23 = -p2*p3;
    ft coeffs[] = {
      c00, c01, c02, c03,
      c01, c11, c12, c13,
      c02, c12, c22, c23,
      c03, c13, c23, c33};
    return af::tiny<FloatType, 4*4>(coeffs, coeffs+4*4) / n3;
  }

  template <typename FloatType>
  struct alignment_base
  {
    //! global frame -> body frame
    rotr3<FloatType> cb_0b;
    //! body frame -> global frame
    rotr3<FloatType> cb_b0;

    alignment_base() {}

    alignment_base(
      rotr3<FloatType> const& cb_0b_,
      rotr3<FloatType> const& cb_b0_)
    :
      cb_0b(cb_0b_),
      cb_b0(cb_b0_)
    {}
  };

  template <typename FloatType>
  struct joint_base
  {
    typedef FloatType ft;

    unsigned degrees_of_freedom;
    unsigned q_size;

    //! Xj = cb_as_spatial_transform(cb_ps)
    rotr3<ft> cb_ps;
    rotr3<ft> cb_sp;

    joint_base(
      unsigned degrees_of_freedom_,
      unsigned q_size_)
    :
      degrees_of_freedom(degrees_of_freedom_),
      q_size(q_size_)
    {}

    virtual
    ~joint_base() {}

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
      af::const_ref<ft> const& value) const = 0;

    virtual
    boost::shared_ptr<joint_base>
    time_step_position(
      af::const_ref<ft> const& qd,
      ft delta_t) const = 0;

    virtual
    af::small<ft, 6>
    time_step_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& qdd,
      ft delta_t) const = 0;

    virtual
    af::small<ft, 7>
    tau_as_d_pot_d_q(
      af::const_ref<ft> const& tau) const = 0;

    virtual
    af::small<ft, 7>
    get_q() const = 0;

    virtual
    boost::shared_ptr<joint_base>
    new_q(
      af::const_ref<ft> const& q) const = 0;
  };

  template <typename FloatType=double>
  struct zero_dof_alignment : alignment_base<FloatType>
  {
    typedef FloatType ft;

    zero_dof_alignment()
    :
      alignment_base<ft>(rotr3<ft>::identity(), rotr3<ft>::identity())
    {}
  };

  template <typename FloatType=double>
  struct zero_dof : joint_base<FloatType>
  {
    typedef FloatType ft;

    zero_dof()
    :
      joint_base<ft>(0, 0)
    {
      this->cb_ps = rotr3<ft>::identity();
      this->cb_sp = this->cb_ps;
    }

    virtual
    af::const_ref<ft>
    qd_zero() const { return af::const_ref<ft>(0, 0); }

    virtual
    af::const_ref<ft>
    qdd_zero() const { return af::const_ref<ft>(0, 0); }

    virtual
    af::const_ref<ft, af::mat_grid>
    motion_subspace() const
    {
      return af::const_ref<ft, af::mat_grid>(0, af::mat_grid(6, 0));
    }

    virtual
    boost::optional<vec3<ft> >
    get_linear_velocity(
      af::const_ref<ft> const& qd) const
    {
      return boost::optional<vec3<ft> >();
    }

    virtual
    af::small<ft, 6>
    new_linear_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& value) const
    {
      return af::small<ft, 6>(0);
    }

    virtual
    boost::shared_ptr<joint_base<ft> >
    time_step_position(
      af::const_ref<ft> const& qd,
      ft delta_t) const
    {
      return boost::shared_ptr<joint_base<ft> >(new zero_dof());
    }

    virtual
    af::small<ft, 6>
    time_step_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& qdd,
      ft delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 0);
      SCITBX_ASSERT(qdd.size() == 0);
      return af::small<ft, 6>(0);
    }

    virtual
    af::small<ft, 7>
    tau_as_d_pot_d_q(
      af::const_ref<ft> const& tau) const
    {
      return af::small<ft, 7>(0);
    }

    virtual
    af::small<ft, 7>
    get_q() const
    {
      return af::small<ft, 7>(0);
    }

    virtual
    boost::shared_ptr<joint_base<ft> >
    new_q(
      af::const_ref<ft> const& q) const
    {
      return boost::shared_ptr<joint_base<ft> >(new zero_dof());
    }
  };

  template <typename FloatType=double>
  struct six_dof_alignment : alignment_base<FloatType>
  {
    typedef FloatType ft;

    six_dof_alignment(
      vec3<ft> const& center_of_mass)
    :
      alignment_base<ft>(
        rotr3<ft>(mat3<ft>(1,1,1), -center_of_mass),
        rotr3<ft>(mat3<ft>(1,1,1), center_of_mass))
    {}
  };

  template <typename FloatType=double>
  struct six_dof : joint_base<FloatType>
  {
    typedef FloatType ft;

    af::tiny<ft, 4> qe;
    vec3<ft> qr;
    af::tiny<ft, 4> unit_quaternion;
    mat3<ft> e;

    six_dof(
      af::tiny<ft, 4> const& qe_,
      vec3<ft> const& qr_)
    :
      joint_base<ft>(6, 7),
      qe(qe_),
      qr(qr_),
      unit_quaternion(vec4_normalize(qe_)), // RBDA, bottom of p. 86
      e(rbda_eq_4_12(unit_quaternion))
    {
      this->cb_ps = rotr3<ft>(e, -e * qr); // RBDA Eq. 2.28
      this->cb_sp = rotr3<ft>(e.transpose(), qr);
    }

    virtual
    af::const_ref<ft>
    qd_zero() const
    {
      static af::tiny<ft, 6> zeros(0,0,0,0,0,0);
      return af::const_ref<ft>(zeros.begin(), 6);
    }

    virtual
    af::const_ref<ft>
    qdd_zero() const { return qd_zero(); }

    virtual
    af::const_ref<ft, af::mat_grid>
    motion_subspace() const
    {
      return af::const_ref<ft, af::mat_grid>(0, af::mat_grid(6, 6));
    }

    virtual
    boost::optional<vec3<ft> >
    get_linear_velocity(
      af::const_ref<ft> const& qd) const
    {
      SCITBX_ASSERT(qd.size() == 6);
      return boost::optional<vec3<ft> >(vec3<ft>(&qd[3]));
    }

    virtual
    af::small<ft, 6>
    new_linear_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& value) const
    {
      SCITBX_ASSERT(qd.size() == 6);
      SCITBX_ASSERT(value.size() == 3);
      af::small<ft, 6> result(&qd[0], &qd[3]);
      for(unsigned i=0;i<3;i++) result.push_back(value[i]);
      return result;
    }

    virtual
    boost::shared_ptr<joint_base<ft> >
    time_step_position(
      af::const_ref<ft> const& qd,
      ft delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 6);
      vec3<ft> w_body_frame(&qd[0]);
      vec3<ft> v_body_frame(&qd[3]);
      af::tiny<FloatType, 4> new_qe = mat4x3_mul_vec3(
        rbda_eq_4_13(unit_quaternion), w_body_frame);
      new_qe *= delta_t;
      new_qe += qe;
      ft den = af::sum_sq(new_qe);
      if (den == 0) {
        throw std::runtime_error(
          "scitbx::rigid_body::joint_lib::six_dof::time_step_position():"
          " failure computing unit quaternion for angular position:"
          " zero norm.");
      }
      new_qe /= den;
      vec3<ft> new_qr = e.transpose() * v_body_frame;
      new_qr *= delta_t;
      new_qr += qr;
      return boost::shared_ptr<joint_base<ft> >(new six_dof(new_qe, new_qr));
    }

    virtual
    af::small<ft, 6>
    time_step_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& qdd,
      ft delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 6);
      SCITBX_ASSERT(qdd.size() == 6);
      af::small<ft, 6> result(qdd.begin(), qdd.end());
      result *= delta_t;
      for(unsigned i=0;i<6;i++) result[i] += qd[i];
      return result;
    }

    virtual
    af::small<ft, 7>
    tau_as_d_pot_d_q(
      af::const_ref<ft> const& tau) const
    {
      SCITBX_ASSERT(tau.size() == 6);
      af::tiny<ft, 4*4> d = d_unit_quaternion_d_qe_matrix(qe);
      af::tiny<ft, 4*3> c = mat4x4_mul_mat4x3(
        d * ft(4), rbda_eq_4_13(unit_quaternion));
      vec3<ft> n(&tau[0]);
      vec3<ft> f(&tau[3]);
      af::tiny<ft, 4> cn = mat4x3_mul_vec3(c, n);
      vec3<ft> etf = e.transpose() * f;
      af::small<ft, 7> result(cn.begin(), cn.end());
      for(unsigned i=0;i<3;i++) result.push_back(etf[i]);
      return result;
    }

    virtual
    af::small<ft, 7>
    get_q() const
    {
      af::small<ft, 7> result(qe.begin(), qe.end());
      for(unsigned i=0;i<3;i++) result.push_back(qr[i]);
      return result;
    }

    virtual
    boost::shared_ptr<joint_base<ft> >
    new_q(
      af::const_ref<ft> const& q) const
    {
      SCITBX_ASSERT(q.size() == 7);
      return boost::shared_ptr<joint_base<ft> >(new six_dof(
        af::tiny<ft, 4>(&q[0], &q[4]),
        vec3<ft>(&q[4])));
    }
  };

}}} // namespace scitbx::rigid_body::joint_lib

#endif // GUARD
