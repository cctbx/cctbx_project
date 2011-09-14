#ifndef SCITBX_RIGID_BODY_JOINT_LIB_H
#define SCITBX_RIGID_BODY_JOINT_LIB_H

#include <scitbx/rigid_body/joint_t.h>
#include <scitbx/rigid_body/matrix_helpers.h>
#include <scitbx/math/r3_rotation.h>
#include <boost/numeric/conversion/cast.hpp>

namespace scitbx { namespace rigid_body {

//! See essence/joint_lib.py
namespace joint_lib {

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

  //! See essence/joint_lib.py
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

  //! See code.
  template <typename FloatType=double>
  struct zero_dof_alignment : alignment_t<FloatType>
  {
    typedef FloatType ft;

    zero_dof_alignment()
    :
      alignment_t<ft>(rotr3<ft>::identity(), rotr3<ft>::identity())
    {}
  };

  //! Zero degree-of-freedom joint model.
  template <typename FloatType=double>
  struct zero_dof : joint_t<FloatType>
  {
    typedef FloatType ft;

    zero_dof()
    :
      joint_t<ft>(0, 0)
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
      // empty list does not work with some compilers (e.g. Visual C++ 8.0)
      static ft const coeffs[] = { -1e20 };
      return af::const_ref<ft, af::mat_grid>(coeffs, af::mat_grid(6, 0));
    }

    virtual
    boost::optional<vec3<ft> >
    get_linear_velocity(
      af::const_ref<ft> const& qd) const
    {
      SCITBX_ASSERT(qd.size() == 0);
      return boost::optional<vec3<ft> >();
    }

    virtual
    af::small<ft, 6>
    new_linear_velocity(
      af::const_ref<ft> const& qd,
      vec3<ft> const& value) const
    {
      SCITBX_ASSERT(qd.size() == 0);
      return af::small<ft, 6>(0);
    }

    virtual
    shared_ptr<joint_t<ft> >
    time_step_position(
      af::const_ref<ft> const& qd,
      ft const& delta_t) const
    {
      return shared_ptr<joint_t<ft> >(new zero_dof());
    }

    virtual
    af::small<ft, 6>
    time_step_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& qdd,
      ft const& delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 0);
      SCITBX_ASSERT(qdd.size() == 0);
      return af::small<ft, 6>(0);
    }

    virtual
    af::small<ft, 7>
    tau_as_d_e_pot_d_q(
      af::small<ft, 6> const& tau) const
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
    shared_ptr<joint_t<ft> >
    new_q(
      af::const_ref<ft> const& q) const
    {
      return shared_ptr<joint_t<ft> >(new zero_dof());
    }
  };

  //! See code.
  template <typename FloatType=double>
  struct six_dof_alignment : alignment_t<FloatType>
  {
    typedef FloatType ft;

    six_dof_alignment(
      vec3<ft> const& center_of_mass)
    :
      alignment_t<ft>(
        rotr3<ft>(mat3<ft>(1,1,1), -center_of_mass),
        rotr3<ft>(mat3<ft>(1,1,1), center_of_mass))
    {}
  };

  //! Six degree-of-freedom joint model (see Featherstone RBDA 2007).
  template <typename FloatType=double>
  struct six_dof : joint_t<FloatType>
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
      joint_t<ft>(6, 7),
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
      return zeros.const_ref();
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
      vec3<ft> const& value) const
    {
      SCITBX_ASSERT(qd.size() == 6);
      af::small<ft, 6> result(&qd[0], &qd[3]);
      for(unsigned i=0;i<3;i++) result.push_back(value[i]);
      return result;
    }

    virtual
    shared_ptr<joint_t<ft> >
    time_step_position(
      af::const_ref<ft> const& qd,
      ft const& delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 6);
      vec3<ft> w_body_frame(&qd[0]);
      vec3<ft> v_body_frame(&qd[3]);
      af::tiny<FloatType, 4> new_qe = mat4x3_mul_vec3(
        rbda_eq_4_13(unit_quaternion), w_body_frame);
      new_qe *= delta_t;
      new_qe += qe;
      ft den = std::sqrt(af::sum_sq(new_qe));
      if (den == 0) {
        throw std::runtime_error(
          "scitbx::rigid_body::joint_lib::six_dof::time_step_position():"
          " failure computing unit quaternion for angular position:"
          " zero norm.");
      }
      new_qe /= den;
      vec3<ft> new_qr = v_body_frame * e;
      new_qr *= delta_t;
      new_qr += qr;
      return shared_ptr<joint_t<ft> >(new six_dof(new_qe, new_qr));
    }

    virtual
    af::small<ft, 6>
    time_step_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& qdd,
      ft const& delta_t) const
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
    tau_as_d_e_pot_d_q(
      af::small<ft, 6> const& tau) const
    {
      SCITBX_ASSERT(tau.size() == 6);
      af::tiny<ft, 4*4> d = d_unit_quaternion_d_qe_matrix(qe);
      af::tiny<ft, 4*3> c = mat4x4_mul_mat4x3(
        d * ft(4), rbda_eq_4_13(unit_quaternion));
      vec3<ft> n(&tau[0]);
      vec3<ft> f(&tau[3]);
      af::tiny<ft, 4> cn = mat4x3_mul_vec3(c, n);
      vec3<ft> etf = f * e;
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
    shared_ptr<joint_t<ft> >
    new_q(
      af::const_ref<ft> const& q) const
    {
      SCITBX_ASSERT(q.size() == 7);
      return shared_ptr<joint_t<ft> >(new six_dof(
        af::tiny<ft, 4>(&q[0], &q[4]),
        vec3<ft>(&q[4])));
    }
  };

  /*! \brief Simplification of aja (featherstone.h) for
      six_dof_alignment and six_dof above.
   */
  template <typename FloatType>
  rotr3<FloatType>
  six_dof_aja_simplified(
    vec3<FloatType> const& center_of_mass,
    af::const_ref<FloatType> const& q)
  {
    SCITBX_ASSERT(q.size() == 7);
    typedef FloatType ft;
    af::tiny<ft, 4> qe(&q[0], &q[4]);
    vec3<ft> qr(&q[4]);
    af::tiny<ft, 4> unit_quaternion(vec4_normalize(qe));
    mat3<ft> et(rbda_eq_4_12(unit_quaternion).transpose());
    return rotr3<ft>(et, center_of_mass + qr - et * center_of_mass);
  }

  //! See code.
  template <typename FloatType=double>
  struct spherical_alignment : alignment_t<FloatType>
  {
    typedef FloatType ft;

    spherical_alignment(
      vec3<ft> const& pivot)
    :
      alignment_t<ft>(
        rotr3<ft>(mat3<ft>(1,1,1), -pivot),
        rotr3<ft>(mat3<ft>(1,1,1), pivot))
    {}
  };

  /*! \brief Spherical (three degrees of freedom) joint model
      (see Featherstone RBDA 2007).
   */
  template <typename FloatType=double>
  struct spherical : joint_t<FloatType>
  {
    typedef FloatType ft;

    af::tiny<ft, 4> qe;
    af::tiny<ft, 4> unit_quaternion;

    spherical(
      af::tiny<ft, 4> const& qe_)
    :
      joint_t<ft>(3, 4),
      qe(qe_),
      unit_quaternion(vec4_normalize(qe_)) // RBDA, bottom of p. 86
    {
      mat3<ft> e(rbda_eq_4_12(unit_quaternion));
      this->cb_ps = rotr3<ft>(e, vec3<ft>(0,0,0));
      this->cb_sp = rotr3<ft>(e.transpose(), vec3<ft>(0,0,0));
    }

    virtual
    af::const_ref<ft>
    qd_zero() const
    {
      static af::tiny<ft, 3> zeros(0,0,0);
      return zeros.const_ref();
    }

    virtual
    af::const_ref<ft>
    qdd_zero() const { return qd_zero(); }

    virtual
    af::const_ref<ft, af::mat_grid>
    motion_subspace() const
    {
      static ft const coeffs[] = {
        1,0,0,
        0,1,0,
        0,0,1,
        0,0,0,
        0,0,0,
        0,0,0};
      return af::const_ref<ft, af::mat_grid>(coeffs, af::mat_grid(6, 3));
    }

    virtual
    boost::optional<vec3<ft> >
    get_linear_velocity(
      af::const_ref<ft> const& qd) const
    {
      SCITBX_ASSERT(qd.size() == 3);
      return boost::optional<vec3<ft> >();
    }

    virtual
    af::small<ft, 6>
    new_linear_velocity(
      af::const_ref<ft> const& qd,
      vec3<ft> const& value) const
    {
      SCITBX_ASSERT(qd.size() == 0);
      return af::small<ft, 6>(0);
    }

    virtual
    shared_ptr<joint_t<ft> >
    time_step_position(
      af::const_ref<ft> const& qd,
      ft const& delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 3);
      vec3<ft> w_body_frame(&qd[0]);
      af::tiny<FloatType, 4> new_qe = mat4x3_mul_vec3(
        rbda_eq_4_13(unit_quaternion), w_body_frame);
      new_qe *= delta_t;
      new_qe += qe;
      ft den = std::sqrt(af::sum_sq(new_qe));
      if (den == 0) {
        throw std::runtime_error(
          "scitbx::rigid_body::joint_lib::spherical::time_step_position():"
          " failure computing unit quaternion for angular position:"
          " zero norm.");
      }
      new_qe /= den;
      return shared_ptr<joint_t<ft> >(new spherical(new_qe));
    }

    virtual
    af::small<ft, 6>
    time_step_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& qdd,
      ft const& delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 3);
      SCITBX_ASSERT(qdd.size() == 3);
      af::small<ft, 6> result(qdd.begin(), qdd.end());
      result *= delta_t;
      for(unsigned i=0;i<3;i++) result[i] += qd[i];
      return result;
    }

    virtual
    af::small<ft, 7>
    tau_as_d_e_pot_d_q(
      af::small<ft, 6> const& tau) const
    {
      SCITBX_ASSERT(tau.size() == 3);
      af::tiny<ft, 4*4> d = d_unit_quaternion_d_qe_matrix(qe);
      af::tiny<ft, 4*3> c = mat4x4_mul_mat4x3(
        d * ft(4), rbda_eq_4_13(unit_quaternion));
      vec3<ft> n(&tau[0]);
      af::tiny<ft, 4> cn = mat4x3_mul_vec3(c, n);
      return af::small<ft, 7>(cn.begin(), cn.end());
    }

    virtual
    af::small<ft, 7>
    get_q() const
    {
      return af::small<ft, 7>(qe.begin(), qe.end());
    }

    virtual
    shared_ptr<joint_t<ft> >
    new_q(
      af::const_ref<ft> const& q) const
    {
      SCITBX_ASSERT(q.size() == 4);
      return shared_ptr<joint_t<ft> >(new spherical(
        af::tiny<ft, 4>(q.begin(), q.end())));
    }
  };

  //! See code.
  template <typename FloatType=double>
  struct revolute_alignment : alignment_t<FloatType>
  {
    typedef FloatType ft;

    revolute_alignment(
      vec3<ft> const& pivot,
      vec3<ft> const& normal)
    {
      mat3<ft> r = math::r3_rotation::vector_to_001(normal);
      this->cb_0b = rotr3<ft>(r, -r * pivot);
      this->cb_b0 = rotr3<ft>(r.transpose(), pivot);
    }
  };

  /*! \brief Revolute (one degree of freedom) joint model
      (see Featherstone RBDA 2007).
   */
  template <typename FloatType=double>
  struct revolute : joint_t<FloatType>
  {
    typedef FloatType ft;

    ft qe;

    revolute(
      af::tiny<ft, 1> const& qe_)
    :
      joint_t<ft>(1, 1),
      qe(qe_[0])
    {
      ft c = boost::numeric_cast<ft>(std::cos(qe));
      ft s = boost::numeric_cast<ft>(std::sin(qe));
      mat3<ft> e(c, s, 0, -s, c, 0, 0, 0, 1); // RBDA Tab. 2.2
      this->cb_ps = rotr3<ft>(e, vec3<ft>(0,0,0));
      this->cb_sp = rotr3<ft>(e.transpose(), vec3<ft>(0,0,0));
    }

    virtual
    af::const_ref<ft>
    qd_zero() const
    {
      static af::tiny<ft, 1> zeros(0);
      return zeros.const_ref();
    }

    virtual
    af::const_ref<ft>
    qdd_zero() const { return qd_zero(); }

    virtual
    af::const_ref<ft, af::mat_grid>
    motion_subspace() const
    {
      static ft const coeffs[] = {0, 0, 1, 0, 0, 0};
      return af::const_ref<ft, af::mat_grid>(coeffs, af::mat_grid(6, 1));
    }

    virtual
    boost::optional<vec3<ft> >
    get_linear_velocity(
      af::const_ref<ft> const& qd) const
    {
      SCITBX_ASSERT(qd.size() == 1);
      return boost::optional<vec3<ft> >();
    }

    virtual
    af::small<ft, 6>
    new_linear_velocity(
      af::const_ref<ft> const& qd,
      vec3<ft> const& value) const
    {
      SCITBX_ASSERT(qd.size() == 1);
      return af::small<ft, 6>(0);
    }

    virtual
    shared_ptr<joint_t<ft> >
    time_step_position(
      af::const_ref<ft> const& qd,
      ft const& delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 1);
      af::tiny<ft, 1> new_qe(qe + qd[0] * delta_t);
      return shared_ptr<joint_t<ft> >(new revolute(new_qe));
    }

    virtual
    af::small<ft, 6>
    time_step_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& qdd,
      ft const& delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 1);
      SCITBX_ASSERT(qdd.size() == 1);
      af::small<ft, 6> result(qdd.begin(), qdd.end());
      result[0] *= delta_t;
      result[0] += qd[0];
      return result;
    }

    virtual
    af::small<ft, 7>
    tau_as_d_e_pot_d_q(
      af::small<ft, 6> const& tau) const
    {
      SCITBX_ASSERT(tau.size() == 1);
      return af::small<ft, 7>(tau.begin(), tau.end());
    }

    virtual
    af::small<ft, 7>
    get_q() const
    {
      return af::small<ft, 7>(1U, qe);
    }

    virtual
    shared_ptr<joint_t<ft> >
    new_q(
      af::const_ref<ft> const& q) const
    {
      SCITBX_ASSERT(q.size() == 1);
      return shared_ptr<joint_t<ft> >(new revolute(
        af::tiny<ft, 1>(q.begin(), q.end())));
    }
  };

  //! See code.
  template <typename FloatType=double>
  struct translational_alignment : six_dof_alignment<FloatType>
  {
    typedef FloatType ft;

    translational_alignment(
      vec3<ft> const& center_of_mass)
    :
      six_dof_alignment<ft>(center_of_mass)
    {}
  };

  /*! \brief Translational (three degrees of freedom) joint model
      (see Featherstone RBDA 2007).
   */
  template <typename FloatType=double>
  struct translational : joint_t<FloatType>
  {
    typedef FloatType ft;

    vec3<ft> qr;

    translational(
      vec3<ft> const& qr_)
    :
      joint_t<ft>(3, 3),
      qr(qr_)
    {
      this->cb_ps = rotr3<ft>(mat3<ft>(1,1,1), -qr);
      this->cb_sp = rotr3<ft>(mat3<ft>(1,1,1), qr);
    }

    virtual
    af::const_ref<ft>
    qd_zero() const
    {
      static af::tiny<ft, 3> zeros(0,0,0);
      return zeros.const_ref();
    }

    virtual
    af::const_ref<ft>
    qdd_zero() const { return qd_zero(); }

    virtual
    af::const_ref<ft, af::mat_grid>
    motion_subspace() const
    {
      static ft const coeffs[] = {
        0,0,0,
        0,0,0,
        0,0,0,
        1,0,0,
        0,1,0,
        0,0,1};
      return af::const_ref<ft, af::mat_grid>(coeffs, af::mat_grid(6, 3));
    }

    virtual
    boost::optional<vec3<ft> >
    get_linear_velocity(
      af::const_ref<ft> const& qd) const
    {
      SCITBX_ASSERT(qd.size() == 3);
      return boost::optional<vec3<ft> >(&qd[0]);
    }

    virtual
    af::small<ft, 6>
    new_linear_velocity(
      af::const_ref<ft> const& qd,
      vec3<ft> const& value) const
    {
      SCITBX_ASSERT(qd.size() == 3);
      return af::small<ft, 6>(&value[0], &value[3]);
    }

    virtual
    shared_ptr<joint_t<ft> >
    time_step_position(
      af::const_ref<ft> const& qd,
      ft const& delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 3);
      vec3<ft> new_qr(&qd[0]);
      new_qr *= delta_t;
      new_qr += qr;
      return shared_ptr<joint_t<ft> >(new translational(new_qr));
    }

    virtual
    af::small<ft, 6>
    time_step_velocity(
      af::const_ref<ft> const& qd,
      af::const_ref<ft> const& qdd,
      ft const& delta_t) const
    {
      SCITBX_ASSERT(qd.size() == 3);
      SCITBX_ASSERT(qdd.size() == 3);
      af::small<ft, 6> result(qdd.begin(), qdd.end());
      result *= delta_t;
      for(unsigned i=0;i<3;i++) result[i] += qd[i];
      return result;
    }

    virtual
    af::small<ft, 7>
    tau_as_d_e_pot_d_q(
      af::small<ft, 6> const& tau) const
    {
      SCITBX_ASSERT(tau.size() == 3);
      return af::small<ft, 7>(tau.begin(), tau.end());
    }

    virtual
    af::small<ft, 7>
    get_q() const
    {
      return af::small<ft, 7>(qr.begin(), qr.end());
    }

    virtual
    shared_ptr<joint_t<ft> >
    new_q(
      af::const_ref<ft> const& q) const
    {
      SCITBX_ASSERT(q.size() == 3);
      return shared_ptr<joint_t<ft> >(new translational(
        vec3<ft>(&q[0])));
    }
  };

}}} // namespace scitbx::rigid_body::joint_lib

#endif // GUARD
