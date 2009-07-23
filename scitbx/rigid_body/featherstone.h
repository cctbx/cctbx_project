#ifndef SCITBX_RIGID_BODY_FEATHERSTONE_H
#define SCITBX_RIGID_BODY_FEATHERSTONE_H

#include <scitbx/rigid_body/joint_lib.h>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/shared.h>
#include <boost/numeric/conversion/cast.hpp>

namespace scitbx { namespace rigid_body { namespace featherstone {

  //! RBDA Tab. 2.2, p. 23.
  /*! Spatial coordinate transform (rotation around origin).
      Calculates the coordinate transform matrix from A to B coordinates
      for spatial motion vectors, in which frame B is rotated relative to
      frame A.
   */
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  xrot(
    mat3<FloatType> const& e)
  {
    FloatType coeffs[] = {
      e[0], e[1], e[2], 0, 0, 0,
      e[3], e[4], e[5], 0, 0, 0,
      e[6], e[7], e[8], 0, 0, 0,
      0, 0, 0,          e[0], e[1], e[2],
      0, 0, 0,          e[3], e[4], e[5],
      0, 0, 0,          e[6], e[7], e[8]};
    return af::versa_mat_grid(coeffs, 6, 6);
  }

  //! RBDA Tab. 2.2, p. 23.
  /*! Spatial coordinate transform (translation of origin).
      Calculates the coordinate transform matrix from A to B coordinates
      for spatial motion vectors, in which frame B is translated by an
      amount r (3D vector) relative to frame A.
   */
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  xtrans(
    vec3<FloatType> const& r)
  {
    FloatType coeffs[] = {
          1,     0,     0, 0, 0, 0,
          0,     1,     0, 0, 0, 0,
          0,     0,     1, 0, 0, 0,
          0,  r[2], -r[1], 1, 0, 0,
      -r[2],     0,  r[0], 0, 1, 0,
       r[1], -r[0],     0, 0, 0, 1};
    return af::versa_mat_grid(coeffs, 6, 6);
  }

  //! RBDA Eq. 2.28, p. 22.
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  cb_as_spatial_transform(
    rotr3<FloatType> const& cb)
  {
    return af::matrix_multiply(
      xrot(cb.r).const_ref(),
      xtrans(-cb.r.transpose() * cb.t).const_ref());
  }

  //! RBDA Eq. 2.31, p. 25.
  /*! Spatial cross-product operator (motion).
      Calculates the 6x6 matrix such that the expression crm(v)*m is the
      cross product of the spatial motion vectors v and m.
   */
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  crm(
    af::tiny<FloatType, 6> const& v)
  {
    FloatType coeffs[] = {
          0, -v[2],  v[1],     0,     0,     0,
       v[2],     0, -v[0],     0,     0,     0,
      -v[1],  v[0],     0,     0,     0,     0,
          0, -v[5],  v[4],     0, -v[2],  v[1],
       v[5],     0, -v[3],  v[2],     0, -v[0],
      -v[4],  v[3],     0, -v[1],  v[0],     0};
    return af::versa_mat_grid(coeffs, 6, 6);
  }

  //! RBDA Eq. 2.32, p. 25.
  /*! Spatial cross-product operator (force).
      Calculates the 6x6 matrix such that the expression crf(v)*f is the
      cross product of the spatial motion vector v with the spatial force
      vector f.
   */
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  crf(
    af::tiny<FloatType, 6> const& v)
  {
    return -af::matrix_transpose(crm(v).const_ref());
  }

  template <typename FloatType>
  struct body_t
  {
    typedef FloatType ft;

    unsigned number_of_sites;
    ft sum_of_masses;
    boost::shared_ptr<joint_lib::alignment_base<ft> > alignment;
    af::versa<ft, af::mat_grid> i_spatial;
    boost::shared_ptr<joint_lib::joint_base<ft> > joint;
    rotr3<ft> cb_tree;
    int parent;

    virtual
    ~body_t() {}

    virtual
    af::const_ref<ft>
    qd() const = 0;
  };

  //! RBDA Eq. 2.67, p. 35.
  template <typename FloatType>
  FloatType
  kinetic_energy(
    af::const_ref<FloatType, af::mat_grid> const& i_spatial,
    af::tiny<FloatType, 6> const& v_spatial)
  {
    af::tiny<FloatType, 6> iv;
    matrix_mul(iv, i_spatial, v_spatial.const_ref());
    return 0.5 * dot_product(v_spatial, iv);
  }

  //! RBDA Tab. 4.3, p. 87.
  template <typename FloatType=double>
  struct system_model
  {
    typedef FloatType ft;

    af::shared<boost::shared_ptr<body_t<ft> > > bodies;
    af::shared<rotr3<ft> > cb_up_array_;
    af::shared<af::versa<ft, af::mat_grid> > xup_array_;

    unsigned
    bodies_size() const
    {
      return boost::numeric_cast<unsigned>(bodies.size());
    }

    system_model() {}

    /*! \brief Stores bodies and computes body.cb_tree (RBDA Fig. 4.7, p. 74)
        for all bodies.
     */
    system_model(
      af::shared<boost::shared_ptr<body_t<ft> > > const& bodies_)
    :
      bodies(bodies_)
    {
      unsigned nb = bodies_size();
      // XXX TODO move out computation of cb_tree
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        int p = body->parent;
        if (p == -1) {
          body->cb_tree = body->alignment->cb_0b;
        }
        else {
          body->cb_tree = body->alignment->cb_0b * bodies[p]->alignment->cb_b0;
        }
      }
    }

    //! RBDA Example 4.4, p. 80.
    af::shared<rotr3<ft> >
    cb_up_array() const
    {
      if (cb_up_array_.size() == 0) {
        af::shared<rotr3<ft> > cb_up_array_(af::reserve(bodies.size()));
        unsigned nb = bodies_size();
        for(unsigned ib=0;ib<nb;ib++) {
          body_t<ft> const* body = bodies[ib].get();
          cb_up_array_.push_back(body->joint->cb_ps * body->cb_tree);
        }
      }
      return cb_up_array_;
    }

    //! RBDA Example 4.4, p. 80.
    af::shared<af::versa<ft, af::mat_grid> >
    xup_array() const
    {
      if (xup_array_.size() == 0) {
        af::shared<rotr3<ft> > cb_up_array = this->cb_up_array();
        af::shared<af::versa<ft, af::mat_grid> >
          xup_array_(af::reserve(bodies.size()));
        unsigned nb = bodies_size();
        for(unsigned ib=0;ib<nb;ib++) {
          xup_array_.push_back(cb_as_spatial_transform(cb_up_array[ib]));
        }
      }
      return xup_array_;
    }

    //! RBDA Example 4.4, p. 80.
    af::shared<af::tiny<ft, 6> >
    spatial_velocities() const
    {
      af::shared<af::tiny<ft, 6> > result(bodies.size());
      unsigned nb = bodies_size();
      af::shared<rotr3<ft> > cb_up_array = this->cb_up_array();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
        af::const_ref<ft> qd = body->qd();
        af::tiny<ft, 6>& res_ib = result[ib];
        if (s.begin() == 0) {
          SCITBX_ASSERT(qd.size() == 6);
          std::copy(qd.begin(), qd.end(), res_ib.begin()); // vj = qd
        }
        else {
          matrix_mul(res_ib, s, qd); // vj = s * qd
        }
        if (body->parent == -1) {
          // result[ib] = vj, already set
        }
        else {
          // result[ib] = xup_array[i] * result[body->parent] + vj
          rotr3<ft> const& cb_up = cb_up_array[ib];
          af::tiny<ft, 6> const& vp = result[body->parent];
          vec3<ft> r_va = cb_up.r * vec3<ft>(&vp[0]);
          vec3<ft> vl = cb_up.r * vec3<ft>(&vp[3]) + cb_up.t.cross(r_va);
          for(unsigned i=0;i<3;i++) res_ib[i] += r_va[i];
          for(unsigned i=0;i<3;i++) res_ib[i+3] += vl[i];
        }
      }
      return result;
    }

    //! RBDA Eq. 2.67, p. 35.
    ft
    e_kin() const
    {
      ft result = 0;
      af::shared<af::tiny<ft, 6> > sv = spatial_velocities();
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        result += kinetic_energy(body->i_spatial.const_ref(), sv[ib]);
      }
      return result;
    }

    af::shared<af::versa<ft, af::mat_grid> >
    accumulated_spatial_inertia() const
    {
      af::shared<af::versa<ft, af::mat_grid> >
        result(af::reserve(bodies.size()));
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        result.push_back(body->i_spatial);
      }
      af::shared<af::versa<ft, af::mat_grid> > xup_array = this->xup_array();
      for(unsigned ib=nb;ib!=0;) {
        ib--;
        body_t<ft> const* body = bodies[ib].get();
        if (body->parent != -1) {
          result[body->parent] += a_transpose_mul_b_mul_a(
            xup_array[ib].const_ref(),
            result[ib].const_ref());
        }
      }
      return result;
    }

    af::shared<ft>
    qd_e_kin_scales(
      ft e_kin_epsilon=1e-12) const
    {
      af::shared<ft> result(af::reserve(bodies.size()));
      af::shared<af::versa<ft, af::mat_grid> >
        accumulated_spatial_inertia = this->accumulated_spatial_inertia();
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
        unsigned j_dof = body->joint->degrees_of_freedom;
        af::small<ft, 6> qd(j_dof, 0);
        for(unsigned i=0;i<j_dof;i++) {
          qd[i] = 1;
          af::tiny<ft, 6> vj;
          if (s.begin() == 0) {
            SCITBX_ASSERT(j_dof == 6);
            std::copy(qd.begin(), qd.end(), vj.begin()); // vj = qd
          }
          else {
            matrix_mul(vj, s, qd.const_ref()); // vj = s * qd
          }
          qd[i] = 0;
          ft e_kin = kinetic_energy(
            accumulated_spatial_inertia[ib].const_ref(), vj);
          if (e_kin < e_kin_epsilon) {
            result.push_back(1);
          }
          else {
            result.push_back(1 / std::sqrt(e_kin));
          }
        }
      }
      return result;
    }

    //! RBDA Tab. 5.1, p. 96.
    /*! Inverse Dynamics of a kinematic tree via
        Recursive Newton-Euler Algorithm.
        qdd_array is a vector of joint acceleration variables.
        The return value (tau) is a vector of joint force variables.
        f_ext_array specifies external forces acting on the bodies.
        If f_ext_array is None then there are no external forces; otherwise,
        f_ext_array[i] is a spatial force vector giving the force acting on
        body i, expressed in body i coordinates.
        grav_accn is a 6D vector expressing the linear acceleration
        due to gravity.
     */
    af::shared<af::small<ft, 6> >
    inverse_dynamics(
      af::const_ref<af::small<ft, 6> > const& qdd_array,
      af::const_ref<af::tiny<ft, 6> > const& f_ext_array,
      af::const_ref<ft> const& grav_accn) const
    {
      SCITBX_ASSERT(qdd_array.size() == bodies.size());
      SCITBX_ASSERT(
        f_ext_array.size() == (f_ext_array.begin() == 0 ? 0 : bodies.size()));
      SCITBX_ASSERT(grav_accn.size() == (grav_accn.begin() == 0 ? 0 : 6));
      unsigned nb = bodies_size();
      af::shared<af::versa<ft, af::mat_grid> > xup_array = this->xup_array();
      af::shared<af::tiny<ft, 6> > v = spatial_velocities();
      boost::scoped_array<af::tiny<ft, 6> > a(new af::tiny<ft, 6>[nb]);
      boost::scoped_array<af::tiny<ft, 6> > f(new af::tiny<ft, 6>[nb]);
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
        af::const_ref<ft> qd = body->qd();
        af::const_ref<ft> qdd = qdd_array[ib].const_ref();
        af::tiny<ft, 6> vj;
        af::tiny<ft, 6> aj;
        if (s.begin() == 0) {
          SCITBX_ASSERT(qd.size() == 6);
          SCITBX_ASSERT(qdd.size() == 6);
          std::copy(qd.begin(), qd.end(), vj.begin()); // vj = qd
          std::copy(qdd.begin(), qdd.end(), aj.begin()); // aj = qdd
        }
        else {
          matrix_mul(vj, s, qd); // vj = s * qd
          matrix_mul(aj, s, qdd); // aj = s * qdd
        }
        if (body->parent == -1) {
          a[ib] = aj;
          if (grav_accn.begin() != 0) {
            a[ib] -= mat_6xn_mul_vec_n(xup_array[ib].const_ref(), grav_accn);
          }
        }
        else {
          a[ib] =
              mat_6xn_mul_vec_n(
                xup_array[ib].const_ref(),
                a[body->parent].const_ref())
            + aj
            + mat_6xn_mul_vec_n(
                crm(v[ib]).const_ref(),
                vj.const_ref());
        }
        f[ib] =
            mat_6xn_mul_vec_n(
              body->i_spatial.const_ref(),
              a[ib].const_ref())
          + mat_6xn_mul_vec_n(
              crf(v[ib]).const_ref(),
              mat_6xn_mul_vec_n(
                body->i_spatial.const_ref(),
                v[ib].const_ref()).const_ref());
        if (f_ext_array.begin() != 0) {
          f[ib] -= f_ext_array[ib];
        }
      }
      af::shared<af::small<ft, 6> > tau_array(nb);
      for(unsigned ib=nb;ib!=0;) {
        ib--;
        body_t<ft> const* body = bodies[ib].get();
        af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
        if (s.begin() == 0) {
          tau_array[ib] = af::small<ft, 6>(f[ib].begin(), f[ib].end());
        }
        else {
          tau_array[ib] = mat_mxn_transpose_mul_vec_n(s, f[ib].const_ref());
        }
        if (body->parent != -1) {
          f[body->parent] += mat_6x6_transpose_mul_vec6(
            xup_array[ib].const_ref(), f[ib].const_ref());
        }
      }
      return tau_array;
    }

    /*! Simplified version of Inverse Dynamics via Recursive Newton-Euler
        Algorithm, with all qd, qdd zero, but non-zero external forces.
     */
    af::shared<af::small<ft, 6> >
    f_ext_as_tau(
      af::const_ref<af::tiny<ft, 6> > const& f_ext_array) const
    {
      SCITBX_ASSERT(f_ext_array.size() == bodies.size());
      unsigned nb = bodies_size();
      af::shared<af::versa<ft, af::mat_grid> > xup_array = this->xup_array();
      boost::scoped_array<af::tiny<ft, 6> > f(new af::tiny<ft, 6>[nb]);
      for(unsigned ib=0;ib<nb;ib++) {
        f[ib] = -f_ext_array[ib];
      }
      af::shared<af::small<ft, 6> > tau_array(nb);
      for(unsigned ib=nb;ib!=0;) {
        ib--;
        body_t<ft> const* body = bodies[ib].get();
        af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
        if (s.begin() == 0) {
          tau_array[ib] = af::small<ft, 6>(f[ib].begin(), f[ib].end());
        }
        else {
          tau_array[ib] = mat_mxn_transpose_mul_vec_n(s, f[ib].const_ref());
        }
        if (body->parent != -1) {
          f[body->parent] += mat_6x6_transpose_mul_vec6(
            xup_array[ib].const_ref(), f[ib].const_ref());
        }
      }
      return tau_array;
    }

    /*! Gradients of potential energy (defined via f_ext_array) w.r.t.
        positional coordinates q. Uses f_ext_as_tau().
     */
    af::shared<af::small<ft, 7> >
    d_pot_d_q(
      af::const_ref<af::tiny<ft, 6> > const& f_ext_array) const
    {
      af::shared<af::small<ft, 7> > result(af::reserve(bodies.size()));
      af::shared<af::small<ft, 6> >
        tau_array = this->f_ext_as_tau(f_ext_array);
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        result.push_back(bodies[ib]->joint->tau_as_d_pot_d_q(tau_array[ib]));
      }
      return result;
    }

    //! RBDA Tab. 7.1, p. 132.
    /*! Forward Dynamics of a kinematic tree via the Articulated-Body
        Algorithm.
        tau_array is a vector of force variables.
        The return value (qdd_array) is a vector of joint acceleration
        variables.
        f_ext_array specifies external forces acting on the bodies. If
        f_ext_array is None then there are no external forces;
        otherwise, f_ext_array[i] is a spatial force vector giving
        the force acting on body i, expressed in body i coordinates.
        grav_accn is a 6D vector expressing the linear acceleration due
        to gravity.
     */
    af::shared<af::small<ft, 6> >
    forward_dynamics_ab(
      af::const_ref<af::small<ft, 6> > const& tau_array,
      af::const_ref<af::tiny<ft, 6> > const& f_ext_array,
      af::const_ref<ft> const& grav_accn) const
    {
      typedef af::tiny<ft, 6> t6;
      typedef af::small<ft, 6> s6;
      typedef af::versa<ft, af::mat_grid> vmg;
      SCITBX_ASSERT(
        tau_array.size() == (tau_array.begin() == 0 ? 0 : bodies.size()));
      SCITBX_ASSERT(
        f_ext_array.size() == (f_ext_array.begin() == 0 ? 0 : bodies.size()));
      SCITBX_ASSERT(grav_accn.size() == (grav_accn.begin() == 0 ? 0 : 6));
      unsigned nb = bodies_size();
      af::shared<s6> qdd_array(nb); // result
      af::shared<vmg> xup_array = this->xup_array();
      af::shared<t6> v = spatial_velocities();
      boost::scoped_array<t6> c(new t6[nb]);
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
        af::const_ref<ft> qd = body->qd();
        t6 vj;
        if (s.begin() == 0) {
          SCITBX_ASSERT(qd.size() == 6);
          std::copy(qd.begin(), qd.end(), vj.begin()); // vj = qd
        }
        else {
          matrix_mul(vj, s, qd); // vj = s * qd
        }
        if (body->parent == -1) {
          c[ib].fill(0);
        }
        else {
          c[ib] = mat_6xn_mul_vec_n(crm(v[ib]).const_ref(), vj.const_ref());
        }
      }
      boost::scoped_array<vmg> ia(new vmg[nb]);
      for(unsigned ib=0;ib<nb;ib++) {
        ia[ib] = bodies[ib]->i_spatial.deep_copy();
      }
      boost::scoped_array<t6> pa(new t6[nb]);
      for(unsigned ib=0;ib<nb;ib++) {
        pa[ib] = mat_6xn_mul_vec_n(
          crf(v[ib]).const_ref(),
          mat_6xn_mul_vec_n(
            bodies[ib]->i_spatial.const_ref(),
            v[ib].const_ref()).const_ref());
        if (f_ext_array.begin() != 0) {
          pa[ib] -= f_ext_array[ib];
        }
      }
      boost::scoped_array<vmg> u(new vmg[nb]);
      boost::scoped_array<vmg> d_inv(new vmg[nb]);
      boost::scoped_array<s6> u_(new s6[nb]);
      for(unsigned ib=nb;ib!=0;) {
        ib--;
        body_t<ft> const* body = bodies[ib].get();
        af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
        vmg d;
        if (s.begin() == 0) {
          u[ib] = ia[ib].deep_copy(); // XXX deep_copy needed?
          d = u[ib];
          for(unsigned i=0;i<6;i++) {
            u_[ib][i] = -pa[ib][i];
          }
        }
        else {
          u[ib] = af::matrix_multiply(ia[ib].const_ref(), s);
          d = af::matrix_transpose_multiply(s, u[ib].const_ref());
          u_[ib] = -mat_mxn_transpose_mul_vec_n(s, pa[ib].const_ref());
        }
        if (tau_array.begin() != 0) {
          u_[ib] += tau_array[ib];
        }
        // XXX XXX TODO d_inv[ib] = generalized_inverse(d)
        if (body->parent != -1) {
          vmg u_d_inv = af::matrix_multiply(
            u[ib].const_ref(),
            d_inv[ib].const_ref());
          vmg ia_ = ia[ib] - af::matrix_multiply(
            u_d_inv.const_ref(),
            af::matrix_transpose(u[ib].const_ref()).const_ref());
          t6 pa_ = pa[ib]
            + mat_6xn_mul_vec_n(ia_.const_ref(), c[ib].const_ref())
            + mat_6xn_mul_vec_n(u_d_inv.const_ref(), u_[ib].const_ref());
          ia[body->parent] += af::matrix_transpose_multiply(
            xup_array[ib].const_ref(),
            af::matrix_multiply(
              ia_.const_ref(),
              xup_array[ib].const_ref()).const_ref());
          pa[body->parent] += mat_6x6_transpose_mul_vec6(
            xup_array[ib].const_ref(), pa_.const_ref());
        }
      }
      boost::scoped_array<t6> a(new t6[nb]);
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
        a[ib] = c[ib];
        if (body->parent == -1) {
          if (grav_accn.begin() != 0) {
            a[ib] -= mat_6xn_mul_vec_n(xup_array[ib].const_ref(), grav_accn);
          }
        }
        else {
          a[ib] += mat_6xn_mul_vec_n(
            xup_array[ib].const_ref(), a[body->parent].const_ref());
        }
        qdd_array[ib] = mat_mxn_mul_vec_n(
          d_inv[ib].const_ref(),
          (u_[ib] - mat_mxn_transpose_mul_vec_n(
                      u[ib].const_ref(),
                      a[ib].const_ref())).const_ref());
        if (s.begin() == 0) {
          SCITBX_ASSERT(qdd_array.size() == 6);
          for(unsigned i=0;i<6;i++) {
            a[ib][i] += qdd_array[ib][i];
          }
        }
        else {
          a[ib] += mat_6xn_mul_vec_n(s, qdd_array[ib].const_ref());
        }
      }
      return qdd_array;
    }
  };

}}} // namespace scitbx::rigid_body::featherstone

#endif // GUARD
