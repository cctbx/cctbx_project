#ifndef SCITBX_RIGID_BODY_FEATHERSTONE_H
#define SCITBX_RIGID_BODY_FEATHERSTONE_H

#include <scitbx/rigid_body/joint_lib.h>
#include <scitbx/array_family/versa_matrix.h>
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
  };

}}} // namespace scitbx::rigid_body::featherstone

#endif // GUARD
