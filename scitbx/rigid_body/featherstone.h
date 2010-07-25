#ifndef SCITBX_RIGID_BODY_FEATHERSTONE_H
#define SCITBX_RIGID_BODY_FEATHERSTONE_H

#include <scitbx/rigid_body/spatial_lib.h>
#include <scitbx/rigid_body/array_packing.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/shared_algebra.h>
#include <tbxx/optional_copy.hpp>
#include <boost/noncopyable.hpp>

namespace scitbx { namespace rigid_body {

//! See essence/featherstone.py
namespace featherstone {

  //! Helper.
  template <typename FloatType>
  struct random_gauss_adaptor
  {
    virtual
    ~random_gauss_adaptor() {}

    virtual
    FloatType
    operator()(
      FloatType const& mu,
      FloatType const& sigma) = 0;
  };

  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  generalized_inverse(
    af::const_ref<FloatType, af::mat_grid> const& m)
  {
    // assumption to achieve stability: order of magnitude of masses is around 1
    return matrix::packed_u_as_symmetric(
      scitbx::matrix::eigensystem::real_symmetric<FloatType>(
        m,
        /*relative_epsilon*/ 1e-6,
        /*absolute_epsilon*/ 1e-6)
          .generalized_inverse_as_packed_u().const_ref());
  }

  //! RBDA Tab. 4.3, p. 87.
  template <typename FloatType=double>
  struct system_model : boost::noncopyable
  {
    typedef FloatType ft;

    af::shared<shared_ptr<body_t<ft> > > bodies;
    unsigned number_of_trees;
    unsigned degrees_of_freedom;
    unsigned q_packed_size;

    protected:
      boost::optional<af::shared<rotr3<ft> > > aja_array_;
      boost::optional<af::shared<mat3<ft> > > jar_array_;
      boost::optional<af::shared<rotr3<ft> > > cb_up_array_;
      boost::optional<af::shared<af::versa<ft, af::mat_grid> > > xup_array_;
      boost::optional<af::shared<af::tiny<ft, 6> > > spatial_velocities_;
      boost::optional<ft> e_kin_;
    public:

    system_model() {}

    system_model(
      af::shared<shared_ptr<body_t<ft> > > const& bodies_)
    :
      bodies(bodies_),
      number_of_trees(0),
      degrees_of_freedom(0),
      q_packed_size(0)
    {
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        if (body->parent == -1) number_of_trees++;
        degrees_of_freedom += body->joint->degrees_of_freedom;
        q_packed_size += body->joint->q_size;
      }
    }

    virtual
    ~system_model() {}

    unsigned
    bodies_size() const
    {
      return boost::numeric_cast<unsigned>(bodies.size());
    }

    //! For reporting in Python.
    af::shared<std::size_t>
    degrees_of_freedom_each_joint() const
    {
#define SCITBX_LOC(ATTR) \
      unsigned nb = bodies_size(); \
      af::shared<std::size_t> result((af::reserve(nb))); \
      for(unsigned ib=0;ib<nb;ib++) { \
        body_t<ft> const* body = bodies[ib].get(); \
        result.push_back(boost::numeric_cast<std::size_t>(body->joint->ATTR));\
      } \
      return result;
      SCITBX_LOC(degrees_of_freedom)
    }

    //! For reporting in Python.
    af::shared<std::size_t>
    q_size_each_joint() const
    {
      SCITBX_LOC(q_size)
#undef SCITBX_LOC
    }

    virtual
    void
    flag_positions_as_changed()
    {
      aja_array_.reset();
      jar_array_.reset();
      cb_up_array_.reset();
      xup_array_.reset();
      flag_velocities_as_changed();
    }

    virtual
    void
    flag_velocities_as_changed()
    {
      spatial_velocities_.reset();
      e_kin_.reset();
    }

    af::shared<std::size_t>
    root_indices() const
    {
      af::shared<std::size_t> result((af::reserve(number_of_trees)));
      std::size_t nb = bodies.size();
      for(std::size_t ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        if (body->parent == -1) {
          result.push_back(ib);
        }
      }
      SCITBX_ASSERT(result.size() == number_of_trees);
      return result;
    }

    af::shared<ft>
    pack_q()
    {
      af::shared<ft> result;
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        af::small<ft, 7> q = bodies[ib]->joint->get_q();
        result.extend(q.begin(), q.end());
      }
      SCITBX_ASSERT(result.size() == q_packed_size);
      return result;
    }

    void
    unpack_q(
      af::const_ref<ft> const& q_packed)
    {
      SCITBX_ASSERT(q_packed.size() == q_packed_size);
      unsigned i = 0;
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        unsigned n = body->joint->q_size;
        body->joint = body->joint->new_q(af::const_ref<ft>(&q_packed[i], n));
        i += n;
      }
      SCITBX_ASSERT(i == q_packed_size);
      flag_positions_as_changed();
    }

    af::shared<ft>
    pack_qd()
    {
      af::shared<ft> result;
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        af::const_ref<ft> qd = bodies[ib]->qd();
        result.extend(qd.begin(), qd.end());
      }
      SCITBX_ASSERT(result.size() == degrees_of_freedom);
      return result;
    }

    void
    unpack_qd(
      af::const_ref<ft> const& qd_packed)
    {
      SCITBX_ASSERT(qd_packed.size() == degrees_of_freedom);
      unsigned i = 0;
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        unsigned n = body->joint->degrees_of_freedom;
        body->set_qd(af::small<ft, 6>(af::adapt(
          af::const_ref<ft>(&qd_packed[i], n))));
        i += n;
      }
      SCITBX_ASSERT(i == degrees_of_freedom);
      flag_velocities_as_changed();
    }

    af::shared<af::tiny<std::size_t, 2> >
    number_of_sites_in_each_tree() const
    {
      af::shared<af::tiny<std::size_t, 2> >
        result((af::reserve(number_of_trees)));
      unsigned nb = bodies_size();
      boost::scoped_array<unsigned> accu(new unsigned[nb]);
      std::fill_n(accu.get(), nb, unsigned(0));
      for(unsigned ib=nb;ib!=0;) {
        ib--;
        body_t<ft> const* body = bodies[ib].get();
        accu[ib] += body->number_of_sites;
        if (body->parent == -1) {
          result.push_back(af::tiny<std::size_t, 2>(ib, accu[ib]));
        }
        else {
          accu[body->parent] += accu[ib];
        }
      }
      SCITBX_ASSERT(result.size() == number_of_trees);
      return result;
    }

    af::shared<std::pair<int, double> >
    sum_of_masses_in_each_tree() const
    {
      af::shared<std::pair<int, double> >
        result((af::reserve(number_of_trees)));
      unsigned nb = bodies_size();
      boost::scoped_array<ft> accu(new ft[nb]);
      std::fill_n(accu.get(), nb, ft(0));
      for(unsigned ib=nb;ib!=0;) {
        ib--;
        body_t<ft> const* body = bodies[ib].get();
        accu[ib] += body->sum_of_masses;
        if (body->parent == -1) {
          result.push_back(std::pair<int, double>(
            boost::numeric_cast<int>(ib),
            boost::numeric_cast<double>(accu[ib])));
        }
        else {
          accu[body->parent] += accu[ib];
        }
      }
      SCITBX_ASSERT(result.size() == number_of_trees);
      return result;
    }

    boost::optional<vec3<ft> >
    mean_linear_velocity(
      af::const_ref<af::tiny<std::size_t, 2> >
        number_of_sites_in_each_tree) const
    {
      vec3<ft> sum_v(0,0,0);
      unsigned sum_n = 0;
#define SCITBX_LOC \
      tbxx::optional_container< \
        af::shared<af::tiny<std::size_t, 2> > > \
          nosiet; \
      if (number_of_sites_in_each_tree.begin() == 0) { \
        nosiet = this->number_of_sites_in_each_tree(); \
        number_of_sites_in_each_tree = nosiet->const_ref(); \
      } \
      SCITBX_ASSERT(number_of_sites_in_each_tree.size() == number_of_trees); \
      std::size_t nb = bodies.size(); \
      for( \
        af::tiny<std::size_t, 2> const* \
          nosiet_it=number_of_sites_in_each_tree.begin(); \
        nosiet_it!=number_of_sites_in_each_tree.end(); \
        nosiet_it++) \
      { \
        std::size_t ib = (*nosiet_it)[0]; \
        SCITBX_ASSERT(ib < nb);
SCITBX_LOC
        body_t<ft> const* body = bodies[ib].get();
        boost::optional<vec3<ft> >
          v = body->joint->get_linear_velocity(body->qd());
        if (!v) continue;
        unsigned n = boost::numeric_cast<unsigned>((*nosiet_it)[1]);
        sum_v += (*v) * boost::numeric_cast<ft>(n);
        sum_n += n;
      }
      if (sum_n == 0) {
        return boost::optional<vec3<ft> >();
      }
      return boost::optional<vec3<ft> >(
        sum_v / boost::numeric_cast<ft>(sum_n));
    }

    void
    subtract_from_linear_velocities(
      af::const_ref<af::tiny<std::size_t, 2> >
        number_of_sites_in_each_tree,
      vec3<ft> const& value)
    {
SCITBX_LOC // {
#undef SCITBX_LOC
        body_t<ft>* body = bodies[ib].get();
        boost::optional<vec3<ft> >
          v = body->joint->get_linear_velocity(body->qd());
        if (!v) continue;
        body->set_qd(
          body->joint->new_linear_velocity(body->qd(), (*v)-value));
      }
    }

    //! Not available in Python.
    af::shared<rotr3<ft> > const&
    aja_array()
    {
      if (!aja_array_) {
        unsigned nb = bodies_size();
        aja_array_ = af::shared<rotr3<ft> >(af::reserve(nb));
        for(unsigned ib=0;ib<nb;ib++) {
          body_t<ft> const* body = bodies[ib].get();
          rotr3<ft>
            aja = body->alignment->cb_b0
                * body->joint->cb_sp
                * body->alignment->cb_0b;
          if (body->parent != -1) {
            aja = (*aja_array_)[body->parent] * aja;
          }
          aja_array_->push_back(aja);
        }
      }
      return *aja_array_;
    }

    //! Not available in Python.
    af::shared<mat3<ft> > const&
    jar_array()
    {
      if (!jar_array_) {
        aja_array();
        unsigned nb = bodies_size();
        jar_array_ = af::shared<mat3<ft> >(af::reserve(nb));
        for(unsigned ib=0;ib<nb;ib++) {
          body_t<ft> const* body = bodies[ib].get();
          mat3<ft> jar = body->joint->cb_ps.r * body->alignment->cb_0b.r;
          if (body->parent != -1) {
            jar = mul_transpose(jar, (*aja_array_)[body->parent].r);
          }
          jar_array_->push_back(jar);
        }
      }
      return *jar_array_;
    }

    //! RBDA Example 4.4, p. 80.
    /*! Not available in Python.
     */
    af::shared<rotr3<ft> > const&
    cb_up_array()
    {
      if (!cb_up_array_) {
        unsigned nb = bodies_size();
        cb_up_array_ = af::shared<rotr3<ft> >(af::reserve(nb));
        for(unsigned ib=0;ib<nb;ib++) {
          body_t<ft> const* body = bodies[ib].get();
          cb_up_array_->push_back(body->joint->cb_ps * body->cb_tree);
        }
      }
      return *cb_up_array_;
    }

    //! RBDA Example 4.4, p. 80.
    /*! Not available in Python.
     */
    af::shared<af::versa<ft, af::mat_grid> > const&
    xup_array()
    {
      if (!xup_array_) {
        cb_up_array();
        unsigned nb = bodies_size();
        xup_array_ = af::shared<af::versa<ft, af::mat_grid> >(
          af::reserve(nb));
        for(unsigned ib=0;ib<nb;ib++) {
          xup_array_->push_back(
            spatial_lib::cb_as_spatial_transform((*cb_up_array_)[ib]));
        }
      }
      return *xup_array_;
    }

    //! RBDA Example 4.4, p. 80.
    /*! Not available in Python.
     */
    af::shared<af::tiny<ft, 6> > const&
    spatial_velocities()
    {
      if (!spatial_velocities_) {
        unsigned nb = bodies_size();
        spatial_velocities_ = af::shared<af::tiny<ft, 6> >((nb));
        af::shared<rotr3<ft> > cb_up_array = this->cb_up_array();
        for(unsigned ib=0;ib<nb;ib++) {
          body_t<ft> const* body = bodies[ib].get();
          af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
          af::const_ref<ft> qd = body->qd();
          af::tiny<ft, 6>& res_ib = (*spatial_velocities_)[ib];
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
            af::tiny<ft, 6> const& vp = (*spatial_velocities_)[body->parent];
            vec3<ft> r_va = cb_up.r * vec3<ft>(&vp[0]);
            vec3<ft> vl = cb_up.r * vec3<ft>(&vp[3]) + cb_up.t.cross(r_va);
            res_ib += spatial_lib::as_tiny_6(r_va, vl);
          }
        }
      }
      return *spatial_velocities_;
    }

    //! RBDA Eq. 2.67, p. 35.
    ft const&
    e_kin()
    {
      if (!e_kin_) {
        ft result = 0;
        af::shared<af::tiny<ft, 6> > sv = spatial_velocities();
        unsigned nb = bodies_size();
        for(unsigned ib=0;ib<nb;ib++) {
          body_t<ft> const* body = bodies[ib].get();
          result += spatial_lib::kinetic_energy(
            body->i_spatial.const_ref(), sv[ib]);
        }
        e_kin_ = result;
      }
      return *e_kin_;
    }

    void
    reset_e_kin(
      ft const& e_kin_target,
      ft const& e_kin_epsilon=1e-12)
    {
      SCITBX_ASSERT(e_kin_target >= 0);
      SCITBX_ASSERT(e_kin_epsilon > 0);
      ft e_kin = this->e_kin();
      if (e_kin >= e_kin_epsilon) {
        ft factor = std::sqrt(e_kin_target / e_kin);
        unsigned nb = bodies_size();
        for(unsigned ib=0;ib<nb;ib++) {
          body_t<ft>* body = bodies[ib].get();
          af::ref<ft> body_qd = body->qd();
          for(std::size_t i=0;i<body_qd.size();i++) {
            body_qd[i] *= factor;
          }
        }
      }
      flag_velocities_as_changed();
    }

    void
    assign_zero_velocities()
    {
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        af::ref<ft> body_qd = body->qd();
        af::const_ref<ft> joint_qd_zero = body->joint->qd_zero();
        SCITBX_ASSERT(joint_qd_zero.size() == body_qd.size());
        std::copy(
          joint_qd_zero.begin(),
          joint_qd_zero.end(), body_qd.begin());
      }
      flag_velocities_as_changed();
    }

    af::shared<af::versa<ft, af::mat_grid> >
    accumulated_spatial_inertia()
    {
      af::shared<af::versa<ft, af::mat_grid> >
        result(af::reserve(bodies.size()));
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        result.push_back(body->i_spatial.deep_copy());
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
      ft e_kin_epsilon=1e-12)
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
          ft e_kin = spatial_lib::kinetic_energy(
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

    boost::optional<af::shared<ft> >
    assign_random_velocities(
      random_gauss_adaptor<ft>& random_gauss,
      boost::optional<ft> const& e_kin_target=boost::optional<ft>(),
      ft const& e_kin_epsilon=1e-12)
    {
      ft work_e_kin_target;
      if (!e_kin_target) {
        work_e_kin_target = 1;
      }
      else if (*e_kin_target == 0) {
        assign_zero_velocities();
        return boost::optional<af::shared<ft> >();
      }
      else {
        SCITBX_ASSERT(*e_kin_target >= 0);
        work_e_kin_target = *e_kin_target;
      }
      af::shared<ft> qd_e_kin_scales = this->qd_e_kin_scales(e_kin_epsilon);
      if (degrees_of_freedom != 0) {
        qd_e_kin_scales *= boost::numeric_cast<ft>(
          std::sqrt(
              work_e_kin_target
            / boost::numeric_cast<ft>(degrees_of_freedom)));
      }
      unsigned i_qd = 0;
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        af::small<ft, 6> qd_new(af::adapt(body->joint->qd_zero()));
        unsigned n = boost::numeric_cast<unsigned>(qd_new.size());
        for(unsigned i=0;i<n;i++,i_qd++) {
          qd_new[i] += random_gauss(/*mu*/ 0, /*sigma*/ qd_e_kin_scales[i_qd]);
        }
        body->set_qd(qd_new);
      }
      SCITBX_ASSERT(i_qd == degrees_of_freedom);
      flag_velocities_as_changed();
      if (e_kin_target) {
        reset_e_kin(*e_kin_target, e_kin_epsilon);
      }
      return boost::optional<af::shared<ft> >(qd_e_kin_scales);
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
      af::const_ref<ft> const& grav_accn)
    {
      SCITBX_ASSERT(
        qdd_array.size() == (qdd_array.begin() == 0 ? 0 : bodies.size()));
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
        af::tiny<ft, 6> vj;
        af::tiny<ft, 6> aj;
        if (s.begin() == 0) {
          SCITBX_ASSERT(qd.size() == 6);
          std::copy(qd.begin(), qd.end(), vj.begin()); // vj = qd
          if (qdd_array.begin() == 0) {
            aj.fill(0);
          }
          else {
            SCITBX_ASSERT(qdd_array[ib].size() == 6);
            ft const* qdd = qdd_array[ib].begin();
            std::copy(qdd, qdd+6, aj.begin()); // aj = qdd
          }
        }
        else {
          matrix_mul(vj, s, qd); // vj = s * qd
          if (qdd_array.begin() == 0) {
            aj.fill(0);
          }
          else {
            matrix_mul(aj, s, qdd_array[ib].const_ref()); // aj = s * qdd
          }
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
                spatial_lib::crm(v[ib]).const_ref(),
                vj.const_ref());
        }
        f[ib] =
            mat_6xn_mul_vec_n(
              body->i_spatial.const_ref(),
              a[ib].const_ref())
          + mat_6xn_mul_vec_n(
              spatial_lib::crf(v[ib]).const_ref(),
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

    //! Returns a packed array of joint force variables (tau_packed).
    af::shared<ft>
    inverse_dynamics_packed(
      af::const_ref<ft> const& qdd_packed=af::const_ref<ft>(0,0),
      af::const_ref<ft> const& f_ext_packed=af::const_ref<ft>(0,0),
      af::const_ref<ft> const& grav_accn=af::const_ref<ft>(0,0))
    {
      af::shared<ft> tau_packed((af::reserve(degrees_of_freedom)));
      af::shared<af::small<ft, 6> >
        tau_array = inverse_dynamics(
          array_packing::unpack_ref_small_6(
            bodies.const_ref(), degrees_of_freedom, qdd_packed).const_ref(),
          array_packing::unpack_ref_tiny<ft, 6>(
            f_ext_packed, bodies.size()).const_ref(),
          grav_accn);
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        tau_packed.extend(tau_array[ib].begin(), tau_array[ib].end());
      }
      SCITBX_ASSERT(tau_packed.size() == degrees_of_freedom);
      return tau_packed;
    }

    /*! Simplified version of Inverse Dynamics via Recursive Newton-Euler
        Algorithm, with all qd, qdd zero, but non-zero external forces.
     */
    af::shared<af::small<ft, 6> >
    f_ext_as_tau(
      af::const_ref<af::tiny<ft, 6> > const& f_ext_array)
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

    /*! Simplified version of Inverse Dynamics via Recursive Newton-Euler
        Algorithm, with all qd, qdd zero, but non-zero external forces.
     */
    af::shared<ft>
    f_ext_as_tau_packed(
      af::const_ref<ft> const& f_ext_packed)
    {
      SCITBX_ASSERT(f_ext_packed.begin() != 0);
      af::shared<ft> tau_packed((af::reserve(degrees_of_freedom)));
      af::shared<af::small<ft, 6> >
        tau_array = f_ext_as_tau(
          array_packing::unpack_ref_tiny<ft, 6>(
            f_ext_packed, bodies.size()).const_ref());
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        tau_packed.extend(tau_array[ib].begin(), tau_array[ib].end());
      }
      SCITBX_ASSERT(tau_packed.size() == degrees_of_freedom);
      return tau_packed;
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
      af::const_ref<ft> const& grav_accn)
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
          c[ib] = mat_6xn_mul_vec_n(
            spatial_lib::crm(v[ib]).const_ref(),
            vj.const_ref());
        }
      }
      boost::scoped_array<vmg> ia(new vmg[nb]);
      for(unsigned ib=0;ib<nb;ib++) {
        ia[ib] = bodies[ib]->i_spatial.deep_copy();
      }
      boost::scoped_array<t6> pa(new t6[nb]);
      for(unsigned ib=0;ib<nb;ib++) {
        pa[ib] = mat_6xn_mul_vec_n(
          spatial_lib::crf(v[ib]).const_ref(),
          mat_6xn_mul_vec_n(
            bodies[ib]->i_spatial.const_ref(),
            v[ib].const_ref()).const_ref());
        if (f_ext_array.begin() != 0) {
          pa[ib] -= f_ext_array[ib];
        }
      }
      v = af::shared<t6>();
      boost::scoped_array<vmg> u(new vmg[nb]);
      boost::scoped_array<vmg> d_inv(new vmg[nb]);
      boost::scoped_array<s6> u_(new s6[nb]);
      for(unsigned ib=nb;ib!=0;) {
        ib--;
        body_t<ft> const* body = bodies[ib].get();
        af::const_ref<ft, af::mat_grid> s = body->joint->motion_subspace();
        vmg d;
        if (s.begin() == 0) {
          u[ib] = ia[ib];
          d = u[ib];
          for(unsigned i=0;i<6;i++) {
            u_[ib].push_back(-pa[ib][i]);
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
        d_inv[ib] = generalized_inverse(d.const_ref());
        if (body->parent != -1) {
          vmg u_d_inv = af::matrix_multiply(
            u[ib].const_ref(),
            d_inv[ib].const_ref());
          vmg ia_ = ia[ib] - af::matrix_multiply_transpose(
            u_d_inv.const_ref(),
            u[ib].const_ref());
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
      pa.reset();
      ia.reset();
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
          SCITBX_ASSERT(qdd_array[ib].size() == 6);
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

    //! Returns a packed array of joint acceleration variables (qdd_packed).
    af::shared<ft>
    forward_dynamics_ab_packed(
      af::const_ref<ft> const& tau_packed=af::const_ref<ft>(0,0),
      af::const_ref<ft> const& f_ext_packed=af::const_ref<ft>(0,0),
      af::const_ref<ft> const& grav_accn=af::const_ref<ft>(0,0))
    {
      af::shared<ft> qdd_packed((af::reserve(degrees_of_freedom)));
      af::shared<af::small<ft, 6> >
        qdd_array = forward_dynamics_ab(
          array_packing::unpack_ref_small_6(
            bodies.const_ref(), degrees_of_freedom, tau_packed).const_ref(),
          array_packing::unpack_ref_tiny<ft, 6>(
            f_ext_packed, bodies.size()).const_ref(),
          grav_accn);
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        qdd_packed.extend(qdd_array[ib].begin(), qdd_array[ib].end());
      }
      SCITBX_ASSERT(qdd_packed.size() == degrees_of_freedom);
      return qdd_packed;
    }
  };

}}} // namespace scitbx::rigid_body::featherstone

#endif // GUARD
