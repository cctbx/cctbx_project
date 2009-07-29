#ifndef SCITBX_RIGID_BODY_TARDY_H
#define SCITBX_RIGID_BODY_TARDY_H

#include <boost/python/import.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/object.hpp>

#include <scitbx/rigid_body/body_lib.h>
#include <scitbx/rigid_body/featherstone.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/selections.h>
#include <scitbx/optional_copy.h>

namespace scitbx { namespace rigid_body { namespace tardy {

  namespace bp = boost::python;

  template <typename ElementType>
  af::shared<ElementType>
  python_sequence_as_af_shared(
    bp::object const& seq)
  {
    bp::ssize_t n = bp::len(seq);
    af::shared<ElementType>
      result((af::reserve(boost::numeric_cast<std::size_t>(n))));
    for(bp::ssize_t i=0;i<n;i++) {
      result.push_back(bp::extract<ElementType>(seq[i])());
    }
    return result;
  }

  template <typename FloatType>
  struct is_singular_revolute
  {
    typedef FloatType ft;

    vec3<ft> const& pivot;
    vec3<ft> axis;
    bool is_singular;

    is_singular_revolute(
      vec3<ft> const& normal_sites_0,
      vec3<ft> const& normal_sites_1,
      af::const_ref<vec3<ft> > const& body_sites,
      ft const& abs_cos_limit)
    :
      pivot(normal_sites_1),
      axis(pivot - normal_sites_0),
      is_singular(true)
    {
      ft axis_length = axis.length();
      if (axis_length != 0) {
        axis /= axis_length;
        for(std::size_t i=0;i<body_sites.size();i++) {
          vec3<ft> diff = body_sites[i] - pivot;
          ft diff_length = diff.length();
          if (diff_length != 0) {
            diff /= diff_length;
            ft abs_cos = axis * diff;
            if (abs_cos < abs_cos_limit) {
              is_singular = false;
              return;
            }
          }
        }
      }
    }
  };

  template <typename FloatType>
  af::shared<shared_ptr<body_t<FloatType> > >
  construct_bodies(
    af::const_ref<vec3<FloatType> > const& sites,
    af::const_ref<FloatType> const& masses,
    bp::object const& cluster_manager,
    FloatType const& near_singular_hinges_angular_tolerance_deg=5)
  {
    SCITBX_ASSERT(masses.size() == sites.size());
    bp::object none;
    typedef FloatType ft;
    af::shared<shared_ptr<body_t<ft> > > result;
    ft abs_cos_limit = fn::absolute(std::cos(
      near_singular_hinges_angular_tolerance_deg * constants::pi_180));
    bp::object fvgci = cluster_manager.attr(
      "fixed_vertices_given_cluster_index_dict")();
    bp::object clusters = cluster_manager.attr("clusters");
    unsigned nc = boost::numeric_cast<unsigned>(bp::len(clusters));
    for(unsigned ic=0;ic<nc;ic++) {
      af::shared<unsigned>
        cluster = python_sequence_as_af_shared<unsigned>(clusters[ic]);
      af::shared<vec3<ft> >
        body_sites = af::select(sites, cluster.const_ref());
      af::shared<ft>
        body_masses = af::select(masses, cluster.const_ref());
      bp::object he = cluster_manager.attr("hinge_edges")[ic];
      SCITBX_ASSERT(bp::len(he) == 2);
      int he_0 = bp::extract<int>(he[0])();
      int he_1 = bp::extract<int>(he[1])();
      SCITBX_ASSERT(he_0 >= -1);
      SCITBX_ASSERT(he_1 >= 0);
      bp::object fixed_vertices_ = fvgci.attr("get")(ic);
      shared_ptr<body_t<ft> > body;
      if (fixed_vertices_.ptr() != none.ptr()) {
        af::shared<unsigned>
          fixed_vertices = python_sequence_as_af_shared<unsigned>(
            fixed_vertices_);
        if (   fixed_vertices.size() > 2
            || fixed_vertices.size() == cluster.size()) {
          body = shared_ptr<body_t<ft> >(new
            body_lib::zero_dof<ft>(
              body_sites.const_ref(),
              body_masses.const_ref()));
        }
        else if (fixed_vertices.size() == 1) {
          body = shared_ptr<body_t<ft> >(new
            body_lib::spherical<ft>(
              body_sites.const_ref(),
              body_masses.const_ref(),
              /*pivot*/ sites[fixed_vertices[0]]));
        }
        else if (fixed_vertices.size() == 2) {
          is_singular_revolute<ft> decision(
            /*normal_sites_0*/ sites[fixed_vertices[0]],
            /*normal_sites_1*/ sites[fixed_vertices[1]],
            body_sites.const_ref(),
            abs_cos_limit);
          if (decision.is_singular) {
            body = shared_ptr<body_t<ft> >(new
              body_lib::zero_dof<ft>(
                body_sites.const_ref(),
                body_masses.const_ref()));
          }
          else {
            body = shared_ptr<body_t<ft> >(new
              body_lib::revolute<ft>(
                body_sites.const_ref(),
                body_masses.const_ref(),
                decision.pivot,
                decision.axis));
          }
        }
        else {
          throw SCITBX_INTERNAL_ERROR(); // unreachable
        }
        body->parent = -1;
      }
      else if (he_0 == -1) {
        if (body_sites.size() == 1) {
          body = shared_ptr<body_t<ft> >(new
            body_lib::translational<ft>(
              body_sites.const_ref(),
              body_masses.const_ref()));
        }
        else {
          body = shared_ptr<body_t<ft> >(new
            body_lib::six_dof<ft>(
              body_sites.const_ref(),
              body_masses.const_ref()));
        }
        body->parent = -1;
      }
      else {
        vec3<ft> const& normal_sites_0 = sites[he_0];
        vec3<ft> const& normal_sites_1 = sites[he_1];
        body = shared_ptr<body_t<ft> >(new
          body_lib::revolute<ft>(
            body_sites.const_ref(),
            body_masses.const_ref(),
            /*pivot*/ normal_sites_1,
            /*normal*/ (normal_sites_1-normal_sites_0).normalize()));
        body->parent = bp::extract<int>(
          cluster_manager.attr("cluster_indices")[he_1])();
      }
      result.push_back(body);
    }
    return result;
  }

  template <typename FloatType=double>
  struct model : boost::noncopyable
  {
    typedef FloatType ft;

    // constructor arguments
    bp::object labels;
    af::shared<vec3<ft> > sites;
    af::shared<ft> masses;
    bp::object tardy_tree;
    bp::object potential_obj;
    ft near_singular_hinges_angular_tolerance_deg;

    // set in constructor
    af::shared<shared_ptr<body_t<FloatType> > > bodies;
    unsigned degrees_of_freedom;

    // dynamically maintained
    protected:
      boost::optional<featherstone::system_model<ft> >
        featherstone_system_model_;
      boost::optional<af::shared<rotr3<ft> > > aja_array_;
      boost::optional<af::shared<mat3<ft> > > jar_array_;
      boost::optional<af::shared<vec3<ft> > > sites_moved_;
      boost::optional<ft> e_pot_;
      boost::optional<af::shared<vec3<ft> > > d_e_pot_d_sites_;
      boost::optional<af::shared<af::tiny<ft, 6> > > f_ext_array_;
      boost::optional<af::shared<af::small<ft, 6> > > qdd_array_;
      boost::optional<ft> e_kin_;
    public:

    unsigned
    bodies_size() const
    {
      return boost::numeric_cast<unsigned>(bodies.size());
    }

    model() {}

    model(
      bp::object const& labels_,
      af::shared<vec3<ft> > const& sites_,
      af::shared<ft> const& masses_,
      bp::object const& tardy_tree_,
      bp::object const& potential_obj_,
      ft const& near_singular_hinges_angular_tolerance_deg_=5)
    :
      labels(labels_),
      sites(sites_),
      masses(masses_),
      tardy_tree(tardy_tree_),
      potential_obj(potential_obj_),
      near_singular_hinges_angular_tolerance_deg(
        near_singular_hinges_angular_tolerance_deg_),
      bodies(construct_bodies(
        sites.const_ref(),
        masses.const_ref(),
        tardy_tree.attr("cluster_manager"),
        near_singular_hinges_angular_tolerance_deg)),
      degrees_of_freedom(0)
    {
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        degrees_of_freedom += body->joint->degrees_of_freedom;
      }
      flag_positions_as_changed();
    }

    void
    flag_positions_as_changed()
    {
      featherstone_system_model_.reset();
      aja_array_.reset();
      jar_array_.reset();
      sites_moved_.reset();
      e_pot_.reset();
      d_e_pot_d_sites_.reset();
      f_ext_array_.reset();
      flag_velocities_as_changed();
    }

    void
    flag_velocities_as_changed()
    {
      qdd_array_.reset();
      e_kin_.reset();
    }

    af::shared<std::size_t>
    root_indices() const
    {
      af::shared<std::size_t> result;
      std::size_t nb = bodies.size();
      for(std::size_t ib=0;ib<nb;ib++) {
        body_t<ft> const* body = bodies[ib].get();
        if (body->parent == -1) {
          result.push_back(ib);
        }
      }
      return result;
    }

    af::shared<af::tiny<std::size_t, 2> >
    number_of_sites_in_each_tree() const
    {
      af::shared<af::tiny<std::size_t, 2> > result;
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
      return result;
    }

    af::shared<std::pair<int, double> >
    sum_of_masses_in_each_tree() const
    {
      af::shared<std::pair<int, double> > result;
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
      optional_copy<af::shared<af::tiny<std::size_t, 2> > > nosiet; \
      if (number_of_sites_in_each_tree.begin() == 0) { \
        nosiet = this->number_of_sites_in_each_tree(); \
        number_of_sites_in_each_tree = nosiet->const_ref(); \
      } \
      SCITBX_ASSERT(number_of_sites_in_each_tree.size() == bodies.size()); \
      std::size_t nb = bodies.size(); \
      for( \
        af::tiny<std::size_t, 2> const* nosiet_it=nosiet->begin(); \
        nosiet_it!=nosiet->end(); \
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
    featherstone::system_model<ft>&
    featherstone_system_model()
    {
      if (!featherstone_system_model_) {
        featherstone_system_model_ = featherstone::system_model<ft>(bodies);
      }
      return *featherstone_system_model_;
    }

    //! Not available in Python.
    af::shared<rotr3<ft> > const&
    aja_array()
    {
      if (!aja_array_) {
        unsigned nb = bodies_size();
        aja_array_ = af::shared<rotr3<ft> >();
        aja_array_->reserve(nb);
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
        jar_array_ = af::shared<mat3<ft> >();
        jar_array_->reserve(nb);
        for(unsigned ib=0;ib<nb;ib++) {
          body_t<ft> const* body = bodies[ib].get();
          mat3<ft> jar = body->joint->cb_ps.r * body->alignment->cb_0b.r;
          if (body->parent != -1) {
            jar = jar * (*aja_array_)[body->parent].r.transpose();
          }
          jar_array_->push_back(jar);
        }
      }
      return *jar_array_;
    }

    af::shared<vec3<ft> > const&
    sites_moved()
    {
      if (!sites_moved_) {
        aja_array();
        sites_moved_ = af::shared<vec3<ft> >(sites.size());
        unsigned n_done = 0;
        bp::object
          clusters = tardy_tree.attr("cluster_manager").attr("clusters");
        unsigned nb = bodies_size();
        for(unsigned ib=0;ib<nb;ib++) {
          rotr3<ft> const& aja = (*aja_array_)[ib];
          af::shared<unsigned>
            cluster = python_sequence_as_af_shared<unsigned>(clusters[ib]);
          unsigned n = boost::numeric_cast<unsigned>(cluster.size());
          for(unsigned i=0;i<n;i++) {
            unsigned i_seq = cluster[i];
            (*sites_moved_)[i_seq] = aja * sites[i_seq];
            n_done++;
          }
        }
        SCITBX_ASSERT(n_done == sites.size());
      }
      return *sites_moved_;
    }

    ft const&
    e_pot()
    {
      if (!e_pot_) {
        bp::object none;
        if (potential_obj.ptr() == none.ptr()) {
          e_pot_ = 0;
        }
        else {
          e_pot_ = bp::extract<ft>(potential_obj.attr("e_pot")(
            sites_moved()))();
        }
      }
      return *e_pot_;
    }

    af::shared<vec3<ft> > const&
    d_e_pot_d_sites()
    {
      if (!d_e_pot_d_sites_) {
        bp::object none;
        if (potential_obj.ptr() == none.ptr()) {
          d_e_pot_d_sites_ = af::shared<vec3<ft> >(
            sites.size(), vec3<ft>(0,0,0));
        }
        else {
          d_e_pot_d_sites_ = bp::extract<af::shared<vec3<ft> > >(
            potential_obj.attr("d_e_pot_d_sites")(
              sites_moved()))();
        }
      }
      return *d_e_pot_d_sites_;
    }

    //! Not available in Python.
    af::shared<af::tiny<ft, 6> > const&
    f_ext_array()
    {
      if (!f_ext_array_) {
        jar_array();
        d_e_pot_d_sites();
        f_ext_array_ = af::shared<af::tiny<ft, 6> >();
        bp::object
          clusters = tardy_tree.attr("cluster_manager").attr("clusters");
        unsigned nb = bodies_size();
        for(unsigned ib=0;ib<nb;ib++) {
          rotr3<ft> const& cb_0b = bodies[ib]->alignment->cb_0b;
          mat3<ft> const& jar = (*jar_array_)[ib];
          vec3<ft> f(0,0,0);
          vec3<ft> nc(0,0,0);
          af::shared<unsigned>
            cluster = python_sequence_as_af_shared<unsigned>(clusters[ib]);
          unsigned n = boost::numeric_cast<unsigned>(cluster.size());
          for(unsigned i=0;i<n;i++) {
            unsigned i_seq = cluster[i];
            vec3<ft> const& s = sites[i_seq];
            vec3<ft> force_bf = -(jar * (*d_e_pot_d_sites_)[i_seq]);
            f += force_bf;
            nc += (cb_0b * s).cross(force_bf);
          }
          f_ext_array_->push_back(spatial_lib::as_tiny_6(nc, f));
        }
      }
      return *f_ext_array_;
    }

    //! Not available in Python.
    af::shared<af::small<ft, 6> > const&
    qdd_array()
    {
      if (!qdd_array_) {
        qdd_array_ = featherstone_system_model().forward_dynamics_ab(
          /*tau_array*/ af::const_ref<af::small<ft, 6> >(0, 0),
          f_ext_array().const_ref(),
          /*grav_accn*/ af::const_ref<ft>(0, 0));
      }
      return *qdd_array_;
    }

    ft const&
    e_kin()
    {
      if (!e_kin_) {
        e_kin_ = featherstone_system_model().e_kin();
      }
      return *e_kin_;
    }

    ft
    e_tot()
    {
      return e_kin() + e_pot();
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
        std::copy(joint_qd_zero.begin(), joint_qd_zero.end(), body_qd.begin());
      }
      flag_velocities_as_changed();
    }

    boost::optional<af::shared<ft> >
    assign_random_velocities(
      boost::optional<ft> const& e_kin_target=boost::optional<ft>(),
      ft const& e_kin_epsilon=1e-12,
      bp::object random_gauss=bp::object())
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
      af::shared<ft> qd_e_kin_scales =
        featherstone_system_model().qd_e_kin_scales(e_kin_epsilon);
      if (degrees_of_freedom != 0) {
        qd_e_kin_scales *= boost::numeric_cast<ft>(
          std::sqrt(
              work_e_kin_target
            / boost::numeric_cast<ft>(degrees_of_freedom)));
      }
      bp::object none;
      if (random_gauss.ptr() == none.ptr()) {
        random_gauss = bp::import("random").attr("gauss");
      }
      unsigned i_qd = 0;
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        af::small<ft, 6> qd_new(af::adapt(body->joint->qd_zero()));
        unsigned n = boost::numeric_cast<unsigned>(qd_new.size());
        for(unsigned i=0;i<n;i++) {
          qd_new[i] += bp::extract<ft>(
            random_gauss(/*mu*/ 0, /*sigma*/ qd_e_kin_scales[i_qd]))();
        }
        i_qd += n;
        body->set_qd(qd_new);
      }
      SCITBX_ASSERT(i_qd == degrees_of_freedom);
      flag_velocities_as_changed();
      if (e_kin_target) {
        reset_e_kin(*e_kin_target, e_kin_epsilon);
      }
      return boost::optional<af::shared<ft> >(qd_e_kin_scales);
    }

    void
    dynamics_step(
      ft const& delta_t)
    {
      qdd_array();
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        body->joint = body->joint->time_step_position(
          body->qd(), delta_t);
      }
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        body->set_qd(body->joint->time_step_velocity(
          body->qd(), (*qdd_array_)[ib].const_ref(), delta_t));
      }
      flag_positions_as_changed();
    }

    //! Not available in Python.
    af::shared<af::small<ft, 7> >
    d_pot_d_q()
    {
      return featherstone_system_model().d_pot_d_q(f_ext_array().const_ref());
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
      return result;
    }

    void
    unpack_q(
      af::const_ref<ft> const& packed_q)
    {
      unsigned i = 0;
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        unsigned n = body->joint->q_size;
        body->joint = body->joint->new_q(af::const_ref<ft>(&packed_q[i], n));
        i += n;
      }
      SCITBX_ASSERT(i == packed_q.size());
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
      return result;
    }

    void
    unpack_qd(
      af::const_ref<ft> const& packed_qd)
    {
      unsigned i = 0;
      unsigned nb = bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = bodies[ib].get();
        unsigned n = body->joint->degrees_of_freedom;
        body->set_qd(af::small<ft, 6>(af::adapt(
          af::const_ref<ft>(&packed_qd[i], n))));
        i += n;
      }
      SCITBX_ASSERT(i == packed_qd.size());
      flag_velocities_as_changed();
    }
  };

}}} // namespace scitbx::rigid_body::tardy

#endif // GUARD
