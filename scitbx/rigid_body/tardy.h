#ifndef SCITBX_RIGID_BODY_TARDY_H
#define SCITBX_RIGID_BODY_TARDY_H

#include <scitbx/boost_python/sequence_as_array.h>

#include <scitbx/rigid_body/featherstone.h>
#include <scitbx/rigid_body/body_lib.h>
#include <scitbx/array_family/selections.h>
#include <tbxx/error_utils.hpp>

namespace scitbx { namespace rigid_body {

//! See essence/tardy.py
namespace tardy {

  namespace bp = boost::python;

  //! Helper.
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
            ft abs_cos = fn::absolute(axis * diff);
            if (abs_cos < abs_cos_limit) {
              is_singular = false;
              return;
            }
          }
        }
      }
    }
  };

  //! Helper.
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
        cluster = boost_python::sequence_as<af::shared<unsigned> >(
          clusters[ic]);
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
          fixed_vertices = boost_python::sequence_as<af::shared<unsigned> >(
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
          throw TBXX_UNREACHABLE_ERROR();
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
        SCITBX_ASSERT(body->parent >= 0);
        SCITBX_ASSERT(body->parent < nc);
      }
      result.push_back(body);
    }
    body_lib::set_cb_tree(result.ref());
    return result;
  }

  //! Torsion-angle refinement and dynamics model.
  template <typename FloatType=double>
  struct model : featherstone::system_model<FloatType>
  {
    typedef FloatType ft;

    // constructor arguments
    bp::object labels;
    af::shared<vec3<ft> > sites;
    af::shared<ft> masses;
    bp::object tardy_tree;
    bp::object potential_obj;
    ft near_singular_hinges_angular_tolerance_deg;

    protected:
      boost::optional<af::shared<vec3<ft> > > sites_moved_;
      boost::optional<ft> e_pot_;
      boost::optional<af::shared<vec3<ft> > > d_e_pot_d_sites_;
      boost::optional<af::shared<af::tiny<ft, 6> > > f_ext_array_;
      boost::optional<af::shared<af::small<ft, 6> > > qdd_array_;
    public:

    model() {}

    model(
      bp::object const& labels_,
      af::shared<vec3<ft> > const& sites_,
      af::shared<ft> const& masses_,
      bp::object const& tardy_tree_,
      bp::object const& potential_obj_,
      ft const& near_singular_hinges_angular_tolerance_deg_=5)
    :
      featherstone::system_model<ft>(construct_bodies(
        sites_.const_ref(),
        masses_.const_ref(),
        tardy_tree_.attr("cluster_manager"),
        near_singular_hinges_angular_tolerance_deg_)),
      labels(labels_),
      sites(sites_),
      masses(masses_),
      tardy_tree(tardy_tree_),
      potential_obj(potential_obj_),
      near_singular_hinges_angular_tolerance_deg(
        near_singular_hinges_angular_tolerance_deg_)
    {}

    virtual
    void
    flag_positions_as_changed()
    {
      sites_moved_.reset();
      e_pot_.reset();
      d_e_pot_d_sites_.reset();
      f_ext_array_.reset();
      featherstone::system_model<ft>::flag_positions_as_changed();
    }

    virtual
    void
    flag_velocities_as_changed()
    {
      qdd_array_.reset();
      featherstone::system_model<ft>::flag_velocities_as_changed();
    }

    //! For testing.
    bool
    sites_moved_is_cached() const { return sites_moved_; }

    //! For testing.
    bool
    qdd_array_is_cached() const { return qdd_array_; }

    af::shared<vec3<ft> > const&
    sites_moved()
    {
      if (!sites_moved_) {
        this->aja_array();
        sites_moved_ = af::shared<vec3<ft> >(sites.size());
        unsigned n_done = 0;
        bp::object
          clusters = tardy_tree.attr("cluster_manager").attr("clusters");
        unsigned nb = this->bodies_size();
        for(unsigned ib=0;ib<nb;ib++) {
          rotr3<ft> const& aja = (*this->aja_array_)[ib];
          af::shared<unsigned>
            cluster = boost_python::sequence_as<af::shared<unsigned> >(
              clusters[ib]);
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
        this->jar_array();
        d_e_pot_d_sites();
        unsigned nb = this->bodies_size();
        f_ext_array_ = af::shared<af::tiny<ft, 6> >(af::reserve(nb));
        bp::object
          clusters = tardy_tree.attr("cluster_manager").attr("clusters");
        for(unsigned ib=0;ib<nb;ib++) {
          rotr3<ft> const& cb_0b = this->bodies[ib]->alignment->cb_0b;
          mat3<ft> const& jar = (*this->jar_array_)[ib];
          vec3<ft> f(0,0,0);
          vec3<ft> nc(0,0,0);
          af::shared<unsigned>
            cluster = boost_python::sequence_as<af::shared<unsigned> >(
              clusters[ib]);
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

    /*! \brief Gradients of potential energy (defined via f_ext_array())
        w.r.t. positional coordinates q.
     */
    /*! Not available in Python.
     */
    af::shared<af::small<ft, 7> >
    d_e_pot_d_q()
    {
      unsigned nb = this->bodies_size();
      af::shared<af::small<ft, 7> > result((af::reserve(nb)));
      af::shared<af::small<ft, 6> >
        tau_array = this->f_ext_as_tau(f_ext_array().const_ref());
      for(unsigned ib=0;ib<nb;ib++) {
        result.push_back(
          this->bodies[ib]->joint->tau_as_d_e_pot_d_q(tau_array[ib]));
      }
      return result;
    }

    af::shared<ft>
    d_e_pot_d_q_packed()
    {
      af::shared<ft> result((af::reserve(this->q_packed_size)));
      af::shared<af::small<ft, 7> > unpacked = d_e_pot_d_q();
      SCITBX_ASSERT(unpacked.size() == this->bodies.size());
      unsigned nb = this->bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        result.extend(unpacked[ib].begin(), unpacked[ib].end());
      }
      SCITBX_ASSERT(result.size() == this->q_packed_size);
      return result;
    }

    ft
    e_tot()
    {
      return this->e_kin() + e_pot();
    }

    //! Not available in Python.
    af::shared<af::small<ft, 6> > const&
    qdd_array()
    {
      if (!qdd_array_) {
        qdd_array_ = this->forward_dynamics_ab(
          /*tau_array*/ af::const_ref<af::small<ft, 6> >(0, 0),
          f_ext_array().const_ref(),
          /*grav_accn*/ af::const_ref<ft>(0, 0));
      }
      return *qdd_array_;
    }

    //! For testing.
    af::shared<ft>
    qdd_packed()
    {
      af::shared<ft> result((af::reserve(this->degrees_of_freedom)));
      qdd_array();
      unsigned nb = this->bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        af::small<ft, 6> const& qdd = (*qdd_array_)[ib];
        result.extend(qdd.begin(), qdd.end());
      }
      SCITBX_ASSERT(result.size() == this->degrees_of_freedom);
      return result;
    }

    void
    dynamics_step(
      ft const& delta_t)
    {
      qdd_array();
      unsigned nb = this->bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = this->bodies[ib].get();
        body->joint = body->joint->time_step_position(
          body->qd(), delta_t);
      }
      for(unsigned ib=0;ib<nb;ib++) {
        body_t<ft>* body = this->bodies[ib].get();
        body->set_qd(body->joint->time_step_velocity(
          body->qd(), (*qdd_array_)[ib].const_ref(), delta_t));
      }
      flag_positions_as_changed();
    }
  };

}}} // namespace scitbx::rigid_body::tardy

#endif // GUARD
