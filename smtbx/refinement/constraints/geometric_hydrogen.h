#ifndef SMTBX_REFINEMENT_CONSTRAINTS_GEOMETRIC_HYDROGEN_H
#define SMTBX_REFINEMENT_CONSTRAINTS_GEOMETRIC_HYDROGEN_H

#include <scitbx/constants.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_tiny_io.h>

#include <cctbx/coordinates.h>

#include <smtbx/refinement/parameter_map.h>
#include <smtbx/import_scitbx_af.h>
#include <smtbx/import_cctbx.h>

namespace smtbx { namespace refinement { namespace constraints {

namespace constants {
  using namespace scitbx::constants;
  static double const tetrahedral_angle = std::acos(-1./3.);
  static double const sin_tetrahedral_angle = std::sin(tetrahedral_angle);
}


template<typename FloatType, class XrayScattererType,
         class HydrogenIndicesArray,
         template<class> class SharedArray1D=af::shared>
class geometrical_hydrogens
{
  public:
    typedef XrayScattererType xray_scatterer_type;
    typedef FloatType float_type;
    typedef parameter_map<xray_scatterer_type> parameter_map_type;

    geometrical_hydrogens(int pivot,
                          HydrogenIndicesArray hydrogens,
                          float_type bond_length,
                          bool stretching=false
                          )
      : on_(true),
        i_pivot(pivot),
        i_hydrogens(hydrogens),
        stretching_(stretching),
        l(bond_length)
    {}

    int pivot() { return i_pivot; }

    HydrogenIndicesArray hydrogens() { return i_hydrogens; }

    bool stretching() { return stretching_; }
    void set_stretching(bool f) { stretching_ = f; }

    float_type bond_length() { return l; }
    void set_bond_length(float_type l_) { l = l_; }

    void initialise_in_context(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::const_ref<xray_scatterer_type> const &scatterers,
      af::ref<xray::scatterer_flags> const &constraint_flags,
      std::map<int, xray::scatterer_flags> &already_constrained)
    {
      for(int i=0; i < i_hydrogens.size(); ++i) {
        int i_h = i_hydrogens[i];
        xray::scatterer_flags f = constraint_flags[i_h];
        if (!f.grad_site()) {
          already_constrained[i_h] = f;
          on_ = false;
        }
        constraint_flags[i_hydrogens[i]].set_grad_site(false);
      }
    }

    /// Riding
    void ride(uctbx::unit_cell const &unit_cell,
              sgtbx::site_symmetry_table const &site_symmetry_table,
              af::const_ref<xray_scatterer_type> const &scatterers,
              parameter_map_type const &crystallographic_parameter_map,
              af::ref<float_type> const &crystallographic_gradients)
  {
      for(int i=0; i < i_hydrogens.size(); ++i) {
        int i_h = i_hydrogens[i];
        SMTBX_ASSERT(scatterers[i_h].flags.grad_site())(i_h);
        int i_grad_site_pivot = crystallographic_parameter_map[i_pivot].site;
        int i_grad_site_h = crystallographic_parameter_map[i_h].site;
        for(int j=0; j < 3; ++j) {
          crystallographic_gradients[i_grad_site_pivot + j]
            += crystallographic_gradients[i_grad_site_h + j];
        }
      }
    }

  protected:
    bool on_;

    int i_pivot;
    HydrogenIndicesArray i_hydrogens;

    bool stretching_;
    float_type l;

    int i_reparametrization_begin;
};



/// Model of Y-XH3 with tetrahedral angles
/**
  X is referred to as the "pivot" and Y as the "pivot neighbour".

  All angles Hi-X-Hj and Hi-X-Y are tetrahedral.
  All distances X-Hi are equal. That unique distance may be a variable
  parameter if stretching is allowed.
  A free rotation around the bond Y-X is allowed.

  The Hydrogen sites ride on the pivot site.
*/
template<typename FloatType, class XrayScattererType,
         template<class> class SharedArray1D=af::shared>
class terminal_X_Hn
  : public geometrical_hydrogens<FloatType, XrayScattererType,
                                 af::small<int, 3>,
                                 SharedArray1D>
{
  public:
    typedef geometrical_hydrogens<FloatType, XrayScattererType,
                                 af::small<int, 3>,
                                 SharedArray1D>
            base_t;
    using base_t::i_pivot;
    using base_t::i_hydrogens;
    using base_t::l;
    using base_t::on_;
    using base_t::i_reparametrization_begin;

    typedef XrayScattererType xray_scatterer_type;
    typedef FloatType float_type;
    typedef parameter_map<xray_scatterer_type> parameter_map_type;
    typedef cartesian<float_type> cart_t;
    typedef fractional<float_type> frac_t;

    terminal_X_Hn(
      int pivot, int pivot_neighbour,
      af::small<int, 3> hydrogens,
      float_type azimuth_, //degrees
      float_type bond_length,
      bool rotating=true,
      bool stretching=false
      )
      : base_t(pivot, hydrogens, bond_length, stretching),
        i_pivot_neighbour(pivot_neighbour),
        rotating_(rotating),
        phi(azimuth_*constants::pi/180)
    {}

    bool rotating() { return rotating_; }
    void set_rotating(bool f) { rotating_ = f; }

    boost::tuple<cart_t, cart_t, cart_t> local_cartesian_frame() {
      return boost::make_tuple(e0, e1, e2);
    }

    float_type azimuth() { return phi*180/constants::pi; }
    void set_azimuth(float_type phi_) { phi = phi_*constants::pi/180; }

    void initialise_in_context(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::const_ref<xray_scatterer_type> const &scatterers,
      af::ref<xray::scatterer_flags> const &constraint_flags,
      std::map<int, xray::scatterer_flags> &already_constrained)
    {
      base_t::initialise_in_context(unit_cell,
                                    site_symmetry_table,
                                    scatterers,
                                    constraint_flags,
                                    already_constrained);
      if (!on_) return;
      cart_t x_pn = unit_cell.orthogonalize(scatterers[i_pivot_neighbour].site);
      cart_t x_p  = unit_cell.orthogonalize(scatterers[i_pivot].site);
      e2 = (x_p - x_pn).normalize();
      e1 = e2.ortho(true);
      e0 = e1.cross(e2);
    }

    void place_constrained_scatterers(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::ref<xray_scatterer_type> const &scatterers)
    {
      using namespace constants;
      /* As the X-Y bond direction changes, we need to update the local frame.
         Our method ensure a smooth rotation whereas
         recomputing e0,e1,e2 as in initialise_in_context may result in
         a sudden jump (c.f. implementation of member function "ortho").
      */
      cart_t x_pn = unit_cell.orthogonalize(scatterers[i_pivot_neighbour].site);
      cart_t x_p  = unit_cell.orthogonalize(scatterers[i_pivot].site);
      cart_t e_bond = x_p - x_pn;
      cart_t f2 = e_bond.normalize();
      e1 = f2.cross(e0);
      e2 = f2;
      e0 = e1.cross(e2);

      // Compute cosines and sines
      af::small<float_type, 3> cos_phi, sin_phi;
      switch (i_hydrogens.size()) {
        case 3:
          cos_phi[2] = std::cos(phi + 4*pi/3);
          sin_phi[2] = std::sin(phi + 4*pi/3);
        case 2:
          cos_phi[1] = std::cos(phi + 2*pi/3);
          sin_phi[1] = std::sin(phi + 2*pi/3);
        case 1:
          cos_phi[0] = std::cos(phi);
          sin_phi[0] = std::sin(phi);
          break;
        default:
          throw std::runtime_error("smtbx: number of geometric hydrogens "
                                   "must be 1,2 or 3");
      }

      // Place hydrogen's
      for (int i=0; i < i_hydrogens.size(); ++i) {
        int i_h = i_hydrogens[i];
        cart_t x_h = x_p
                     + l*(sin_tetrahedral_angle*(cos_phi[i]*e0 + sin_phi[i]*e1)
                          + e2/3);
        scatterers[i_h].site = unit_cell.fractionalize(x_h);
      }

      // Compute derivatives
      for (int i=0; i < i_hydrogens.size(); ++i) {
        if (rotating()) {
          dx_over_dphi[i]
            = l*sin_tetrahedral_angle*(-sin_phi[i]*e0 + cos_phi[i]*e1);
        }
        if (base_t::stretching()) {
          dx_over_dl[i]
            = sin_tetrahedral_angle*(cos_phi[i]*e0 + sin_phi[i]*e1) + e2/3;
        }
      }
    }

    void compute_gradients(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::const_ref<xray_scatterer_type> const &scatterers,
      parameter_map_type const &crystallographic_parameter_map,
      af::ref<float_type> const &crystallographic_gradients,
      SharedArray1D<float_type> reparametrization_gradients)
    {
      if (!on_) return;

      base_t::ride(unit_cell,
                   site_symmetry_table,
                   scatterers,
                   crystallographic_parameter_map,
                   crystallographic_gradients);

      if (!rotating() && !base_t::stretching()) return;

      using namespace constants;
      i_reparametrization_begin = reparametrization_gradients.size();
      af::tiny<cart_t, 3> dF_over_dx;
      for (int i=0; i < i_hydrogens.size(); ++i) {
        int i_h = i_hydrogens[i];
        int i_grad_site_h = crystallographic_parameter_map[i_h].site;
        frac_t dF_over_dx_frac(&crystallographic_gradients[i_grad_site_h]);
        dF_over_dx[i] = unit_cell.orthogonalize_gradient(dF_over_dx_frac);
      }

      // azimuthal rotation
      if (rotating()) {
        float_type dF_over_dphi = 0;
        for (int i=0; i < i_hydrogens.size(); ++i) {
          dF_over_dphi += dF_over_dx[i] * dx_over_dphi[i];
        }
        dF_over_dphi *= pi/180;
        reparametrization_gradients.push_back(dF_over_dphi);
      }

      // stretching
      if (base_t::stretching()) {
        float_type dF_over_dl = 0;
        for (int i=0; i < i_hydrogens.size(); ++i) {
          dF_over_dl += dF_over_dx[i] * dx_over_dl[i];
        }
        reparametrization_gradients.push_back(dF_over_dl);
      }
    }

    void apply_shifts(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::ref<xray_scatterer_type> const &scatterers,
      parameter_map_type const &crystallographic_parameter_map,
      af::const_ref<float_type> const &crystallographic_shifts,
      af::const_ref<float_type> const &reparametrization_shifts)
    {
      if (!on_) return;

      using namespace constants;
      if (rotating()) {
        float_type delta_phi = reparametrization_shifts[i_dF_over_dphi()];
        phi += delta_phi;
      }
      if (base_t::stretching()) {
        float_type delta_l = reparametrization_shifts[i_dF_over_dl()];
        l += delta_l;
      }
      place_constrained_scatterers(unit_cell, site_symmetry_table, scatterers);
    }

  private:
    int i_pivot_neighbour;
    bool rotating_;

    cart_t e0, e1, e2;
    float_type phi;

    af::small<cart_t, 3> dx_over_dphi, dx_over_dl;

    int i_dF_over_dphi() { return i_reparametrization_begin; }
    int i_dF_over_dl() { return i_reparametrization_begin + 1; }
};


/// Model of X-CH2-Y
/**
  C is referred to as the "pivot" and X and Y as pivot's neighbour 1 and 2.

  All angles Hi-C-X and Hi-C-Y are equal.
  The angle H-C-H depends on XY^2 in a simple linear manner as ShelXL does it.
*/
template<typename FloatType, class XrayScattererType,
         template<class> class SharedArray1D=af::shared>
class secondary_CH2
  : public geometrical_hydrogens<FloatType, XrayScattererType,
                                 af::tiny<int, 2>,
                                 SharedArray1D>
{
  public:
    typedef geometrical_hydrogens<FloatType, XrayScattererType,
                                 af::tiny<int, 2>,
                                 SharedArray1D>
            base_t;
    using base_t::i_pivot;
    using base_t::i_hydrogens;
    using base_t::l;
    using base_t::on_;
    using base_t::i_reparametrization_begin;

    typedef XrayScattererType xray_scatterer_type;
    typedef FloatType float_type;
    typedef parameter_map<xray_scatterer_type> parameter_map_type;
    typedef cartesian<float_type> cart_t;
    typedef fractional<float_type> frac_t;

    static float_type theta0, dtheta_over_dXY_sq;

    secondary_CH2(int pivot,
                  af::tiny<int, 2> pivot_neighbours,
                  af::tiny<int, 2> hydrogens,
                  float_type bond_length,
                  bool stretching=false)
      : base_t(pivot, hydrogens, bond_length, stretching),
        i_pivot_neighbours(pivot_neighbours)
    {}

    void place_constrained_scatterers(
      uctbx::unit_cell const &uc,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::ref<xray_scatterer_type> const &scatterers)
    {
      /* Recompute local frame: (C,e0,e1) is the bisecting plane
         of the angle X-C-Y with e0 bisecting X-C-Y
      */
      cart_t x_p  = uc.orthogonalize(scatterers[i_pivot].site);
      cart_t x_pn_1 = uc.orthogonalize(scatterers[i_pivot_neighbours[0]].site),
             x_pn_2 = uc.orthogonalize(scatterers[i_pivot_neighbours[1]].site);
      cart_t u_pn_1 = (x_p - x_pn_1).normalize(),
             u_pn_2 = (x_p - x_pn_2).normalize();
      e0 = (u_pn_1 + u_pn_2).normalize();
      cart_t e2 = (u_pn_2 - u_pn_1).normalize();
      e1 = e2.cross(e0);

      // Compute cosine and sine
      float_type d_XY_sq = (x_pn_2 - x_pn_1).length_sq();
      float_type theta = theta0 - dtheta_over_dXY_sq*d_XY_sq;
      float_type c = std::cos(theta), s = std::sin(theta);

      // Place hydrogen's
      cart_t site_h_1 = x_p + l*(c*e0 + s*e1);
      cart_t site_h_2 = x_p + l*(c*e0 - s*e1);
      scatterers[i_hydrogens[0]].site = uc.fractionalize(site_h_1);
      scatterers[i_hydrogens[1]].site = uc.fractionalize(site_h_2);

      // Compute derivatives
      dx_over_dl[0] = c*e0 + s*e1;
      dx_over_dl[1] = c*e0 - s*e1;
    }

    void compute_gradients(
      uctbx::unit_cell const &uc,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::const_ref<xray_scatterer_type> const &scatterers,
      parameter_map_type const &crystallographic_parameter_map,
      af::ref<float_type> const &crystallographic_gradients,
      SharedArray1D<float_type> reparametrization_gradients)
    {
      if (!on_) return;

      base_t::ride(uc,
                   site_symmetry_table,
                   scatterers,
                   crystallographic_parameter_map,
                   crystallographic_gradients);

      if (!base_t::stretching()) return;

      i_reparametrization_begin = reparametrization_gradients.size();
      float_type dF_over_dl = 0;
      for (int i=0; i < 2; ++i) {
        int i_h = i_hydrogens[i];
        int i_grad_site_h = crystallographic_parameter_map[i_h].site;
        frac_t dF_over_dx_frac(&crystallographic_gradients[i_grad_site_h]);
        cart_t dF_over_dx = uc.orthogonalize_gradient(dF_over_dx_frac);
        dF_over_dl += dF_over_dx * dx_over_dl[i];
      }
      reparametrization_gradients.push_back(dF_over_dl);
    }

    void apply_shifts(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::ref<xray_scatterer_type> const &scatterers,
      parameter_map_type const &crystallographic_parameter_map,
      af::const_ref<float_type> const &crystallographic_shifts,
      af::const_ref<float_type> const &reparametrization_shifts)
    {
      if (!on_) return;
      place_constrained_scatterers(unit_cell, site_symmetry_table, scatterers);
    }

  private:
    cart_t e0, e1;
    af::tiny<cart_t, 2> dx_over_dl;
    af::tiny<int, 2> i_pivot_neighbours;
};


// Numbers from ShelXL (file xl.f, line 8411)
template<typename FloatType, class XrayScattererType,
         template<class> class SharedArray1D>
FloatType
secondary_CH2<FloatType, XrayScattererType, SharedArray1D>::
theta0 = 1.0376;

template<typename FloatType, class XrayScattererType,
         template<class> class SharedArray1D>
FloatType
secondary_CH2<FloatType, XrayScattererType, SharedArray1D>::
dtheta_over_dXY_sq = -0.0349;


/// Model of tertiary CH
/** All angles Hi-C-X are equal.
*/
template<typename FloatType, class XrayScattererType,
         template<class> class SharedArray1D=af::shared>
class tertiary_CH
  : public geometrical_hydrogens<FloatType, XrayScattererType,
                                 af::tiny<int, 1>,
                                 SharedArray1D>
{
  public:
    typedef geometrical_hydrogens<FloatType, XrayScattererType,
                                 af::tiny<int, 1>,
                                 SharedArray1D>
            base_t;
    using base_t::i_pivot;
    using base_t::i_hydrogens;
    using base_t::l;
    using base_t::on_;
    using base_t::i_reparametrization_begin;

    typedef XrayScattererType xray_scatterer_type;
    typedef FloatType float_type;
    typedef parameter_map<xray_scatterer_type> parameter_map_type;
    typedef cartesian<float_type> cart_t;
    typedef fractional<float_type> frac_t;


    tertiary_CH(int pivot,
                af::tiny<int, 3> pivot_neighbours,
                int hydrogen,
                float_type bond_length,
                bool stretching=false)
      : base_t(pivot, af::tiny<int,1>(hydrogen), bond_length, stretching),
        i_pivot_neighbours(pivot_neighbours)
    {}

    void place_constrained_scatterers(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::ref<xray_scatterer_type> const &scatterers)
    {

      cart_t x_p = unit_cell.orthogonalize(scatterers[i_pivot].site);
      e0 = cart_t(0,0,0);
      cart_t
        x_X = unit_cell.orthogonalize(scatterers[i_pivot_neighbours[0]].site),
        x_Y = unit_cell.orthogonalize(scatterers[i_pivot_neighbours[1]].site),
        x_Z = unit_cell.orthogonalize(scatterers[i_pivot_neighbours[2]].site);
      cart_t u_XC = (x_p - x_X).normalize(),
             u_YC = (x_p - x_Y).normalize(),
             u_ZC = (x_p - x_Z).normalize();
      cart_t u = u_XC - u_YC;
      cart_t v = u_YC - u_ZC;
      e0 = u.cross(v).normalize();
      if (e0*(u_XC + u_YC + u_ZC) < 0) e0 = -e0;
      cart_t x_h = x_p + l*e0;
      scatterers[i_hydrogens[0]].site = unit_cell.fractionalize(x_h);
    }

    void compute_gradients(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::const_ref<xray_scatterer_type> const &scatterers,
      parameter_map_type const &crystallographic_parameter_map,
      af::ref<float_type> const &crystallographic_gradients,
      SharedArray1D<float_type> reparametrization_gradients)
    {
      if (!on_) return;

      base_t::ride(unit_cell,
                   site_symmetry_table,
                   scatterers,
                   crystallographic_parameter_map,
                   crystallographic_gradients);

      if (base_t::stretching()) {
        i_reparametrization_begin = reparametrization_gradients.size();
        int i_grad_site_h = crystallographic_parameter_map[i_hydrogens[0]].site;
        frac_t dF_over_dx_frac(&crystallographic_gradients[i_grad_site_h]);
        cart_t dF_over_dx = unit_cell.orthogonalize_gradient(dF_over_dx_frac);
        float_type dF_over_dl = dF_over_dx * e0;
        reparametrization_gradients.push_back(dF_over_dl);
      }
    }

    void apply_shifts(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::ref<xray_scatterer_type> const &scatterers,
      parameter_map_type const &crystallographic_parameter_map,
      af::const_ref<float_type> const &crystallographic_shifts,
      af::const_ref<float_type> const &reparametrization_shifts)
    {
      if (!on_) return;
      place_constrained_scatterers(unit_cell, site_symmetry_table, scatterers);
    }

  private:
    cart_t e0;
    af::tiny<int, 3> i_pivot_neighbours;

};

}}} // namespace smtbx::refinement::constraints

#endif // GUARD
