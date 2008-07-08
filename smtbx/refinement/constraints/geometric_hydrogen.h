#include <scitbx/constants.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_tiny_io.h>

#include <cctbx/coordinates.h>

namespace smtbx { namespace refinement { namespace constraints {

namespace constants {
  using namespace scitbx::constants;
  static double const tetrahedral_angle = std::acos(-1./3.);
  static double const sin_tetrahedral_angle = std::sin(tetrahedral_angle);
}


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
class stretchable_rotatable_riding_terminal_X_Hn
{
  public:
    typedef XrayScattererType xray_scatterer_type;
    typedef FloatType float_type;
    typedef parameter_map<xray_scatterer_type> parameter_map_type;
    typedef cartesian<float_type> cart_t;
    typedef fractional<float_type> frac_t;

    stretchable_rotatable_riding_terminal_X_Hn(
      int pivot, int pivot_neighbour,
      af::small<int, 3> hydrogens_,
      float_type azimuth_, //degrees
      float_type bond_length,
      bool rotating=true,
      bool stretching=false
      )
      : on_(true),
        i_pivot(pivot), i_pivot_neighbour(pivot_neighbour),
        i_hydrogens(hydrogens_),
        rotating_(rotating),
        stretching_(stretching),
        phi(azimuth_*constants::pi/180),
        l(bond_length)
    {}

    int pivot() { return i_pivot; }

    af::small<int, 3> hydrogens() { return i_hydrogens; }

    bool rotating() { return rotating_; }
    void set_rotating(bool f) { rotating_ = f; }

    bool stretching() { return stretching_; }
    void set_stretching(bool f) { stretching_ = f; }

    boost::tuple<cart_t, cart_t, cart_t> local_cartesian_frame() {
      return boost::make_tuple(e0, e1, e2);
    }

    float_type azimuth() { return phi*180/constants::pi; }
    float_type set_azimuth(float_type phi_) { phi = phi_*constants::pi/180; }

    float_type bond_length() { return l; }
    float_type set_bond_length(float_type l_) { l = l_; }

    void initialise_in_context(
      uctbx::unit_cell const &unit_cell,
      af::const_ref<xray_scatterer_type> const &scatterers,
      af::ref<xray::scatterer_flags> const &constraint_flags,
      std::map<int, xray::scatterer_flags> &already_constrained)
    {
      SMTBX_ASSERT(scatterers.size() == constraint_flags.size())
                  (scatterers.size())
                  (constraint_flags.size());
      cart_t x_pn = unit_cell.orthogonalize(scatterers[i_pivot_neighbour].site);
      cart_t x_p  = unit_cell.orthogonalize(scatterers[i_pivot].site);
      e2 = (x_p - x_pn).normalize();
      e1 = e2.ortho(true);
      e0 = e1.cross(e2);
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

    void compute_gradients(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::const_ref<xray_scatterer_type> const &scatterers,
      parameter_map_type const &crystallographic_parameter_map,
      af::ref<float_type> const &crystallographic_gradients,
      SharedArray1D<float_type> reparametrization_gradients)
    {
      if (!on_) return;

      using namespace constants;
      for(int i=0; i < i_hydrogens.size(); ++i) {
        int i_h = i_hydrogens[i];
        SMTBX_ASSERT(scatterers[i_h].flags.grad_site())(i_h);
        int i_grad_site_pivot = crystallographic_parameter_map[i_pivot].site;
        int i_grad_site_h = crystallographic_parameter_map[i_h].site;
        // riding
        for(int j=0; j < 3; ++j) {
          crystallographic_gradients[i_grad_site_pivot + j]
            += crystallographic_gradients[i_grad_site_h + j];
        }
      }
      if (rotating() || stretching()) {
        i_reparametrization_begin = reparametrization_gradients.size();
        cosines_and_sines phi_(phi, i_hydrogens.size());
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
            cart_t dx_over_dphi
              = l*sin_tetrahedral_angle*(-phi_.sin[i]*e0 + phi_.cos[i]*e1);
            dF_over_dphi += dF_over_dx[i] * dx_over_dphi;
          }
          dF_over_dphi *= pi/180;
          reparametrization_gradients.push_back(dF_over_dphi);
        }

        // stretching
        if (stretching()) {
          float_type dF_over_dl = 0;
          for (int i=0; i < i_hydrogens.size(); ++i) {
            cart_t dx_over_dl
              = sin_tetrahedral_angle*(phi_.cos[i]*e0 + phi_.sin[i]*e1) + e2/3;
            dF_over_dl += dF_over_dx[i] * dx_over_dl;
          }
          reparametrization_gradients.push_back(dF_over_dl);
        }
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
      if (stretching()) {
        float_type delta_l = reparametrization_shifts[i_dF_over_dl()];
        l += delta_l;
      }
      place_constrained_scatterers(unit_cell, site_symmetry_table, scatterers);
    }

    void place_constrained_scatterers(
      uctbx::unit_cell const &unit_cell,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::ref<xray_scatterer_type> const &scatterers)
    {
      using namespace constants;
      cart_t x_pn = unit_cell.orthogonalize(scatterers[i_pivot_neighbour].site);
      cart_t x_p  = unit_cell.orthogonalize(scatterers[i_pivot].site);
      update_local_cartesian_frame(x_p - x_pn);
      cosines_and_sines phi_(phi, i_hydrogens.size());
      for (int i=0; i < i_hydrogens.size(); ++i) {
        int i_h = i_hydrogens[i];
        cart_t x_h = x_p
                     + l*(sin_tetrahedral_angle*(phi_.cos[i]*e0 + phi_.sin[i]*e1)
                          + e2/3);
        scatterers[i_h].site = unit_cell.fractionalize(x_h);
      }
    }

  private:
    bool on_;

    int i_pivot, i_pivot_neighbour;
    af::small<int, 3> i_hydrogens;
    bool rotating_, stretching_;

    cart_t e0, e1, e2;
    float_type phi, l;
    int i_reparametrization_begin;

    //! As the X-Y bond direction changes, we need to update the local frame.
    /** Our method ensure a smooth rotation whereas
        recomputing e0,e1,e2 as in the constructor may result in
        a sudden jump (c.f. implementation of member function "ortho").
    */
    void update_local_cartesian_frame(cart_t e_bond) {
      cart_t f2 = e_bond.normalize();
      e1 = f2.cross(e0);
      e2 = f2;
      e0 = e1.cross(e2);
    }

    struct cosines_and_sines
    {
      af::small<float_type, 3> cos, sin;
      cosines_and_sines(float_type phi, int n) {
        using constants::pi;
        switch (n) {
          case 3:
            cos[2] = std::cos(phi + 4*pi/3);
            sin[2] = std::sin(phi + 4*pi/3);
          case 2:
            cos[1] = std::cos(phi + 2*pi/3);
            sin[1] = std::sin(phi + 2*pi/3);
          case 1:
            cos[0] = std::cos(phi);
            sin[0] = std::sin(phi);
            break;
          default:
            throw std::runtime_error("smtbx: number of geometric hydrogens "
                                     "must be 1,2 or 3");
        }
      }
    };

    int i_dF_over_dphi() { return i_reparametrization_begin; }
    int i_dF_over_dl() { return i_reparametrization_begin + 1; }
};

}}} // namespace smtbx::refinement::constraints
