#ifndef CCTBX_REFINEMENT_PARAMETER_MAP_H
#define CCTBX_REFINEMENT_PARAMETER_MAP_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/xray/twin_component.h>

namespace cctbx { namespace xray {

/// Position of the various derivatives or shifts in an array
struct parameter_indices
{
  static int const invariable = -1;

  int site, u_iso, u_aniso, occupancy, fp, fdp;

  parameter_indices()
    : site(invariable), u_iso(invariable), u_aniso(invariable),
      occupancy(invariable), fp(invariable), fdp(invariable)
  {}

  int top() {
    return   fdp       != invariable ? fdp
           : fp        != invariable ? fp
           : occupancy != invariable ? occupancy
           : u_aniso   != invariable ? u_aniso
           : u_iso     != invariable ? u_iso
           : site      != invariable ? site
           :                           invariable;
  }
};

/// An array of parameter_indices
/** It is associated to an array of scatterers: the i-th parameter_indices
    provides the position of the derivatives or shifts of the crystallographic
    parameters of the i-th scatterer.
*/
template <class XRayScattererType>
class parameter_map
{
    af::shared<parameter_indices> map;
    int params;

  public:
    typedef XRayScattererType xray_scatterer_type;
    typedef parameter_indices const *iterator;

    parameter_map() : map() {}

    parameter_map(af::const_ref<xray_scatterer_type> const &scatterers)
      : map(scatterers.size()),
        params(0)
    {
      init_scatterer_part(scatterers);
    }

    int add_independent_scalar() { return ++params; }
private:
    void init_scatterer_part(af::const_ref<xray_scatterer_type> const &scatterers)
    {
      for(int i_sc=0; i_sc < scatterers.size(); ++i_sc) {
        xray_scatterer_type const &sc = scatterers[i_sc];
        parameter_indices &ids = map[i_sc];
        if (sc.flags.grad_site()) {
          ids.site = params;
          params += 3;
        }
        if (sc.flags.use_u_iso() && sc.flags.grad_u_iso()) ids.u_iso = params++;
        if (sc.flags.use_u_aniso() && sc.flags.grad_u_aniso()) {
          ids.u_aniso = params;
          params += 6;
        }
        if (sc.flags.grad_occupancy()) ids.occupancy = params++;
        if (sc.flags.grad_fp()) ids.fp = params++;
        if (sc.flags.grad_fdp()) ids.fdp = params++;
      }
    }

public:
    parameter_indices const &operator[](int i_sc) const { return map[i_sc]; }

    iterator begin() const { return map.begin(); }

    iterator end() const { return map.end(); }

    int size() const { return map.size(); }

    int n_scatterers() const { return size(); }

    int n_parameters() const { return params; }
};

}} // cctbx::xray

#endif // GUARD
