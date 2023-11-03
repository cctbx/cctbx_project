
#include <mmtbx/error.h>
#include <iotbx/pdb/hierarchy.h>
#include <scitbx/vec3.h>
#include <boost/optional.hpp>
#include <limits>
#include <cstring>

#define HB_Q1 0.42
#define HB_Q2 0.2
#define HB_F 332
#define HB_F_Q1_Q2 27.888

namespace mmtbx { namespace secondary_structure { namespace dssp {
  using namespace iotbx::pdb::hierarchy;

  inline double hbond_energy (
      scitbx::vec3<double> c_xyz,
      scitbx::vec3<double> h_xyz,
      scitbx::vec3<double> o_xyz,
      scitbx::vec3<double> n_xyz) {
    double r_ON = (o_xyz - n_xyz).length();
    double r_CH = (c_xyz - h_xyz).length();
    double r_OH = (o_xyz - h_xyz).length();
    double r_CN = (c_xyz - n_xyz).length();
    double E = HB_F_Q1_Q2 * (1./r_ON + 1./r_CH - 1./r_OH - 1./r_CN);
    return E;
  }

  // given a backbone nitrogen atom, calculate the position of the attached
  // hydrogen atom.
  boost::optional< scitbx::vec3<double> >
  get_n_h_position (
    atom const& N,
    double nh_bond_length=1.01)
  {
    bool have_caN = false;
    bool have_cN = false;
    scitbx::vec3<double> N_xyz = N.data->xyz;
    scitbx::vec3<double> caN_xyz(0,0,0);
    scitbx::vec3<double> cN_xyz(0,0,0);
    // find CA attached to N
    boost::optional<atom_group> ag_N = N.parent();
    MMTBX_ASSERT(ag_N);
    unsigned n_ats = ag_N->atoms_size();
    std::vector<atom> agN_atoms = ag_N->atoms();
    for(unsigned i_at=0;i_at<n_ats;i_at++) {
      atom const& a = agN_atoms[i_at];
      if (std::strcmp(a.data->name.elems, " CA ") == 0) {
        caN_xyz = scitbx::vec3<double>(a.data->xyz);
        have_caN = true;
        break;
      }
    }
    if (! have_caN) {
      return boost::optional< scitbx::vec3<double> >();
    }
    // find previous C attached to N
    boost::optional<residue_group> rg_N = ag_N->parent();
    boost::optional<chain> chn_N = rg_N->parent();
    unsigned n_rgs = chn_N->residue_groups_size();
    std::vector<residue_group> rgs = chn_N->residue_groups();
    residue_group prev_rg;
    bool have_prev_rg = false;
    for (unsigned i_rg = 0; i_rg < n_rgs; i_rg++) {
      if (rgs[i_rg].memory_id() == rg_N->memory_id()) {
        break;
      } else {
        prev_rg = rgs[i_rg];
        have_prev_rg = true;
      }
    }
    if (! have_prev_rg) {
      return boost::optional< scitbx::vec3<double> >();
    }
    unsigned n_ags = prev_rg.atom_groups_size();
    std::vector<atom_group> const& ags = prev_rg.atom_groups();
    for (unsigned i_ag = 0; i_ag < n_ags; i_ag++) {
      if (ags[i_ag].data->altloc[0] == ag_N->data->altloc[0]) {
        n_ats = ags[i_ag].atoms_size();
        std::vector<atom> prev_atoms = ags[i_ag].atoms();
        for (unsigned i_at = 0; i_at < n_ats; i_at++) {
          atom const& a = prev_atoms[i_at];
          if (std::strcmp(a.data->name.elems, " C  ") == 0) {
            cN_xyz = scitbx::vec3<double>(a.data->xyz);
            have_cN = true;
            break;
          }
        }
        break;
      }
    }
    if (! have_cN) {
      return boost::optional< scitbx::vec3<double> >();
    }
    scitbx::vec3<double> midpoint = (caN_xyz + cN_xyz) / 2;
    scitbx::vec3<double> vec_nm = N_xyz - midpoint;
    scitbx::vec3<double>hN_xyz = N_xyz + (vec_nm.normalize() * nh_bond_length);
    return boost::optional< scitbx::vec3<double> >(hN_xyz);
  }


  // given carbonyl oxygen and amino nitrogen atoms, calculate the hydrogen
  // bond energy between them.
  boost::optional<double>
  get_o_n_hbond_energy (
    atom const& O,
    atom const& N,
    double nh_bond_length=1.01)
  {
    scitbx::vec3<double> N_xyz = N.data->xyz;
    scitbx::vec3<double> O_xyz = O.data->xyz;
    bool have_cO = false;
    scitbx::vec3<double> cO_xyz(0,0,0);
    boost::optional< scitbx::vec3<double> > hN_xyz;
    // find C attached to O
    boost::optional<atom_group> ag_O = O.parent();
    MMTBX_ASSERT(ag_O);
    unsigned n_ats = ag_O->atoms_size();
    std::vector<atom> const& atoms = ag_O->atoms();
    for(unsigned i_at=0;i_at<n_ats;i_at++) {
      atom const& a = atoms[i_at];
      if (std::strcmp(a.data->name.elems, " C  ") == 0) {
        cO_xyz = scitbx::vec3<double>(a.data->xyz);
        have_cO = true;
        break;
      }
    }
    if (! have_cO) {
      return boost::optional<double>();
    }
    double rCN = (cO_xyz - N_xyz).length_sq();
    if (rCN > 49.0) {
      return boost::optional<double>();
    }
    // now calculate the hydrogen position
    hN_xyz = get_n_h_position(N, nh_bond_length);
    if (! hN_xyz) {
      return boost::optional<double>();
    }
    // finally we can calculate the energy
    return hbond_energy(cO_xyz, *hN_xyz, O_xyz, N_xyz);
  }

}}} // namespace mmtbx::secondary_structure::dssp
