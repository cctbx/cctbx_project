#ifndef MMTBX_GEOMETRY_RESTRAINTS_ROTAMER_H
#define MMTBX_GEOMETRY_RESTRAINTS_ROTAMER_H

#include <scitbx/array_family/shared.h>
#include <mmtbx/error.h>
#include <cctbx/geometry_restraints/dihedral.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>

#include <string>

#define PI_OVER_180 0.017453292519943295

namespace mmtbx { namespace geometry_restraints {
  using cctbx::geometry_restraints::dihedral;
  using cctbx::geometry_restraints::dihedral_proxy;
  namespace af = scitbx::af;

  struct rotamer_proxy
  {
    std::string residue_name;
    std::string residue_type;
    size_t n_angles;
    bool have_phi_psi;
    af::tiny<unsigned, 16> chi_i_seqs;
    af::tiny<unsigned, 5> phi_psi_i_seqs;
    size_t residue_index;

    // default initializer
    rotamer_proxy () {}

    // both phi-psi and chi angles defined
    rotamer_proxy (
      std::string const& residue_name_,
      std::string const& residue_type_,
      size_t n_angles_,
      af::tiny<unsigned, 5> const& phi_psi_i_seqs_,
      af::tiny<unsigned, 4> const& chi1_i_seqs,
      af::tiny<unsigned, 4> const& chi2_i_seqs,
      af::tiny<unsigned, 4> const& chi3_i_seqs,
      af::tiny<unsigned, 4> const& chi4_i_seqs,
      size_t residue_index_=1)
    :
      residue_name(residue_name_),
      residue_type(residue_type_),
      n_angles(n_angles_),
      have_phi_psi(true),
      phi_psi_i_seqs(phi_psi_i_seqs_),
      residue_index(residue_index_)
    {
      MMTBX_ASSERT(residue_index > 0);
      initialize_i_seqs(chi1_i_seqs, chi2_i_seqs, chi3_i_seqs, chi4_i_seqs);
    }

    // phi-psi only
    rotamer_proxy(
      std::string const& residue_name_,
      std::string const& residue_type_,
      af::tiny<unsigned, 5> const& phi_psi_i_seqs_,
      size_t residue_index_=1)
    :
      residue_name(residue_name_),
      residue_type(residue_type_),
      n_angles(0),
      have_phi_psi(true),
      phi_psi_i_seqs(phi_psi_i_seqs_),
      residue_index(residue_index_)
    {
      MMTBX_ASSERT(residue_index > 0);
    }

    // chi angles only
    rotamer_proxy(
      std::string const& residue_name_,
      size_t n_angles_,
      af::tiny<unsigned, 4> const& chi1_i_seqs,
      af::tiny<unsigned, 4> const& chi2_i_seqs,
      af::tiny<unsigned, 4> const& chi3_i_seqs,
      af::tiny<unsigned, 4> const& chi4_i_seqs,
      size_t residue_index_=1)
    :
      residue_name(residue_name_),
      n_angles(n_angles_),
      have_phi_psi(false),
      residue_index(residue_index_)
    {
      MMTBX_ASSERT(residue_index > 0);
      initialize_i_seqs(chi1_i_seqs, chi2_i_seqs, chi3_i_seqs, chi4_i_seqs);
    }

    size_t first_i_seq () {
      return chi_i_seqs[0];
    }

    void initialize_i_seqs (
      af::tiny<unsigned, 4> const& chi1_i_seqs,
      af::tiny<unsigned, 4> const& chi2_i_seqs,
      af::tiny<unsigned, 4> const& chi3_i_seqs,
      af::tiny<unsigned, 4> const& chi4_i_seqs)
    {
      for (unsigned i = 0; i < 4; i++) {
        chi_i_seqs[i] = chi1_i_seqs[i];
        chi_i_seqs[i+4] = chi2_i_seqs[i];
        chi_i_seqs[i+8] = chi3_i_seqs[i];
        chi_i_seqs[i+12] = chi4_i_seqs[i];
      }
    }

    double get_rotamer_rmsd (
      af::tiny<double, 4> const& angles,
      af::const_ref<scitbx::vec3<double> > const& sites_cart)
    {
      MMTBX_ASSERT(n_angles > 0);
      using cctbx::geometry_restraints::dihedral;
      double rmsd = 0.0;
      for (unsigned j = 0; j < n_angles; j++) {
        af::tiny<scitbx::vec3<double>, 4> chi_sites;
        for (unsigned k = 0; k < 4; k++) {
          chi_sites[k] = sites_cart[chi_i_seqs[(j*4)+k]];
        }
        dihedral chi(chi_sites, 0, 1.0);
        double angle_rad = angles[j] * PI_OVER_180;
        double chi_rad = chi.angle_model * PI_OVER_180;
        rmsd += std::pow(angle_rad - chi_rad, 2);
      }
      return rmsd / n_angles;
    }

    int find_dihedral_proxy (
      dihedral_proxy const& dihedral_proxy)
    {
      MMTBX_ASSERT(n_angles > 0);
      for (unsigned j = 0; j < n_angles; j++) {
        af::tiny<unsigned, 4> j_seqs;
        af::tiny<unsigned, 4> i_seqs = dihedral_proxy.i_seqs;
        for (unsigned k = 0; k < 4; k++) {
          j_seqs[k] = chi_i_seqs[(j*4)+k];
        }
        if ((j_seqs[0] == i_seqs[0]) && (j_seqs[3] == i_seqs[3])) {
          // order: 0 1 2 3
          if ((j_seqs[1] == i_seqs[1]) && (j_seqs[2] == i_seqs[2])) {
            return j + 1;
          // order: 0 2 1 3
          } else if ((j_seqs[1] == i_seqs[2]) && (j_seqs[2] == i_seqs[1])) {
            return - (j + 1);
          }
        } else if ((j_seqs[0] == i_seqs[3]) && (j_seqs[3] == i_seqs[0])) {
          // order: 3 1 2 0
          if ((j_seqs[1] == i_seqs[1]) && (j_seqs[2] == i_seqs[2])) {
            return j + 1;
          // order: 3 2 1 0
          } else if ((j_seqs[1] == i_seqs[2]) && (j_seqs[2] == i_seqs[1])) {
            return - (j + 1);
          }
        }
      }
      return 0;
    }

  };

}} // namespace mmtbx::geometry_restraints

#endif
