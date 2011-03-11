#ifndef MMTBX_GEOMETRY_RESTRAINTS_ROTAMER_H
#define MMTBX_GEOMETRY_RESTRAINTS_ROTAMER_H

#include <scitbx/array_family/shared.h>
#include <mmtbx/error.h>

#include <string>

namespace mmtbx { namespace geometry_restraints {
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

  };

}} // namespace mmtbx::geometry_restraints

#endif
