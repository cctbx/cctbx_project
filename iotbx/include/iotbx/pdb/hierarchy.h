#ifndef IOTBX_PDB_HIERARCHY_H
#define IOTBX_PDB_HIERARCHY_H

#include <iotbx/pdb/small_str.h>
#include <scitbx/sym_mat3.h>
#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>
#include <map>
#include <string>

namespace iotbx {

//! Handling of files in PDB format.
namespace pdb {

  namespace af = scitbx::af;

  typedef scitbx::vec3<double> vec3;
  typedef scitbx::sym_mat3<double> sym_mat3;
  typedef boost::int32_t int32_t;

  struct atom
  {
    str4 name;
    str4 segid;
    str2 element;
    str2 charge;
    vec3 xyz;
    vec3 sigxyz;
    double occ;
    double sigocc;
    double b;
    double sigb;
    sym_mat3 uij;
    sym_mat3 siguij;

    atom() {}

    atom(
      str4 name_, str4 segid_, str2 element_, str2 charge_,
      vec3 const& xyz_, vec3 const& sigxyz_,
      double occ_, double sigocc_,
      double b_, double sigb_,
      sym_mat3 const& uij_, sym_mat3 const& siguij_)
    :
      name(name_), segid(segid_), element(element_), charge(charge_),
      xyz(xyz_), sigxyz(sigxyz_),
      occ(occ_), sigocc(sigocc_),
      b(b_), sigb(sigb_),
      uij(uij_), siguij(siguij_)
    {}
  };

  typedef boost::shared_ptr<atom> atom_shptr;

  struct residue
  {
    str4 name;
    int32_t seq;
    str4 icode;
    typedef std::map<str4, atom_shptr> atoms_t;
    atoms_t atoms;
    bool hetero;

    residue() {}

    residue(residue const& other)
    :
      name(other.name),
      seq(other.seq),
      icode(other.icode),
      hetero(other.hetero)
    {
      typedef atoms_t::const_iterator aci;
      aci o_end = other.atoms.end();
      for(aci o=other.atoms.begin();o!=o_end;o++) {
        atoms[o->first] = atom_shptr(new atom(*(o->second)));
      }
    }
  };

  struct conformer
  {
    std::string id;
    std::vector<residue> residues;
  };

  struct chain
  {
    std::string id;
    std::vector<conformer> conformers;
  };

  struct model
  {
    int id;
    std::vector<chain> chains;
  };

}} // namespace iotbx::pdb

#endif // IOTBX_PDB_HIERARCHY_H
