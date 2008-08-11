#ifndef CCTBX_GEOMETRY_RESTRAINTS_MOTIF_H
#define CCTBX_GEOMETRY_RESTRAINTS_MOTIF_H

#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared.h>
#include <string>

namespace cctbx { namespace geometry_restraints {

  struct motif {

    struct atom
    {
      std::string name;
      std::string scattering_type;
      std::string nonbonded_type;
      double partial_charge;

      atom() : partial_charge(0) {}

      atom(
        const char* name_,
        const char* scattering_type_="",
        const char* nonbonded_type_="",
        double partial_charge_=0)
      :
        name(name_),
        scattering_type(scattering_type_),
        nonbonded_type(nonbonded_type_),
        partial_charge(partial_charge_)
      {}
    };

    struct bond
    {
      af::tiny<std::string, 2> atom_names;
      std::string type;
      double distance_ideal;
      double weight;
      std::string id;

      bond() : distance_ideal(0), weight(0) {}

      bond(
        af::tiny<std::string, 2> atom_names_,
        const char* type_="",
        double distance_ideal_=0,
        double weight_=0,
        const char* id_="")
      :
        atom_names(atom_names_),
        type(type_),
        distance_ideal(distance_ideal_),
        weight(weight_),
        id(id_)
      {}
    };

    struct angle
    {
      af::tiny<std::string, 3> atom_names;
      double angle_ideal;
      double weight;
      std::string id;

      angle() : angle_ideal(0), weight(0) {}

      angle(
        af::tiny<std::string, 3> atom_names_,
        double angle_ideal_=0,
        double weight_=0,
        const char* id_="")
      :
        atom_names(atom_names_),
        angle_ideal(angle_ideal_),
        weight(weight_),
        id(id_)
      {}
    };

    struct dihedral
    {
      af::tiny<std::string, 4> atom_names;
      double angle_ideal;
      double weight;
      int periodicity;
      std::string id;

      dihedral() : angle_ideal(0), weight(0), periodicity(0) {}

      dihedral(
        af::tiny<std::string, 4> atom_names_,
        double angle_ideal_=0,
        double weight_=0,
        int periodicity_=0,
        const char* id_="")
      :
        atom_names(atom_names_),
        angle_ideal(angle_ideal_),
        weight(weight_),
        periodicity(periodicity_),
        id(id_)
      {}
    };

    struct chirality
    {
      af::tiny<std::string, 4> atom_names;
      std::string volume_sign;
      bool both_signs;
      double volume_ideal;
      double weight;
      std::string id;

      chirality() : both_signs(false), volume_ideal(0), weight(0) {}

      chirality(
        af::tiny<std::string, 4> atom_names_,
        const char* volume_sign_="",
        bool both_signs_=false,
        double volume_ideal_=0,
        double weight_=0,
        const char* id_="")
      :
        atom_names(atom_names_),
        volume_sign(volume_sign_),
        both_signs(both_signs_),
        volume_ideal(volume_ideal_),
        weight(weight_),
        id(id_)
      {}
    };

    struct planarity
    {
      af::shared<std::string> atom_names;
      af::shared<double> weights;
      std::string id;

      planarity() {}

      planarity(
        af::shared<std::string> const& atom_names_,
        af::shared<double> const& weights_,
        const char* id_="")
      :
        atom_names(atom_names_),
        weights(weights_),
        id(id_)
      {
        CCTBX_ASSERT(weights.size() == atom_names.size());
      }
    };

    std::string id;
    std::string description;
    af::shared<std::string> info;
    af::shared<std::string> manipulation_ids;
    af::shared<atom> atoms;
    af::shared<bond> bonds;
    af::shared<angle> angles;
    af::shared<dihedral> dihedrals;
    af::shared<chirality> chiralities;
    af::shared<planarity> planarities;

    struct alteration
    {
      struct action_type
      {
        action_type(std::string const& description="")
        {
          if      (description == "")       code_ = 0;
          else if (description == "add")    code_ = 1;
          else if (description == "delete") code_ = 2;
          else if (description == "change") code_ = 3;
          else {
            throw std::runtime_error(
              "Unknown cctbx::geometry_restraints::motif::alteration"
              "::action_type: \"" + description + "\"\n"
              "  Possible action types are:"
              " \"\", \"add\", \"delete\", \"change\"");
          }
        }

        int
        code() const { return code_; }

        bool is_add() const { return code_ == 1; }
        bool is_delete() const { return code_ == 2; }
        bool is_change() const { return code_ == 3; }

        std::string
        description() const
        {
          if (is_add())    return "add";
          if (is_delete()) return "delete";
          if (is_change()) return "change";
          return "";
        }

        private:
          int code_;
      };

      struct operand_type
      {
        operand_type(std::string const& description="")
        {
          if      (description == "")          code_ = 0;
          else if (description == "atom")      code_ = 1;
          else if (description == "bond")      code_ = 2;
          else if (description == "angle")     code_ = 3;
          else if (description == "dihedral")  code_ = 4;
          else if (description == "chirality") code_ = 5;
          else if (description == "planarity") code_ = 6;
          else {
            throw std::runtime_error(
              "Unknown cctbx::geometry_restraints::motif::alteration"
              "::operand_type: \"" + description + "\"\n"
              "  Possible operand types are:"
              " \"\", \"atom\", \"bond\", \"angle\","
              " \"dihedral\", \"chirality\", \"planarity\"");
          }
        }

        int
        code() const { return code_; }

        bool is_atom() const { return code_ == 1; }
        bool is_bond() const { return code_ == 2; }
        bool is_angle() const { return code_ == 3; }
        bool is_dihedral() const { return code_ == 4; }
        bool is_chirality() const { return code_ == 5; }
        bool is_planarity() const { return code_ == 6; }

        std::string
        description() const
        {
          if (is_atom())      return "atom";
          if (is_bond())      return "bond";
          if (is_angle())     return "angle";
          if (is_dihedral())  return "dihedral";
          if (is_chirality()) return "chirality";
          if (is_planarity()) return "planarity";
          return "";
        }

        private:
          int code_;
      };

      action_type action;
      operand_type operand;
      af::shared<std::string> motif_ids;
      std::string motif_atom_name;
      motif::atom atom;
      motif::bond bond;
      motif::angle angle;
      motif::dihedral dihedral;
      motif::chirality chirality;
      motif::planarity planarity;
      std::string planarity_motif_id;
      af::shared<action_type> planarity_atom_actions;
      unsigned change_bits;

      alteration(
        std::string const& action_="",
        std::string const& operand_="")
      :
        action(action_),
        operand(operand_),
        change_bits(0)
      {}

      static const unsigned change_partial_charge_bit = 0x00000001U;
      static const unsigned change_distance_ideal_bit = 0x00000002U;
      static const unsigned change_weight_bit =         0x00000004U;
      static const unsigned change_angle_ideal_bit =    0x00000008U;
      static const unsigned change_periodicity_bit =    0x00000010U;
      static const unsigned change_both_signs_bit =     0x00000020U;
      static const unsigned change_volume_ideal_bit =   0x00000040U;

#define CCTBX_GEOMETRY_RESTRAINTS_MOTIF_MANIP_CHNG_BITS_GET_SET(attr) \
      bool \
      change_##attr() const { return change_bits & change_##attr##_bit; } \
\
      alteration& \
      set_change_##attr(bool state) \
      { \
        if (state) change_bits |= change_##attr##_bit; \
        else       change_bits &= ~change_##attr##_bit; \
        return *this; \
      }

      CCTBX_GEOMETRY_RESTRAINTS_MOTIF_MANIP_CHNG_BITS_GET_SET(partial_charge)
      CCTBX_GEOMETRY_RESTRAINTS_MOTIF_MANIP_CHNG_BITS_GET_SET(distance_ideal)
      CCTBX_GEOMETRY_RESTRAINTS_MOTIF_MANIP_CHNG_BITS_GET_SET(weight)
      CCTBX_GEOMETRY_RESTRAINTS_MOTIF_MANIP_CHNG_BITS_GET_SET(angle_ideal)
      CCTBX_GEOMETRY_RESTRAINTS_MOTIF_MANIP_CHNG_BITS_GET_SET(periodicity)
      CCTBX_GEOMETRY_RESTRAINTS_MOTIF_MANIP_CHNG_BITS_GET_SET(both_signs)
      CCTBX_GEOMETRY_RESTRAINTS_MOTIF_MANIP_CHNG_BITS_GET_SET(volume_ideal)
    };

    struct manipulation
    {
      std::string id;
      std::string description;
      af::shared<std::string> info;
      af::shared<alteration> alterations;
    };

  }; // struct motif

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_MOTIF_H
