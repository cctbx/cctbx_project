#ifndef MMTBX_VALIDATION_RAMACHANDRAN_RAMA_EVAL_H
#define MMTBX_VALIDATION_RAMACHANDRAN_RAMA_EVAL_H

#include <mmtbx/validation/ramachandran/rama8000_tables.h>
#include <iostream>
#include <algorithm>
#include <scitbx/math/linear_interpolation.h>
#include <math.h>


namespace mmtbx { namespace validation { namespace ramachandran {

class rama_eval {
  public:
    rama_eval() {}

    double get_score(int const& rama_class, double const& phi, double const& psi) {
      scitbx::af::const_ref<double, scitbx::af::c_grid<2> > table;
      table = get_rama_table(rama_class);
      double ranged_phi = phi;
      double ranged_psi = psi;
      bring_angle_to_range(ranged_phi);
      bring_angle_to_range(ranged_psi);
      // std::cout << "ranged phi/psi " << phi << "/" << psi << std::endl;
      int lower_bin_phi, higher_bin_phi, lower_bin_psi, higher_bin_psi;
      double lower_value_phi,higher_value_phi,lower_value_psi,higher_value_psi;
      get_bins_and_values(
          ranged_phi, lower_bin_phi, higher_bin_phi,
          lower_value_phi, higher_value_phi);
      get_bins_and_values(
          ranged_psi, lower_bin_psi, higher_bin_psi,
          lower_value_psi, higher_value_psi);
      double result;
      // std::cout << "arguments to interpolation_2d: \n" <<
      //     "  x1, y1: "<< lower_value_phi << ", " << lower_value_psi <<"\n" <<
      //     "  x2, y2: "<< higher_value_phi <<", " << higher_value_psi <<"\n" <<
      //     "  lower bins: " << lower_bin_phi <<", " << lower_bin_psi <<"\n" <<
      //     "  higher bins: " << higher_bin_phi <<", " << higher_bin_psi <<"\n" <<
      //     "  v1,v2,v3,v4: "<< table(lower_bin_phi, lower_bin_psi) << ", " <<
      //     table(higher_bin_phi, higher_bin_psi) << ", " <<
      //     table(lower_bin_phi, higher_bin_psi) << ", " <<
      //     table(higher_bin_phi, lower_bin_psi) << "\n" <<
      //     "  original phi, psi: " << phi << ", " << psi << std::endl;
      result = scitbx::math::linear_interpolation_2d(
          lower_value_phi, lower_value_psi,
          higher_value_phi, higher_value_psi,
          table(lower_bin_phi, lower_bin_psi),
          table(higher_bin_phi, higher_bin_psi),
          table(lower_bin_phi, higher_bin_psi),
          table(higher_bin_phi, lower_bin_psi),
          ranged_phi, ranged_psi);
      return result;
    }

    double
    get_score(std::string const& rama_class, double const& phi, double const& psi) {
      int int_res_type = convert_rama_class(rama_class);
      return get_score(int_res_type, phi, psi);
    }

    int evaluate_angles(std::string const& rama_class, double const& phi, double const& psi) {
      int int_res_type = convert_rama_class(rama_class);
      return evaluate_angles(int_res_type, phi, psi);
    }

    int evaluate_angles(int const& rama_class, double const& phi, double const& psi) {
      double score = get_score(rama_class, phi, psi);
      return evaluate_score(rama_class, score);
    }

    int evaluate_score(std::string const& residue_type, double const& score) {
      int int_res_type = convert_rama_class(residue_type);
      return evaluate_score(int_res_type, score);
    }

    int evaluate_score(int const& residue_type, double const& score) {
      if (score >= 0.02) return RAMALYZE_FAVORED;
      if (residue_type == RAMA_GENERAL) {
        if (score >=0.0005) return RAMALYZE_ALLOWED;
        else return RAMALYZE_OUTLIER;
      }
      if (residue_type == RAMA_CISPRO) {
        if (score >= 0.0020) return RAMALYZE_ALLOWED;
        else return RAMALYZE_OUTLIER;
      }
      if (score >= 0.0010) return RAMALYZE_ALLOWED;
        else return RAMALYZE_OUTLIER;
    }

    int convert_rama_class(std::string const& rama_class_str) {
      return std::distance(
          res_types, std::find(res_types, res_types + 5, rama_class_str));
    }


  protected:

    inline
    void
    bring_angle_to_range(double& value) {
      while (value > 180) value -= 360;
      while (value < -180) value += 360;
    }

    inline
    void
    get_bins_and_values(
        double& value,
        int& lower_bin,
        int& higher_bin,
        double& lower_value,
        double& higher_value)
    {
      lower_value = floor(value);
      if (int(lower_value) % 2 == 0) lower_value -=1;
      higher_value = ceil(value);
      if (int(higher_value) % 2 == 0) higher_value +=1;
      if (lower_value == higher_value) {
        higher_value += 2;
      }
      lower_bin = get_bin_number(lower_value);
      higher_bin = get_bin_number(higher_value);
    }

    inline
    int get_bin_number(double & value) {
      int bin = int((value + 179) / 2.);
      if (bin > 179) bin -= 180;
      if (bin < 0) bin += 180;
      return bin;
    }


    scitbx::af::const_ref<double, scitbx::af::c_grid<2> >
    get_rama_table(int const& rama_class) {
      if (rama_class == RAMA_GENERAL) return table_general;
      else if (rama_class == RAMA_GLYCINE) return table_glycine;
      else if (rama_class == RAMA_CISPRO) return table_cis_pro;
      else if (rama_class == RAMA_TRANSPRO) return table_trans_pro;
      else if (rama_class == RAMA_PREPRO) return table_pre_pro;
      else if (rama_class == RAMA_ILE_VAL) return table_ile_val;
      else {
        char buf[]= "Unknown Ramachandran type.";
        throw std::runtime_error(buf);
      }

    };
};

}}} // namespace mmtbx::validation::ramachandran

#endif // MMTBX_VALIDATION_RAMACHANDRAN_RAMA_EVAL_H
