#ifndef MMTBX_VALIDATION_RAMACHANDRAN_RAMA_EVAL_H
#define MMTBX_VALIDATION_RAMACHANDRAN_RAMA_EVAL_H

#include <mmtbx/validation/ramachandran/rama8000_tables.h>
#include <iostream>
#include <scitbx/math/linear_interpolation.h>
#include <math.h>


namespace mmtbx { namespace validation { namespace ramachandran {

class rama_eval {
  public:
    rama_eval() {}

    double
    get_value(std::string const& rama_class, double phi, double psi) {
      scitbx::af::const_ref<double, scitbx::af::c_grid<2> > table;
      table = get_rama_table(rama_class);
      bring_angle_to_range(phi);
      bring_angle_to_range(psi);
      // std::cout << "ranged phi/psi " << phi << "/" << psi << std::endl;
      int lower_bin_phi, higher_bin_phi, lower_bin_psi, higher_bin_psi;
      double lower_value_phi,higher_value_phi,lower_value_psi,higher_value_psi;
      get_bins_and_values(
          phi, lower_bin_phi, higher_bin_phi,
          lower_value_phi, higher_value_phi);
      get_bins_and_values(
          psi, lower_bin_psi, higher_bin_psi,
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
          phi, psi);
      return result;
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
    get_rama_table(std::string const& rama_class) {
      if (rama_class == "general") return table_general;
      else if (rama_class == "glycine") return table_glycine;
      else if (rama_class == "cis-proline") return table_cis_pro;
      else if (rama_class == "trans-proline") return table_trans_pro;
      else if (rama_class == "pre-proline") return table_pre_pro;
      else if (rama_class == "isoleucine or valine") return table_ile_val;
      else {
        char buf[]= "Unknown Ramachandran type.";
        throw std::runtime_error(buf);
      }

    };
};

}}} // namespace mmtbx::validation::ramachandran

#endif // MMTBX_VALIDATION_RAMACHANDRAN_RAMA_EVAL_H
