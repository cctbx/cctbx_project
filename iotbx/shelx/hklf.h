#ifndef IOTBX_SHELX_HKLF_H
#define IOTBX_SHELX_HKLF_H

#include <cctbx/miller.h>
#include <scitbx/array_family/shared.h>
#include <fable/fem/read.hpp>
#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/sym_mat3.h>
#include <math.h>

namespace iotbx { namespace shelx {

namespace af = scitbx::af;

class hklf_reader
{
  public:
    typedef char char_t;
    typedef cctbx::miller::index<> miller_t;
    typedef cctbx::miller::index<double> vec3_t;
    typedef scitbx::mat3<double> mat3_t;

  protected:
    af::shared<miller_t> indices_;
    af::shared<double> data_;
    af::shared<double> sigmas_;
    af::shared<int> batch_numbers_or_phases_;
    af::shared<double> wavelengths_;
    af::shared<vec3_t> u_incs_;
    af::shared<vec3_t> u_scats_;
    af::shared<vec3_t> v_scats_;
    af::shared<double> pol_factors_;

    static
    void
    prepare_for_read(
      std::string& line,
      std::size_t target_size)
    {
      std::size_t initial_size = line.size();
      for(std::size_t i=0;i<initial_size;i++) {
        if (line[i] < ' ') line[i] = ' '; // emulates shelxl
      }
      if (initial_size < target_size) {
        line.append(target_size-initial_size, ' ');
      }
    }

    static
    bool
    substr_is_whitespace_only(
      std::string const& line,
      std::size_t start,
      std::size_t stop)
    {
      for(std::size_t i=start;i<stop;i++) {
        if (!fem::utils::is_whitespace(line[i])) return false;
      }
      return true;
    }

  public:
    hklf_reader(
      af::const_ref<std::string> const& lines)
    {
      std::size_t n = lines.size();
      indices_.reserve(n);
      data_.reserve(n);
      sigmas_.reserve(n);
      batch_numbers_or_phases_.reserve(n);
      wavelengths_.reserve(n);
      bool have_bp = false;
      bool have_wa = false;
      for(std::size_t i=0;i<n;i++) {
        std::string line = lines[i];
        miller_t h;
        double datum, sigma, wa;
        int bp;
        prepare_for_read(line, 40);
        fem::read_from_string(line, "(3i4,2f8.0,i4,f8.4)"),
          h[0], h[1], h[2], datum, sigma, bp, wa;
        if (h.is_zero()) break;
        indices_.push_back(h);
        data_.push_back(datum);
        sigmas_.push_back(sigma);
        batch_numbers_or_phases_.push_back(bp);
        wavelengths_.push_back(wa);
        if (!have_bp) have_bp = !substr_is_whitespace_only(line, 28, 32);
        if (!have_wa) have_wa = !substr_is_whitespace_only(line, 32, 40);
      }
      if (indices_.size() == 0) {
        throw std::runtime_error("No data in SHELX hklf file.");
      }
      if (!have_bp) batch_numbers_or_phases_ = af::shared<int>();// free memory
      if (!have_wa) wavelengths_ = af::shared<double>();// free memory
    }

    hklf_reader(
        af::const_ref<std::string> const& hkl_lines,
        af::const_ref<std::string> const& raw_lines,
        mat3_t const& ort_matrix,
        af::const_ref<int> const& offsets,
        scitbx::sym_mat3<double> const& metrical_matrix)
    {
      mat3_t UB_inv = ort_matrix.inverse();
      double len_a = sqrt(metrical_matrix[0]);
      double len_b = sqrt(metrical_matrix[1]);
      double len_c = sqrt(metrical_matrix[2]);
      scitbx::sym_mat3<double> g_star = metrical_matrix.inverse();
      std::size_t n = hkl_lines.size();
      indices_.reserve(n);
      data_.reserve(n);
      sigmas_.reserve(n);
      u_incs_.reserve(n);
      u_scats_.reserve(n);

      std::string raw_fstring;
      if (1) {          //just to fold this huge string
        raw_fstring =
          std::string("(3i4,") + // h_r
          "2f8.0,"    + // i and sigma
          "i4,"       + // batch number
          "6f8.0,"    + // cosines_inc and cosines_scat
          "i3,"       + // status mask
          "2f7.0,"    + // detector coordinates
          "f8.0,"     + // frame number of refln centroid
          "2f7.0,"    + // detector coordinates (reduced)
          "f8.0,"     + // frame number (predicted)
          "f6.0,"     + // LP factor
          "f5.0,"     + // profile correlation
          "3f7.0,"    + // cumul. exp. time, det. angle, scan axis angle
                        // (omega for omega scan, phi for phi scan)
          "i2,"       + // scan axis: 2 for omega, 3 for phi
          "i5,"       + // 10000*sin(th)/lambda
          "i9,"       + // total counts
          "i7,"       + // background
          "f7.2,"     + // scan axis relat. to start of run
          "i4,"       + // crystal number
          "f6.0,"     + // fraction of peak vol. nearest this hkl
          "i11,"      + // sort key
          "i3,"       + // point-group multiplicity of refln
          "f6.0,"     + // Lorentz factor
          "4f8.0,"    + // corrected det. coordinates, chi, "other angle"
                        // (phi for omega scan, omega for phi scan)
          "i3)";        // twin comp. number
      }

      //loop over hkl_lines
      for(std::size_t i=0; i<n; i++) {

        //read line from hkl
        std::string line = hkl_lines[i];
        miller_t h;
        vec3_t e_inc, e_scat;
        double datum, sigma, frame;
        int run;
        prepare_for_read(line, 40);
        fem::read_from_string(line, "(3i4,2f8.0,i4,f8.0)"),
          h[0], h[1], h[2], datum, sigma, run, frame;
        if (h.is_zero()) continue;

        //correct hkl frame number to match raw file
        frame -= offsets[run-1];

        //find corresponding entry in raw file and store direction cosines and
        //goniometer angles
        vec3_t cosines_inc, cosines_scat;
        miller_t h_r;
        double angle1_deg, angle2_deg, omega_deg, chi_deg, phi_deg;
        double omega, chi, phi;
        double frame_r;
        int scan_axis;
        double djunk;
        int ijunk;

        // hkl_lines and raw_lines are both sorted by hkl before we receive them
        // so we can expect to find the matching raw line around i_r=i. We back
        // up a few lines because we haven't performed any sorting by frame #
        // within groups of identical hkls.
        for (
            std::size_t i_r=(i<5 ? 0 : i-5);
            i_r<raw_lines.size();
            i_r++)
        {
          std::string line_r = raw_lines[i_r];
          prepare_for_read(line_r, 255);

          //std::cout<<"\n"<<line<<"\n\n"<<line_r<<std::endl;
          fem::read_from_string(line_r, raw_fstring.c_str()),
            h_r[0], h_r[1], h_r[2],  djunk, djunk, ijunk, cosines_inc[0],
            cosines_scat[0], cosines_inc[1], cosines_scat[1], cosines_inc[2],
            cosines_scat[2], ijunk, djunk, djunk, frame_r, djunk, djunk, djunk,
            djunk, djunk, djunk, djunk, angle1_deg, scan_axis, ijunk, ijunk,
            ijunk, djunk, ijunk, djunk,  ijunk, ijunk, djunk, djunk, djunk,
            chi_deg, angle2_deg, ijunk;

          //found match
          if ( h==h_r && abs(frame-frame_r)<.1 ) {
            break;
          }

          //crash if we didn't find a matching .raw entry
          if (i_r==raw_lines.size()-1) {
            throw std::runtime_error("No matching .raw entry for reflection");

          }
        }

        //Sort out angles and convert to radians
        if (scan_axis==2) {
          omega_deg = angle1_deg;
          phi_deg = angle2_deg;
        }
        else if (scan_axis==3) {
          phi_deg = angle1_deg;
          omega_deg = angle2_deg;
        }
        chi = scitbx::deg_as_rad(chi_deg);
        omega = scitbx::deg_as_rad(omega_deg);
        phi = scitbx::deg_as_rad(phi_deg);

        //compute the hkl indices of the inc and scattered beams
        vec3_t h_inc(
            cosines_inc[0] * len_a,
            cosines_inc[1] * len_b,
            cosines_inc[2] * len_c);
        vec3_t h_scat(
            cosines_scat[0] * len_a,
            cosines_scat[1] * len_b,
            cosines_scat[2] * len_c);

        //compute the polarization vectors
        vec3_t hkl_zaxis = hkl_z(UB_inv, phi, chi, omega);
        vec3_t u_inc = make_perpendicular(hkl_zaxis, h_inc, g_star);
        vec3_t u_scat = plane_projection(u_inc, h_scat, g_star);
        vec3_t v_scat = make_perpendicular(u_scat, h_scat, g_star);
        //std::cout<<"u_inc:"<<' '<<u_inc[0]<<' '<<u_inc[1]<<' '<<u_inc[2]<<std::endl;
        //std::cout<<"u_scat:"<<' '<<u_scat[0]<<' '<<u_scat[1]<<' '<<u_scat[2]<<std::endl;
        //std::cout<<"v_scat:"<<' '<<v_scat[0]<<' '<<v_scat[1]<<' '<<v_scat[2]<<std::endl;

        u_inc = length1(u_inc, g_star);
        u_scat = length1(u_scat, g_star);
        v_scat = length1(v_scat, g_star);
        double pol_factor = u_inc * g_star * u_scat;
        //std::cout<<"Normalizing..."<<std::endl;
        //std::cout<<"u_inc:"<<' '<<u_inc[0]<<' '<<u_inc[1]<<' '<<u_inc[2]<<std::endl;
        //std::cout<<"u_scat:"<<' '<<u_scat[0]<<' '<<u_scat[1]<<' '<<u_scat[2]<<std::endl;
        //std::cout<<"v_scat:"<<' '<<v_scat[0]<<' '<<v_scat[1]<<' '<<v_scat[2]<<std::endl;
        //std::cout<<"pol_factor: "<<pol_factor<<std::endl;

        double pol2 = u_inc * g_star * v_scat;
        //std::cout<<"pol factor, second beam: "<<pol2<<std::endl;

        //done
        indices_.push_back(h);
        data_.push_back(datum);
        sigmas_.push_back(sigma);
        u_incs_.push_back(u_inc);
        u_scats_.push_back(u_scat);
        v_scats_.push_back(v_scat);
        pol_factors_.push_back(pol_factor);
      }
    }

    static
    vec3_t plane_projection(
        vec3_t const &v1, vec3_t const &v2, scitbx::sym_mat3<double> const &g) {
      //Construct the projection of v1 on the plane normal to v2. g is the
      //metric tensor.

      vec3_t v1_on_v2 = (v1 * g * length1(v2, g)) * v2.normalize();
      return v1 - v1_on_v2;
    }

    static
    vec3_t length1(vec3_t const &v1, scitbx::sym_mat3<double> const &g) {
      //Construct a vector along v1 that has length 1 in the system with metric
      //tensor g.
      return v1 / sqrt(v1 * g * v1);
    }

    static
    vec3_t make_perpendicular(
        vec3_t const &v1, vec3_t const &v2, scitbx::sym_mat3<double> const &g) {
      //Construct a unit vector normal to v1 and v2 in a non-orthogonal system
      //described by the metric tensor g.

      vec3_t c = v1.cross(v2);
      mat3_t K = mat3_t(
          v1[0], v2[0], c[0],
          v1[1], v2[1], c[1],
          v1[2], v2[2], c[2]);
      vec3_t perp = vec3_t(0,0,1) * K.inverse() * g.inverse();
      return perp.normalize();

    }


    vec3_t hkl_z(
        mat3_t const &UB_inv,
        double const &phi,
        double const &chi,
        double const &om)
    {
      //Compute the hkl indices of the laboratory z axis. Ort_matrix is similar
      //to the Busing/Levy UB-matrix but conforms to the Bruker convention in
      //which chi is a rotation around the Y axis instead of X.
      const mat3_t r_phi_inv = mat3_t(
          cos(phi),     sin(phi),     0,
          -1*sin(phi),  cos(phi),     0,
          0,            0,            1
          ).inverse();

      const mat3_t r_chi_inv = mat3_t(
          1,            0,            0,
          0,            cos(chi),     sin(chi),
          0,            -1*sin(chi),  cos(chi)
          ).inverse();

      const mat3_t r_omega_inv = mat3_t(
          cos(om),      -1*sin(om),   0,
          sin(om),      cos(om),      0,
          0,            0,            1
          ).inverse();

      vec3_t v_l(0, 1, 0);

      vec3_t v_h = UB_inv * r_phi_inv * r_chi_inv * r_omega_inv * v_l;

      /*
      vec3_t v = vec3_t(
          sin(phi) * sin(chi),
          -1 * cos(phi) * sin(chi),
          cos(chi));
      vec3_t h_z = UB_inv * v;
      */

      return v_h.normalize();
    }


    af::shared<miller_t> indices() { return indices_; };

    af::shared<double> data() { return data_; }

    af::shared<double> sigmas() { return sigmas_; }

    af::shared<int> alphas() { return batch_numbers_or_phases_; }

    af::shared<int> batch_numbers() { return batch_numbers_or_phases_; }

    af::shared<double> wavelengths() { return wavelengths_; }

    af::shared<vec3_t> u_incs() { return u_incs_; }

    af::shared<vec3_t> u_scats() { return u_scats_; }

    af::shared<vec3_t> v_scats() { return v_scats_; }

    af::shared<double> pol_factors() { return pol_factors_; }
};

}} // iotbx::shelx

#endif // GUARD
