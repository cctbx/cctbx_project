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
    typedef scitbx::vec3<double> vec3_t;
    typedef scitbx::mat3<double> mat3_t;

  protected:
    af::shared<miller_t> indices_;
    af::shared<double> data_;
    af::shared<double> sigmas_;
    af::shared<int> batch_numbers_or_phases_;
    af::shared<double> wavelengths_;
    af::shared<vec3_t> e_incidents_;
    af::shared<vec3_t> e_scattereds_;

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
        af::const_ref<std::string> const& lines,
        af::const_ref<std::string> const& raw_lines,
        mat3_t const& ort_matrix,
        af::const_ref<int> const& offsets,
        scitbx::sym_mat3<double> const& metrical_matrix)
    {
      mat3_t UB_inv = ort_matrix.inverse();
      scitbx::sym_mat3<double> g_star = metrical_matrix.inverse();
      std::size_t n = lines.size();
      indices_.reserve(n);
      data_.reserve(n);
      sigmas_.reserve(n);
      e_incidents_.reserve(n);
      e_scattereds_.reserve(n);
      for(std::size_t i=0; i<n; i++) {

        //read line from hkl
        std::string line = lines[i];
        miller_t h;
        vec3_t e_inc, e_scat;
        double datum, sigma, frame;
        int run;
        prepare_for_read(line, 40);
        fem::read_from_string(line, "(3i4,2f8.0,i4,f8.0)"),
          h[0], h[1], h[2], datum, sigma, run, frame;
        if (h.is_zero()) break;

        //correct hkl frame number to match raw file
        for (std::size_t i_offset=0; i<run; i++) {
          frame -= offsets[i_offset];
        }

        //find corresponding entry in raw file and store direction cosines and
        //goniometer angles
        vec3_t cosines_inc, cosines_scat;
        double angle1_deg, angle2_deg;
        double omega_deg, chi_deg, phi_deg;
        double omega, chi, phi;
        int scan_axis;
        for(std::size_t i_r=0; i_r<n; i_r++) {
          std::string line_r = raw_lines[i_r];
          miller_t h_r;
          double frame_r;
          double djunk;
          int ijunk;
          prepare_for_read(line_r, 255);

          /*
          fem::read_from_string(
              line_r, "(3i4,20a,6f8.0,17a,f8.0,47a,f7.0,i2,74a,2f8.0)"),
              h_r[0], h_r[1], h_r[2], junk1, cosines_inc[0], cosines_inc[1],
              cosines_inc[2], cosines_scat[0], cosines_scat[1], cosines_scat[2],
              junk2, frame_r, junk3, angle1_deg, scan_axis, junk4, chi_deg, 
              angle2_deg;
              */

          std::string fstring =
            std::string("(3i4,")     + // h_r
            "2f8.0,"    + // i and sigma
            "i4,"       + // batch number
            "6f8.0,"    + // cosines_inc and cosines_scat
            "i2,"       + // status mask
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
              

          fem::read_from_string(line_r, fstring.c_str()),
            h_r[0], h_r[1], h_r[2],  djunk, djunk, ijunk, cosines_inc[0], 
            cosines_inc[1], cosines_inc[2], cosines_scat[0], cosines_scat[1], 
            cosines_scat[2], ijunk, djunk, djunk, frame_r, djunk, djunk, djunk, 
            djunk, djunk, djunk, djunk, angle1_deg, scan_axis, ijunk, ijunk, 
            ijunk, djunk, ijunk, djunk,  ijunk, ijunk, djunk, djunk, djunk, 
            chi_deg, angle2_deg, ijunk;



          if ( h==h_r && abs(frame-frame_r)<.1 ) {
            break;
          }

          //std::cout<<"No match for reflection "<<h[0]<<" "<<h[1]<<" "<<h[2]<<std::endl;
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
        else {
          std::cout<<"What a nice day for scan axis "<<scan_axis<<std::endl;
          throw "Scan axis in .raw file should be 2 (omega) or 3 (phi)";
        }
        chi = scitbx::deg_as_rad(chi_deg); 
        omega = scitbx::deg_as_rad(omega_deg); 
        phi = scitbx::deg_as_rad(phi_deg); 


        vec3_t hkl_zaxis = hkl_z(UB_inv, phi, chi);

        e_inc = make_perpendicular(hkl_zaxis, cosines_inc, g_star);
        e_inc = e_inc.normalize();
        e_scat = make_perpendicular(hkl_zaxis, cosines_scat, g_star);
        e_scat = e_scat.normalize();



      }
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
        double const &chi) {
      //Compute the hkl indices of the laboratory z axis. Ort_matrix is similar
      //to the Busing/Levy UB-matrix but conforms to the Bruker convention in
      //which chi is a rotation around the Y axis instead of X.

      double sin_phi = sin(phi);
      double cos_phi = cos(phi);
      double sin_chi = sin(chi);
      double cos_chi = cos(chi);
      vec3_t v = vec3_t(
          sin(phi) * sin(chi),
          -1 * cos(phi) * sin(chi),
          cos(chi));
      vec3_t h_z = UB_inv * v;
      return h_z;
    }


    af::shared<miller_t> indices() { return indices_; };

    af::shared<double> data() { return data_; }

    af::shared<double> sigmas() { return sigmas_; }

    af::shared<int> alphas() { return batch_numbers_or_phases_; }

    af::shared<int> batch_numbers() { return batch_numbers_or_phases_; }

    af::shared<double> wavelengths() { return wavelengths_; }

    af::shared<vec3_t> e_incidents() { return e_incidents_; }

    af::shared<vec3_t> e_scattereds() { return e_scattereds_; }
};

}} // iotbx::shelx

#endif // GUARD
