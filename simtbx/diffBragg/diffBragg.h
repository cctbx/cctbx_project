
#include <simtbx/nanoBragg/nanoBragg.h>
#include <vector>

namespace simtbx {
namespace nanoBragg {

class derivative_manager{
  public:
    derivative_manager();
    void initialize(int sdim, int fdim);

    double* floatimage;
    af::flex_double raw_pixels;
    double value; // the value of the parameter
    double dI; // the incremental derivative
};

class rot_manager: public derivative_manager{
  public:
    rot_manager();
    double compute_increment(
        int Na, int Nb, int Nc,
        double hfrac, double kfrac, double lfrac,
        double fudge,
        mat3 U, mat3 A, mat3 B, mat3 C,
        vec3 a, vec3 b, vec3 c,  vec3 q,
        double Hrad, double Fcell, double Flatt,
        double source_I, double capture_fraction, double omega_pixel);

    mat3 XYZ;
}; // end of rot_manager

class diffBragg: public nanoBragg{
  public:
  diffBragg(const dxtbx::model::Detector& detector, const dxtbx::model::Beam& beam,
            int verbose, int panel_id);

  void initialize_managers();
  void vectorize_umats();
  void add_diffBragg_spots();

  mat3 RX, RY, RZ, RXYZ;
  mat3 dRX, dRY, dRZ;
  vec3 a_vec, ap_vec;
  vec3 b_vec, bp_vec;
  vec3 c_vec, cp_vec;
  vec3 q_vec; // scattering vector

  std::vector<mat3> UMATS;
  std::vector<mat3> UMATS_RXYZ;

  /* derivative managers */
  rot_manager rotX_man;
  rot_manager rotY_man;
  rot_manager rotZ_man;

}; // end of diffBragg

} // end of namespace nanoBragg
} // end of namespace simtbx
