
#include <simtbx/nanoBragg/nanoBragg.h>
#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>
//#include <boost/python/numpy.hpp>

namespace simtbx {
namespace nanoBragg {

class derivative_manager{
  public:
    derivative_manager();
    void initialize(int sdim, int fdim);
    af::flex_double raw_pixels;
    af::flex_double raw_pixels2;
    double value; // the value of the parameter
    double dI; // the incremental derivative
    double dI2; // the incremental derivative
    bool refine_me;
    bool second_derivatives;
    void increment_image(int idx, double value, double value2);
};


class rot_manager: public derivative_manager{
  public:
    rot_manager();
    virtual ~rot_manager(){}
    virtual void set_R();
    void increment(
        double fudge,
        mat3 X, mat3 Y, mat3 Z,
        mat3 X2, mat3 Y2, mat3 Z2,
        mat3 N, mat3 U, mat3 UBOt,
        vec3 q, vec3 V,
        double Hrad, double Fcell, double Flatt,
        double source_I, double capture_fraction, double omega_pixel);

    mat3 XYZ, XYZ2;
    mat3 R, dR, dR2;
}; // end of rot_manager

class Ncells_manager: public derivative_manager{
  public:
    Ncells_manager();
    virtual ~Ncells_manager(){}
    void increment(
        vec3 V, vec3 H0_vec, vec3 H_vec,
        double Hrad, double Fcell, double Flatt, double fudge,
        double source_I, double capture_fraction, double omega_pixel);
};


class Fcell_manager: public derivative_manager{
  public:
    Fcell_manager();
    virtual ~Fcell_manager(){}
    void increment(
        double Hrad, double Fcell, double Flatt, double fudge,
        double source_I, double capture_fraction, double omega_pixel);
};


class ucell_manager: public derivative_manager{
  public:
    ucell_manager();
    virtual ~ucell_manager(){}
    void increment(
        vec3 V, mat3 NABC, mat3 UR, vec3 q, mat3 Ot,
        double Hrad, double Fcell, double Flatt, double fudge,
        double source_I, double capture_fraction, double omega_pixel);

    mat3 dB, dB2;
};

class origin_manager: public derivative_manager{
  public:
    origin_manager();
    virtual ~origin_manager(){}
    void increment(
        vec3 V, mat3 N, mat3 UBO, vec3 k_diffracted, vec3 o_vec,
        double air_path, double wavelen, double Hrad, double Fcell, double Flatt, double fudge,
        double source_I, double capture_fraction, double omega_pixel, double pixel_size);
    vec3 dk; /* derivative of the diffracted vector along this origin component */
    double FF, FdF, FdF2, dFdF;
};

class rotX_manager: public rot_manager{
  public:
    rotX_manager();
    void set_R();
};
class rotY_manager: public rot_manager{
  public:
    rotY_manager();
    void set_R();
};
class rotZ_manager: public rot_manager{
  public:
    rotZ_manager();
    void set_R();
};

class diffBragg: public nanoBragg{
  public:
  diffBragg(const dxtbx::model::Detector& detector, const dxtbx::model::Beam& beam,
            int verbose, int panel_id);

  ~diffBragg(){};
  void initialize_managers();
  void vectorize_umats();
  void add_diffBragg_spots();
  void init_raw_pixels_roi();
  void zero_raw_pixel_rois();
  void set_ucell_derivative_matrix(int refine_id, af::shared<double> const& value);
  void set_ucell_second_derivative_matrix(int refine_id, af::shared<double> const& value);
  //void reset_derivative_pixels(int refine_id);

  /* methods for interacting with the derivative managers */
  void refine(int refine_id);
  void update_dxtbx_geoms(const dxtbx::model::Detector& detector, const dxtbx::model::Beam& beam,
        int panel_id);
  void set_value( int refine_id, double value);
  double get_value( int refine_id);
  af::flex_double get_derivative_pixels(int refine_id);
  af::flex_double get_second_derivative_pixels(int refine_id);
  af::flex_double get_raw_pixels_roi();

  /* override to cache some of the polarization calc variables to use in derivatives*/
  double polarization_factor(double kahn_factor, double *incident, double *diffracted, double *axis);

  double psi;  /* for polarization correction */
  double u;
  double kEi;
  double kBi;

  //bool use_omega_pixel_ave;
  double om;
  double omega_pixel_ave;
  double airpath_ave;

  double diffracted_ave[4];
  double pixel_pos_ave[4];
  double Fdet_ave, Sdet_ave, Odet_ave;
  vec3 k_diffracted_ave;
  vec3 k_incident_ave;

  mat3 Umatrix;
  mat3 Bmatrix;
  mat3 Omatrix;
  mat3 UBO;
  mat3 Bmat_realspace, NABC;
  mat3 RXYZ;
  std::vector<mat3> RotMats;
  std::vector<mat3> dRotMats, d2RotMats;
  std::vector<mat3> R3, R3_2;

  vec3 k_diffracted;
  vec3 o_vec;
  vec3 Ei_vec, Bi_vec;
  vec3 H_vec, H0_vec;
  vec3 a_vec, ap_vec;
  vec3 b_vec, bp_vec;
  vec3 c_vec, cp_vec;
  vec3 q_vec; // scattering vector

  std::vector<mat3> UMATS;
  std::vector<mat3> UMATS_RXYZ;
  //bool vectorized_umats;

  /* derivative managers */
  std::vector<boost::shared_ptr<rot_manager> > rot_managers;
  std::vector<boost::shared_ptr<ucell_manager> > ucell_managers;
  std::vector<boost::shared_ptr<Ncells_manager> > Ncells_managers;
  std::vector<boost::shared_ptr<origin_manager> > origin_managers;

  boost::shared_ptr<Fcell_manager> fcell_man;
  double* floatimage_roi;
  af::flex_double raw_pixels_roi;

  bool update_oversample_during_refinement;
  bool oversample_omega;
  double uncorrected_I;
  vec3 max_I_hkl;// the hkl corresponding to the maximum intensity in the array (debug)
  //int max_I_h, max_I_k, max_I_l;

  double D;
  double D2;
  double dD;

  // helpful definitions..
  double per_k ;
  double per_k2 ;
  double per_k3 ;
  double per_k4 ;
  double per_k5;
  double per_k6;
  double per_k7;
  double G ;

  double du;
  double du2;

  double w ;
  double w2;
  double BperE2;
  double dkE ;
  double dkB;
  double v ;
  double dv;
  double dpsi ;
  double dpsi2;

  double c2psi ;
  double s2psi ;
  double gam_cos2psi;
  double gam_sin2psi;

  double dpolar ;

  double dpolar2;

  double pp ;
  double dOmega ;
  double dOmega2;

  double FF ;
  double FdF ;
  double dFdF ;
  double FdF2;

  /* om is the average solid angle in the pixel (average over sub pixels) */
  double origin_dI ;
  double origin_dI2;
  // different scale term here because polar and domega terms depend on originZ
  double scale_term2 ;
  double scale_term ;



}; // end of diffBragg

} // end of namespace nanoBragg
} // end of namespace simtbx
