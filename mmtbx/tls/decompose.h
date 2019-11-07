#ifndef MMTBX_TLS_DECOMPOSE_H
#define MMTBX_TLS_DECOMPOSE_H

#include <string>

// Basic data types
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/sym_mat3.h>

// Allow arrays of the above
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/matrix/eigensystem.h>

namespace mmtbx { namespace tls { namespace decompose {

namespace af = scitbx::af;

typedef scitbx::matrix::eigensystem::real_symmetric<double> Eig;
typedef af::shared<double>                  eigenvalues;
typedef af::versa<double, af::c_grid<2> >   eigenvectors;

// Actual validation functions
class decompose_tls_matrices {

  public:
    // Constructor
    decompose_tls_matrices(scitbx::sym_mat3<double> const& T,
                           scitbx::sym_mat3<double> const& L,
                           scitbx::mat3<double> const& S,
                           bool l_and_s_in_degrees = true,
                           bool verbose = false,
                           double tol = 1.e-6,
                           double eps = 1.e-8,
                           std::string const& t_S_formula = "11",
                           double t_S_value = 0.0);

    bool is_valid();
    bool is_verbose();
    std::string error();

    // Tolerance of the physical parameters
    double precision_tolerance();
    // Tolerance of the numerical parameters
    double floating_point_limit();

    // --------------------
    // Output variables (set throughout and grouped here for clarity)
    // --------------------
    // RMS amplitudes of libration
    af::shared<double> l_amplitudes;
    // Libration axes ... in the M basis
    scitbx::vec3<double> l1_M, l2_M, l3_M;
    // Points on libration axes ... in M basis
    scitbx::vec3<double> w1_M, w2_M, w3_M;
    // Amplitudes of screw motions (along librational axes)
    af::shared<double> s_amplitudes;
    // RMS amplitudes of vibration
    af::shared<double> v_amplitudes;
    // Vibration axes ... in the L basis
    scitbx::vec3<double> v1_L, v2_L, v3_L;
    // Vibration axes ... in the M basis
    scitbx::vec3<double> v1_M, v2_M, v3_M;

    // --------------------
    // Coordinate transformation
    // --------------------
    // Rotation matrix to transform between original and libration coordinate frames (M -> L)
    // R_ML * a_M = a_L
    scitbx::mat3<double> R_ML;
    scitbx::mat3<double> R_ML_t;
    // Rotation matrix to transform between original and vibrational coordinate frames
    // R_MV * a_M = a_V
    scitbx::mat3<double> R_MV;
    scitbx::mat3<double> R_MV_t;

  private:
    bool verbose_;

    // Response variables
    bool valid_;
    std::string error_;

    // Define analysis precision (for determining positive-definiteness, etc).
    double tol;
    // Define floating point boundary (for determining when something is zero)
    double eps;

    // Value for forcing t_S
    double t_S_value;
    // Select which formula to choose for finding t_S
    std::string t_S_formula;

    // printing
    static std::string spacer;

    // --------------------
    // Decomposition functions
    // --------------------
    void run();
    void stepA();
    void stepB();
    void stepC();
    void stepD();

    // --------------------
    // Input variables
    // --------------------
    // Input matrices
    scitbx::sym_mat3<double> T_M;
    scitbx::sym_mat3<double> L_M;
    scitbx::mat3<double>     S_M;
    // Origin -- only used for manipulating the output axes, does not affect decomposition
    // TODO NOT IMPLEMENTED TODO
    //scitbx::vec3<double> const origin;

    // Matrices at various stages during the decomposition and internal variables
    // --------------------
    // Set at step A
    // --------------------
    // Matrices in the basis of the libration matrix
    scitbx::sym_mat3<double> T_L;
    scitbx::sym_mat3<double> L_L;
    scitbx::mat3<double> S_L;
    // --------------------
    // Set at step B
    // --------------------
    // Store information about the position of the libration axes
    struct W {
      double wy_lx;
      double wz_lx;
      double wz_ly;
      double wx_ly;
      double wx_lz;
      double wy_lz;
      scitbx::vec3<double> w_lx,  w_ly,  w_lz;
    } w_15, w;
    // Translational components from offset of librational axes from origin
    scitbx::sym_mat3<double> D_WL;
    // Translational components in the L-basis
    // T_CL = T_L - D_WL
    scitbx::sym_mat3<double> T_CL;
    // --------------------
    // Set at step C
    // --------------------
    // Amount to subtract from the trace of the S-matrix
    double t_S;
    // S-matrix after subtraction of t_S
    // S_C = S_L - diag(t_S)
    scitbx::mat3<double> S_C;
    // Screw contributions to translation matrix
    scitbx::sym_mat3<double> C_L_t_S;
    // Vibrational component in the L-basis
    // V_L = T_CL - C_L_t_S
    scitbx::sym_mat3<double> V_L;
    // --------------------
    // Set at step D
    // --------------------
    // (Pure) Vibrations in the vibration basis
    scitbx::sym_mat3<double> V_V;

    // --------------------
    // Eigen-analysis functions
    // --------------------
    bool is_pd(scitbx::sym_mat3<double> const& matrix);
    // Quadratic equation coefficients for determinant of V-matrix from T+S matrices
    scitbx::vec3<double> as_bs_cs(
            double t,
            double Sxx, double Syy, double Szz,
            double T11, double T22, double T33, double T12, double T13, double T23
            );

    // --------------------
    // Logic/zero functions
    // --------------------
    // Dealing with small values
    bool is_zero(double value);
    bool is_positive(double value);
    bool is_negative(double value);
    double zero_filter(double value);
    // Set small values to zero
    void zero_small_values(scitbx::vec3<double>& vector);
    void zero_small_values(scitbx::mat3<double>& matrix);
    void zero_small_values(scitbx::sym_mat3<double>& matrix);

    // --------------------
    // Matrix Helper functions
    // --------------------
    // Regular matrices
    double xx(scitbx::mat3<double> const& matrix);
    double xy(scitbx::mat3<double> const& matrix);
    double xz(scitbx::mat3<double> const& matrix);
    double yx(scitbx::mat3<double> const& matrix);
    double yy(scitbx::mat3<double> const& matrix);
    double yz(scitbx::mat3<double> const& matrix);
    double zx(scitbx::mat3<double> const& matrix);
    double zy(scitbx::mat3<double> const& matrix);
    double zz(scitbx::mat3<double> const& matrix);
    // Symmetric matrices
    double xx(scitbx::sym_mat3<double> const& matrix);
    double xy(scitbx::sym_mat3<double> const& matrix);
    double xz(scitbx::sym_mat3<double> const& matrix);
    double yx(scitbx::sym_mat3<double> const& matrix);
    double yy(scitbx::sym_mat3<double> const& matrix);
    double yz(scitbx::sym_mat3<double> const& matrix);
    double zx(scitbx::sym_mat3<double> const& matrix);
    double zy(scitbx::sym_mat3<double> const& matrix);
    double zz(scitbx::sym_mat3<double> const& matrix);

    // Print matrices
    void print(std::string const& label, double value);
    void print(std::string const& label, scitbx::mat3<double> const& matrix);
    void print(std::string const& label, scitbx::sym_mat3<double> const& matrix);
    void print(std::string const& label, scitbx::vec3<double> const& vector);
    void print(std::string const& label, af::shared<double> const& vector);

};

}}} // close mmtbx/tls/decompose

#endif

