#ifndef MMTBX_TLS_DECOMPOSE_H
#define MMTBX_TLS_DECOMPOSE_H

#include <string>
#include <iostream>

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/sym_mat3.h>

namespace mmtbx { namespace tls { namespace decompose {

// Actual validation functions
class decompose_tls_matrices {

  public:
    // Constructor with origin
    decompose_tls_matrices(scitbx::sym_mat3<double> const& T,
                           scitbx::sym_mat3<double> const& L_deg,
                           scitbx::mat3<double> const& S_deg,
                           scitbx::vec3<double> const& origin); 

    bool is_valid() { return valid_; }
    std::string error() { return error_; }

  private:

    // Response variables 
    bool valid_{true};
    std::string error_{"none"};

    // Define analysis precision (for determining positive-definiteness, etc).
    double tol {1.e-6};
    // Define floating point boundary (for determining when something is zero)
    double eps {1.e-8};

    // printing
    std::string spacer{"   "};

    // --------------------
    // Decomposition functions
    // --------------------
    void run();
    void stepA();
    void stepB();
    void stepC();
    void stepD();
    // Apply origin shifts, etc
    void tidy();

    // --------------------
    // Input variables
    // --------------------
    // Input matrices 
    scitbx::sym_mat3<double> T_M; 
    scitbx::sym_mat3<double> L_M;
    scitbx::mat3<double>     S_M;
    // Origin -- only used for manipulating the output axes, does not affect decomposition
    // TODO NOT IMPLEMENTED TODO
    scitbx::vec3<double> const origin;

    // --------------------
    // Output variables (set throughout and grouped here for clarity)
    // --------------------
    // Amplitudes of libration
    scitbx::vec3<double> l_amplitudes;
    // Libration axes ... in the M basis
    scitbx::vec3<double> l1_M, l2_M, l3_M;
    // Points on libration axes ... in M basis
    scitbx::vec3<double> w1_M, w2_M, w3_M;
    // Amplitudes of screw motions (along librational axes)
    scitbx::vec3<double> s_amplitudes;
    // Amplitudes of vibration
    scitbx::vec3<double> v_amplitudes;
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
        double wy_lx{0.0}, wz_lx{0.0}; 
        double wz_ly{0.0}, wx_ly{0.0};
        double wx_lz{0.0}, wy_lz{0.0};
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
    double t_S{0.0};
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
            double T11, double T22, double T33, double T12, double T13, double T23);

    // --------------------
    // Logic/zero functions
    // --------------------
    // Dealing with small values
    bool is_zero(double value);
    bool is_positive(double value);
    bool is_negative(double value);
    // Return value or zero
    double zero_filter(double value);
    // Set small values to zero
    void zero_small_values(scitbx::vec3<double>& vector);
    void zero_small_values(scitbx::mat3<double>& matrix);
    void zero_small_values(scitbx::sym_mat3<double>& matrix);

    // --------------------
    // Matrix Helper functions 
    // --------------------
    // Normal matrices
    // xx, xy, xz      0, 1, 2
    // yx, yy, yz  ->  3, 4, 5
    // zx, zy, zz      6, 7, 8
    double xx(scitbx::mat3<double> const& matrix) { return matrix[0]; }
    double xy(scitbx::mat3<double> const& matrix) { return matrix[1]; }
    double xz(scitbx::mat3<double> const& matrix) { return matrix[2]; }
    double yx(scitbx::mat3<double> const& matrix) { return matrix[3]; }
    double yy(scitbx::mat3<double> const& matrix) { return matrix[4]; }
    double yz(scitbx::mat3<double> const& matrix) { return matrix[5]; }
    double zx(scitbx::mat3<double> const& matrix) { return matrix[6]; }
    double zy(scitbx::mat3<double> const& matrix) { return matrix[7]; }
    double zz(scitbx::mat3<double> const& matrix) { return matrix[8]; }
    // Symmetric matrices
    // 0, 3, 4
    // -, 1, 5
    // -, -, 2
    double xx(scitbx::sym_mat3<double> const& matrix) { return matrix[0]; }
    double xy(scitbx::sym_mat3<double> const& matrix) { return matrix[3]; }
    double xz(scitbx::sym_mat3<double> const& matrix) { return matrix[4]; }
    double yx(scitbx::sym_mat3<double> const& matrix) { return xy(matrix); }
    double yy(scitbx::sym_mat3<double> const& matrix) { return matrix[1]; }
    double yz(scitbx::sym_mat3<double> const& matrix) { return matrix[5]; }
    double zx(scitbx::sym_mat3<double> const& matrix) { return xz(matrix); }
    double zy(scitbx::sym_mat3<double> const& matrix) { return yz(matrix); }
    double zz(scitbx::sym_mat3<double> const& matrix) { return matrix[2]; }

    // Print matrices 
    void print(std::string const& label, double value);
    void print(std::string const& label, scitbx::mat3<double> const& matrix);
    void print(std::string const& label, scitbx::sym_mat3<double> const& matrix);
    void print(std::string const& label, scitbx::vec3<double> const& vector);
};
                            
}}} // close mmtbx/tls/decompose

#endif

