#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/constants.h>

#include "decompose.h"

namespace mmtbx { namespace tls { namespace decompose {

namespace af = scitbx::af;
typedef scitbx::matrix::eigensystem::real_symmetric<double> Eig;
typedef af::shared<double>                  eigenvalues;
typedef af::versa<double, af::c_grid<2> >   eigenvectors;

// Template FloatType
//template <typename FloatType>

// Constructor

decompose_tls_matrices::decompose_tls_matrices(
        scitbx::sym_mat3<double> const& T,
        scitbx::sym_mat3<double> const& L_deg,
        scitbx::mat3<double> const& S_deg,
        scitbx::vec3<double> const& origin = (0.0,0.0,0.0) )
{
    double deg2rad = scitbx::deg_as_rad(1.0);
    double deg2radsq = deg2rad * deg2rad;

    // Convert and store the input matrices
    T_M = T;
    L_M = L_deg * deg2radsq;
    S_M = S_deg * deg2rad;

    std::cout << "--------------------------" << std::endl;
    std::cout << "Input Matrices:" << std::endl;
    print("T", T_M);
    print("L", L_M);
    print("S", S_M);

    // Run decomposition
    run();
}

// linalg functions

bool decompose_tls_matrices::is_pd(scitbx::sym_mat3<double> const& matrix) {
    // Create eigensystem object
    Eig es(matrix);
    // Extract eigenvalues
    af::shared<double> ev = es.values();
    // Check all are positive
    for (int i=0; i<3; i++) {
        if (matrix[i] < (-1.0*tol)) { return false; }
    }
    return true;
}

scitbx::vec3<double> decompose_tls_matrices::as_bs_cs(
        double t,
        double Sxx,
        double Syy,
        double Szz,
        double T11,
        double T22,
        double T33,
        double T12,
        double T13,
        double T23) {

    double xx_ = (t-Sxx)*(t-Sxx) - T11;
    double yy_ = (t-Syy)*(t-Syy) - T22;
    double zz_ = (t-Szz)*(t-Szz) - T33;

    return scitbx::vec3<double>(
            xx_ + yy_ + zz_,
            xx_*yy_ + yy_*zz_ + zz_*xx_ - (T12*T12 + T23*T23 + T13*T13),
            xx_*yy_*zz_ - T23*T23*xx_ - T13*T13*yy_ - T12*T12*zz_ - 2.0*T12*T23*T13);

}

// check functions

bool decompose_tls_matrices::is_zero(double value) {
    if (std::abs(value) < tol) { return true; }
    return false;
}

bool decompose_tls_matrices::is_positive(double value) {
    if ( value > tol ) { return true; }
    return false;
}

bool decompose_tls_matrices::is_negative(double value) {
    if ( value < -tol) { return true; }
    return false;
}

// print functions

void decompose_tls_matrices::print(std::string const& label, double value) {
    std::cout << label << ": " << spacer << value << spacer << std::endl;
}

void decompose_tls_matrices::print(std::string const& label, scitbx::mat3<double> const& matrix) {
    std::cout << std::setprecision(6) << std::fixed << std::showpos;
    std::cout << label << ":" << std::endl;
    std::cout << spacer << matrix[0] << spacer << matrix[1] << spacer << matrix[2] << std::endl;
    std::cout << spacer << matrix[3] << spacer << matrix[4] << spacer << matrix[5] << std::endl;
    std::cout << spacer << matrix[6] << spacer << matrix[7] << spacer << matrix[8] << std::endl;
}

void decompose_tls_matrices::print(std::string const& label, scitbx::sym_mat3<double> const& matrix) {
    std::cout << std::setprecision(6) << std::fixed << std::showpos;
    std::cout << label << ":" << std::endl;
    std::cout << spacer << matrix[0]  << spacer << matrix[3]  << spacer << matrix[4] << std::endl;
    std::cout << spacer << "+-.------" << spacer << matrix[1]  << spacer << matrix[5] << std::endl;
    std::cout << spacer << "+-.------" << spacer << "+-.------" << spacer << matrix[2] << std::endl;
}

void decompose_tls_matrices::print(std::string const& label, scitbx::vec3<double> const& vector) {
    std::cout << std::setprecision(6) << std::fixed << std::showpos;
    std::cout << label << ": ";
    std::cout << spacer << vector[0] << spacer << vector[1] << spacer << vector[2] << std::endl;
}

// Small-value functions

double decompose_tls_matrices::zero_filter(double value) {
    if (std::abs(value) < eps) { return 0.0; }
    return value;
}

void decompose_tls_matrices::zero_small_values(scitbx::mat3<double>& matrix) {
    //Iterate through values
    for (int i=0; i<9; i++) {
        if (std::abs(matrix[i]) < eps) { matrix[i] = 0.0; }
    }
}

void decompose_tls_matrices::zero_small_values(scitbx::sym_mat3<double>& matrix) {
    //Iterate through values
    for (int i=0; i<6; i++) {
        if (std::abs(matrix[i]) < eps) { matrix[i] = 0.0; }
    }
}

void decompose_tls_matrices::zero_small_values(scitbx::vec3<double>& vector) {
    //Iterate through values
    for (int i=0; i<3; i++) {
        if (std::abs(vector[i]) < eps) { vector[i] = 0.0; }
    }
}

// Main functions

/**
 * Diagonalise the L matrix
*/
void decompose_tls_matrices::stepA() {

    std::cout << std::endl;
    std::cout << " ######################################" << std::endl;
    std::cout << " <-------        Step A        ------->" << std::endl;
    std::cout << " ######################################" << std::endl << std::endl;

    // Bail if T not positive definite
    if (! is_pd(T_M)) {
        throw std::runtime_error("Step A: Input matrix T[M] is not positive semidefinite.");
    }
    // Bail if L not positive definite
    if (! is_pd(L_M)) {
        throw std::runtime_error("Step A: Input matrix L[M] is not positive semidefinite.");
    }

    // Find eigenvalues and eigenvectors of L-matrix
    Eig l_eig(L_M);
    eigenvectors l_eig_vecs = l_eig.vectors();
    eigenvalues  l_eig_vals = l_eig.values();

    // Store the eigenvalues and eigenvectors of L
    l_amplitudes = scitbx::vec3<double>(l_eig_vals[0], l_eig_vals[1], l_eig_vals[2]);
    l1_M = scitbx::vec3<double>(&l_eig_vecs[0]);
    l2_M = scitbx::vec3<double>(&l_eig_vecs[3]);
    l3_M = scitbx::vec3<double>(&l_eig_vecs[6]);
    // OVERWRITE - Choose l1 to be cross product of l2 and l3 (guaranteed righthanded system)
    l1_M = l2_M.cross(l3_M);

    // Extract the rotation matrix from the eigenvectors
    // Eigenvectors of L_M are ROWS of R_ML
    // R_ML = ( e1_x, e1_y, e1_z )
    //        ( e2_x, e2_y, e2_z )
    //        ( e3_x, e3_y, e3_z )
    // elements are R_ML(row, column)
    //R_ML = static_cast <scitbx::mat3<double>> (l_eig.vectors().begin());
    R_ML = scitbx::mat3<double>(l1_M[0], l1_M[1], l1_M[2], l2_M[0], l2_M[1], l2_M[2], l3_M[0], l3_M[1], l3_M[2]);
    R_ML_t = R_ML.transpose();

    // Check that determinant is 1 and transpose is inverse
    double det = R_ML.determinant();
    if ( std::abs(det-1.0) > tol ) {
        throw std::runtime_error("Step A: Rotation matrix R[ML] does not have unitary determinant");
    }

    // Print Eigenvalues of L
    std::cout << "--------------------------" << std::endl;
    std::cout << "Eigenvalues of L: " << std::endl;
    print("1", l_amplitudes[0]);
    print("2", l_amplitudes[1]);
    print("3", l_amplitudes[2]);
    std::cout << "Eigenvectors of L: " << std::endl;
    print("1", l1_M);
    print("2", l2_M);
    print("3", l3_M);

    // Print eigenvectors/rotation matrix
    std::cout << "--------------------------" << std::endl;
    print("R[M->L]", R_ML);

    // Check eigenvectors are unitary
    double mod1, mod2, mod3;
    mod1 = R_ML(0,0)*R_ML(0,0) + R_ML(0,1)*R_ML(0,1) + R_ML(0,2)*R_ML(0,2);
    mod2 = R_ML(1,0)*R_ML(1,0) + R_ML(1,1)*R_ML(1,1) + R_ML(1,2)*R_ML(1,2);
    mod3 = R_ML(2,0)*R_ML(2,0) + R_ML(2,1)*R_ML(2,1) + R_ML(2,2)*R_ML(2,2);

    // Sanity checks -- verbose only
    std::cout << "--------------------------" << std::endl;
    std::cout << "Modulus of the eigenvectors of the L_M matrix" << std::endl;
    std::cout << "modulus 1st vector: " << mod1 << std::endl;
    std::cout << "modulus 2nd vector: " << mod2 << std::endl;
    std::cout << "modulus 3rd vector: " << mod3 << std::endl;

    std::cout << "--------------------------" << std::endl;
    print("R[ML] * L[M] * t(R[ML])", R_ML * L_M * R_ML_t);

    // Transform the input matrices to the basis of L_M
    scitbx::mat3<double> T_L_ = R_ML * T_M * R_ML_t;
    scitbx::mat3<double> L_L_ = R_ML * L_M * R_ML_t;
    // Make symmetric and recast
    T_L = scitbx::sym_mat3<double>((T_L_ + T_L_.transpose()) / 2.0);
    L_L = scitbx::sym_mat3<double>((L_L_ + L_L_.transpose()) / 2.0);
    // S is mat3 anyway
    S_L = R_ML * S_M * R_ML_t;

    std::cout << "--------------------------" << std::endl;
    std::cout << "Basis[L] Matrices" << std::endl;
    print("T[L]", T_L);
    print("L[L]", L_L);
    print("S[L]", S_L);

}

/**
 * Find the libration axes in the L-basis
*/
void decompose_tls_matrices::stepB() {

    std::cout << std::endl;
    std::cout << " ######################################" << std::endl;
    std::cout << " <-------        Step B        ------->" << std::endl;
    std::cout << " ######################################" << std::endl << std::endl;

    // Lxx component
    double Lxx{xx(L_L)};
    double Sxy{xy(S_L)};
    double Sxz{xz(S_L)};
    if ( ! is_zero(Lxx) ) {
        w.wz_lx = w_15.wz_lx =   Sxy / Lxx;       // xy is     cyclic permutation >  plus sign
        w.wy_lx = w_15.wy_lx = - Sxz / Lxx;       // xz is anticyclic permutation > minus sign
    } else if ( ! (is_zero(Sxy) and is_zero(Sxz)) ) {
        throw std::runtime_error("Step B: Non-zero off-diagonal S[L] and zero L[L] elements.");
    }

    // Lyy Component
    double Lyy{yy(L_L)};
    double Syz{yz(S_L)};
    double Syx{yx(S_L)};
    if ( ! is_zero(Lyy) ) {
        w.wx_ly = w_15.wx_ly =   Syz / Lyy;       // yz is     cyclic permutation >  plus sign
        w.wz_ly = w_15.wz_ly = - Syx / Lyy;       // yx is anticyclic permutation > minus sign
    } else if ( ! (is_zero(Syz) and is_zero(Syx)) ) {
        throw std::runtime_error("Step B: Non-zero off-diagonal S[L] and zero L[L] elements.");
    }

    // Lzz component
    double Lzz{zz(L_L)};
    double Szx{zx(S_L)};
    double Szy{zy(S_L)};
    if ( ! is_zero(Lzz) ) {
        w.wy_lz = w_15.wy_lz =   Szx / Lzz;       // zx is     cyclic permutation >  plus sign
        w.wx_lz = w_15.wx_lz = - Szy / Lzz;       // zy is anticyclic permutation > minus sign
    } else if ( ! (is_zero(Szx) and is_zero(Szy)) ) {
        throw std::runtime_error("Step B: Non-zero off-diagonal S[L] and zero L[L] elements.");
    }

    // Combine to create vectors
    w_15.w_lx = scitbx::vec3<double>(       0.0, w_15.wy_lx, w_15.wz_lx);
    w_15.w_ly = scitbx::vec3<double>(w_15.wx_ly,        0.0, w_15.wz_ly);
    w_15.w_lz = scitbx::vec3<double>(w_15.wx_lz, w_15.wy_lz,        0.0);

    std::cout << "--------------------------" << std::endl;
    print("w_lx (eq.15)", w_15.w_lx);
    print("w_ly (eq.15)", w_15.w_ly);
    print("w_lz (eq.15)", w_15.w_lz);

    // But what is actually wanted:
    double wxx{(w.wx_ly+w.wx_lz)/2.0};
    double wyy{(w.wy_lz+w.wy_lx)/2.0};
    double wzz{(w.wz_lx+w.wz_ly)/2.0};
    w.w_lx = scitbx::vec3<double>(    wxx, w.wy_lx, w.wz_lx);
    w.w_ly = scitbx::vec3<double>(w.wx_ly,     wyy, w.wz_ly);
    w.w_lz = scitbx::vec3<double>(w.wx_lz, w.wy_lz,     wzz);

    std::cout << "--------------------------" << std::endl;
    print("w_lx (eq.16)", w.w_lx);
    print("w_ly (eq.16)", w.w_ly);
    print("w_lz (eq.16)", w.w_lz);

    // Transform points on libration axes ... in M basis
    w1_M = R_ML_t * w.w_lx;
    w2_M = R_ML_t * w.w_ly;
    w3_M = R_ML_t * w.w_lz;

    std::cout << "--------------------------" << std::endl;
    std::cout << "Points on Libration axis (M-basis)" << std::endl;
    print("w_lx[M]", w1_M);
    print("w_ly[M]", w2_M);
    print("w_lz[M]", w3_M);

    // Calculate D matrix
    double d11 = w.wz_ly * w.wz_ly * Lyy + w.wy_lz * w.wy_lz * Lzz;
    double d22 = w.wx_lz * w.wx_lz * Lzz + w.wz_lx * w.wz_lx * Lxx;
    double d33 = w.wy_lx * w.wy_lx * Lxx + w.wx_ly * w.wx_ly * Lyy;
    double d12 = - w.wx_lz * w.wy_lz * Lzz;
    double d13 = - w.wz_ly * w.wx_ly * Lyy;
    double d23 = - w.wy_lx * w.wz_lx * Lxx;

    // Find the apparent translation caused by offset of libration axes
    D_WL = scitbx::sym_mat3<double>(d11,d22,d33,d12,d13,d23);
    // Subtract from translation matrix
    T_CL = T_L - D_WL;
    // Convert small values to zero in-place
    zero_small_values(T_CL);

    // Check PD
    if (! is_pd(T_CL)) {
        throw std::runtime_error("Step B: Matrix T_C[L] is not positive semidefinite.");
    }

    // Print output matrices from this section
    std::cout << "--------------------------" << std::endl;
    print("D_W[L]", D_WL);
    std::cout << "--------------------------" << std::endl;
    print("T_C[L]", T_CL);

}

void decompose_tls_matrices::stepC() {

    std::cout << std::endl;
    std::cout << " ######################################" << std::endl;
    std::cout << " <-------        Step C        ------->" << std::endl;
    std::cout << " ######################################" << std::endl << std::endl;

    // Have not implemented FORCE_t_S or find_t_S_formula=="10"

    double Lxx{zero_filter(xx(L_L))};
    double Lyy{zero_filter(yy(L_L))};
    double Lzz{zero_filter(zz(L_L))};

    double Sxx{xx(S_L)};
    double Syy{yy(S_L)};
    double Szz{zz(S_L)};

    double T_CLxx = xx(T_CL);
    double T_CLyy = yy(T_CL);
    double T_CLzz = zz(T_CL);

    double t11 = T_CLxx * Lxx;
    double t22 = T_CLyy * Lyy;
    double t33 = T_CLzz * Lzz;

    double rx = std::sqrt(t11);
    double ry = std::sqrt(t22);
    double rz = std::sqrt(t33);

    double t12 = xy(T_CL) * std::sqrt(Lxx*Lyy);
    double t13 = xz(T_CL) * std::sqrt(Lxx*Lzz);
    double t23 = yz(T_CL) * std::sqrt(Lyy*Lzz);

    // Coefficients of the quadratic
    scitbx::vec3<double> abc;

    //
    //  LEFT BRANCH
    //
    // All diagonal components of L are non-zero
    if ( ! (is_zero(Lxx) or is_zero(Lyy) or is_zero(Lzz)) ) {

        std::cout << std::endl << "--------------------------" << std::endl;
        std::cout << "TAKING THE LEFT BRANCH" << std::endl;
        std::cout << "--------------------------" << std::endl << std::endl;

        double t_min_C = std::max({Sxx-rx, Syy-ry, Szz-rz});
        double t_max_C = std::min({Sxx+rx, Syy+ry, Szz+rz});

        if (t_min_C > t_max_C) {
            print("t_min_C", t_min_C);
            print("t_max_C", t_max_C);
            throw std::runtime_error("Step C (left branch): Empty (tmin_c,tmax_c) interval.");
        }

        double num = Sxx * Lyy*Lyy * Lzz*Lzz + \
                     Syy * Lzz*Lzz * Lxx*Lxx + \
                     Szz * Lxx*Lxx * Lyy*Lyy;
        double den = Lyy*Lyy * Lzz*Lzz + \
                     Lzz*Lzz * Lxx*Lxx + \
                     Lxx*Lxx * Lyy*Lyy;
        double t_0 = num/den;

        std::cout << "--------------------------" << std::endl;
        print("t_0", t_0);

        // Construct T_lambda matrix
        scitbx::sym_mat3<double> T_lambda{t11, t22, t33, t12, t13, t23};

        std::cout << "--------------------------" << std::endl;
        print("T_lambda", T_lambda);

        // Calculate eigenvalues & vectors of T_lambda
        Eig es(T_lambda);
        eigenvalues  ev = es.values();

        std::cout << "--------------------------" << std::endl;
        print("Eigenvalues", scitbx::vec3<double>(ev.begin()));

        // Eigenvalues should be in decreasing order
        if ( ! (ev[0] >= ev[1] >= ev[2])) {
            throw std::logic_error("Eigenvalues are not in descending size order.");
        }

        double tau_max = ev[0];

        // Largest eigenvalue must be positive
        if(tau_max < 0.0) {
            throw std::runtime_error("Step C (left branch): Eq.(32): tau_max<0.");
        }

        std::cout << "--------------------------" << std::endl;

        // Print here for cleaness (even though calculated earlier)
        print("t_min_C", t_min_C);
        print("t_max_C", t_max_C);

        double t_min_tau = std::max({Sxx,Syy,Szz})-std::sqrt(tau_max);
        double t_max_tau = std::min({Sxx,Syy,Szz})+std::sqrt(tau_max);

        print("t_min_tau", t_min_tau);
        print("t_max_tau", t_max_tau);

        if(t_min_tau > t_max_tau) {
            throw std::runtime_error("Step C (left branch): Empty (tmin_t,tmax_t) interval.");
        }

        double t_a_sq = t_0*t_0 + (t11+t22+t33)/3.0 - (Sxx*Sxx+Syy*Syy+Szz*Szz)/3.0;
        if(t_a_sq < 0.0) {
            throw std::runtime_error("Step C (left branch): Negative argument when estimating tmin_a.");
        }

        double t_a = std::sqrt(t_a_sq);

        print("t_a", t_a);

        double t_min_a = t_0-t_a;
        double t_max_a = t_0+t_a;

        print("t_min_a", t_min_a);
        print("t_max_a", t_max_a);

        // compute t_min, t_max
        double t_min = std::max({t_min_C, t_min_tau, t_min_a});
        double t_max = std::min({t_max_C, t_max_tau, t_max_a});

        print("t_min", t_min);
        print("t_max", t_max);

        // Negative-size interval
        if (t_min > t_max) {
            throw std::runtime_error("Step C (left branch): Intersection of the intervals for t_S is empty.");
        // Zero-sized interval
        } else if (is_zero(t_min-t_max)) {
            abc = as_bs_cs(t_min, Sxx, Syy, Szz, t11, t22, t33, t12, t23, t13);
            if ( ! (is_negative(abc[1]) or is_positive(abc[2])) ) {
                t_S = t_min;
            } else {
                throw std::runtime_error("Step C (left branch): t_min=t_max gives non-positive semidefinite V.");
            }
        // Positive-size interval
        } else {
            // Initialise loop variables
            bool found_solution = false;
            double step = (t_max-t_min)/1.e5;
            double target = 1.e+9;
            double target_test;

            // Find the closest valid t to t_0
            for (double t_test=t_min; (t_test<=t_max); t_test+=step) {
                // How far is this from the target
                target_test = std::abs(t_0-t_test);
                // Check if closer than current best valid solution
                if (target_test < target) {
                    // Extract multipliers for the quadratic expansion
                    abc = as_bs_cs(t_test, Sxx, Syy, Szz, t11, t22, t33, t12, t23, t13);
                    // Check linear component positive and quadratic component negative
                    if ( ! (is_negative(abc[1]) or is_positive(abc[2])) ) {
                        // Store the new values
                        target = target_test;
                        t_S = t_test;
                        // Note that a valid solution has been found
                        found_solution = true;
                    }
                }
            }
            // No valid solution found in interval
            if (! found_solution) {
                throw std::runtime_error("Step C (left branch): Interval (t_min,t_max) has no t giving positive semidefinite V.");
            }
        }
    // One or more diagonal components of L is zero
    } else {

        std::cout << std::endl << "--------------------------" << std::endl;
        std::cout << "TAKING THE RIGHT BRANCH" << std::endl;
        std::cout << "--------------------------" << std::endl << std::endl;

        // whether a solution was found for the cauchy
        bool found_solution(false);

        for (int ii=0; ii<3; ii++) {
            // Generate matrix indices for diagonals (sym and non-sym)
            int s_i = (ii+0)%3;
            int s_j = (ii+1)%3;
            int s_k = (ii+2)%3;
            int m_i = (4*ii+0)%12;
            int m_j = (4*ii+4)%12;
            int m_k = (4*ii+8)%12;

            if ( is_zero(L_L[s_i]) ) {
                double t_test = S_L[m_i];
                // Check cauchy conditions
                double cp1 = (S_L[m_j]-t_test)*(S_L[m_j]-t_test) - T_CL[s_j]*L_L[s_j];
                double cp2 = (S_L[m_k]-t_test)*(S_L[m_k]-t_test) - T_CL[s_k]*L_L[s_k];
                if ( is_positive(cp1) or is_positive(cp2) ) {
                    throw std::runtime_error("Step C (right branch): Cauchy condition failed (23).");
                }
                // Check standard conditions
                abc = as_bs_cs(t_test, Sxx, Syy, Szz, t11, t22, t33, t12, t23, t13);
                if ( is_positive(abc[0]) or is_negative(abc[1]) or is_negative(abc[2]) ) {
                    throw std::runtime_error("Step C (right branch): Conditions 33-35 failed.");
                }
                // Checks passed -- does it agree with previous solutions?
                if ( found_solution and (! is_zero(t_S-t_test)) ) {
                    throw std::runtime_error("Step C (right branch): different solutions do not agree (should not happen?)");
                } else {
                    // Store the found solution
                    t_S = t_test;
                    // Note that a valid solution has been found
                    found_solution = true;
                }
            }
        }
        if ( ! found_solution ) {
            throw std::runtime_error("Step C (right branch): No valid solutions found.");
        }
    }

    // Check the final solution
    abc = as_bs_cs(t_S, Sxx, Syy, Szz, t11, t22, t33, t12, t23, t13);
    if ( is_negative(abc[1]) or is_positive(abc[2]) ) {
        throw std::runtime_error("Step C (left/right branch): Final t_S does not give positive semidefinite V.");
    }

    std::cout << "--------------------------" << std::endl;
    print("Final t_S", t_S);

    // Calculate the trace-corrected form of S_L
    S_C = S_L - scitbx::mat3<double>(t_S);

    std::cout << "--------------------------" << std::endl;
    print("S[L]", S_L);
    print("S_C[L]", S_C);

    double S_Cxx = xx(S_C);
    double S_Cyy = yy(S_C);
    double S_Czz = zz(S_C);

    double sx{0.0}, sy{0.0}, sz{0.0};

    if ( is_zero(Lxx) and (! is_zero(S_Cxx)) ) {
        throw std::runtime_error("Step C: incompatible L_L and S_C matrices.");
    } else {
        sx = S_Cxx/Lxx;
    }

    if ( is_zero(Lyy) and (! is_zero(S_Cyy)) ) {
        throw std::runtime_error("Step C: incompatible L_L and S_C matrices.");
    } else {
        sy = S_Cyy/Lyy;
    }

    if ( is_zero(Lzz) and (! is_zero(S_Czz)) ) {
        throw std::runtime_error("Step C: incompatible L_L and S_C matrices.");
    } else {
        sz = S_Czz/Lzz;
    }

    std::cout << "--------------------------" << std::endl;
    std::cout <<  "Screw parameters:" << std::endl;
    print("sx", sx);
    print("sy", sy);
    print("sz", sz);

    // Store s-amplitudes
    s_amplitudes = scitbx::vec3<double>(sx,sy,sz);

    // Calculate the contribution to the translation matrix
    C_L_t_S = scitbx::sym_mat3<double>(sx*S_Cxx, sy*S_Cyy, sz*S_Czz, 0.0, 0.0, 0.0);

    std::cout << "--------------------------" << std::endl;
    print("C[L](t=t_S)", C_L_t_S);

    // Calculate the vibration matrix
    V_L = T_CL - C_L_t_S;

    std::cout << "--------------------------" << std::endl;
    print("V[L]", V_L);

    if (! is_pd(V_L)) {
        throw std::runtime_error("Step C: Matrix V[L] is not positive semidefinite.");
    }

}

void decompose_tls_matrices::stepD() {

    std::cout << std::endl;
    std::cout << " ######################################" << std::endl;
    std::cout << " <-------        Step D        ------->" << std::endl;
    std::cout << " ######################################" << std::endl << std::endl;

    // Diagonalise V[L] matrix
    Eig v_eig(V_L);
    eigenvectors v_eig_vecs = v_eig.vectors();
    eigenvalues  v_eig_vals = v_eig.values();

    // Store the eigenvectors of V
    v_amplitudes = scitbx::vec3<double>(v_eig_vals[0], v_eig_vals[1], v_eig_vals[2]);
    v1_L = scitbx::vec3<double>(&v_eig_vecs[0]);
    v2_L = scitbx::vec3<double>(&v_eig_vecs[3]);
    v3_L = scitbx::vec3<double>(&v_eig_vecs[6]);
    // OVERWRITE - Choose v1 to be cross product of v2 and v3 (guaranteed righthanded system)
    v1_L = v2_L.cross(v3_L);

    // Print Eigenvalues of V
    std::cout << "--------------------------" << std::endl;
    std::cout << "Eigenvalues of V[L]: " << std::endl;
    print("1", v_amplitudes[0]);
    print("2", v_amplitudes[1]);
    print("3", v_amplitudes[2]);
    std::cout << "Eigenvectors of V[L]: " << std::endl;
    print("1", v1_L);
    print("2", v2_L);
    print("3", v3_L);

    // T values are sqrt of v
    double t1 = std::sqrt(v_amplitudes[0]);
    double t2 = std::sqrt(v_amplitudes[1]);
    double t3 = std::sqrt(v_amplitudes[2]);

    std::cout << "RMS vibrational components: " << std::endl;
    print("1", t1);
    print("2", t2);
    print("3", t3);

    // Extract the rotation matrix from the eigenvectors
    //scitbx::mat3<double> R_LV = static_cast <scitbx::mat3<double>> (v_eig.vectors().begin());
    scitbx::mat3<double> R_LV = scitbx::mat3<double>(v1_L[0], v1_L[1], v1_L[2], v2_L[0], v2_L[1], v2_L[2], v3_L[0], v3_L[1], v3_L[2]);
    scitbx::mat3<double> R_LV_t = R_LV.transpose();

    // Calculate V in own basis
    scitbx::mat3<double> V_V_ = R_LV * V_L * R_LV_t;
    V_V = scitbx::sym_mat3<double>((V_V_ + V_V_.transpose()) / 2.0);

    std::cout << "--------------------------" << std::endl;
    print("V[V]", V_V);

    // Return vibration axes to original basis
    v1_M = R_ML_t * v1_L;
    v2_M = R_ML_t * v2_L;
    v3_M = R_ML_t * v3_L;

    std::cout << "--------------------------" << std::endl;
    std::cout << "Vibrational axes (M-basis)" << std::endl;
    print("1", v1_M);
    print("2", v2_M);
    print("3", v3_M);

    // Calculate full transformation
    R_MV = R_LV * R_ML;

    std::cout << "--------------------------" << std::endl;
    print("R[M->L]", R_ML);
    print("R[L->V]", R_LV);
    print("R[M->V]", R_MV);

}

void decompose_tls_matrices::run() {

    try
    {
        stepA();
        stepB();
        stepC();
        stepD();
    }
    catch (std::runtime_error e)
    {
        std::cout << "Exception Caught: " << e.what() << std::endl;
        error_ = e.what();
    }

    // print result
    std::cout << std::boolalpha << "TLS matrices are valid: " << is_valid() << std::endl;
    if  (! is_valid() ) {
        std::cout << "Result: " << error_ << std::endl;
    }
}

}}} // close mmtbx/tls/decompose

int main()
{
    double one(1.0), nul(0.0), sml(1e-9);

    // Simple test (no S)
//    scitbx::sym_mat3<double> T (+0.280, +0.395, +0.234, +0.063, +0.118, +0.033);
//    scitbx::sym_mat3<double> L (+0.106, +0.107, +0.696, -0.021, -0.014, -0.152);
//    scitbx::mat3<double> S (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

    // Test against python implementation
    scitbx::sym_mat3<double> T (0.0959882450545,0.156509777229,0.0886002086034,-5.19222597403e-05,-0.0172511466857,-0.00869540261855);
    scitbx::sym_mat3<double> L (0.0123001438881,0.0461188308667,0.0218102205791,0.00435192109608,0.00517157269537,-0.0155699937912);
    scitbx::mat3<double> S (-0.00725668016697,-0.0201548281778,0.00431134160979,0.036329322896,-0.00252099816856,-0.00771964293756,-0.00968895922456,0.0170374138375,0.00977767833552);

    // Make matrices larger so can be seen after scaling to radians
    double mult(100.0);
    T *= mult;
    L *= mult;
    S *= mult;

    mmtbx::tls::decompose::decompose_tls_matrices dtm(T, L, S);

    return 0;
}
