#include <cctbx/uctbx/fast_minimum_reduction.h>
#include <cctbx/sgtbx/change_of_basis_op.h>
#include <scitbx/math/unimodular_generator.h>
#include <scitbx/array_family/tiny_algebra.h>

namespace cctbx { namespace uctbx {

  namespace {

    void throw_corrupt_metrical_matrix()
    {
      throw error("Corrupt metrical matrix.");
    }

    double
    dot_g(uc_vec3 const& u, uc_sym_mat3 const& g, uc_vec3 const& v)
    {
      return u * (g * v);
    }

    uc_vec3
    cross_g(double sqrt_det_g, uc_sym_mat3 const& g,
            uc_vec3 const& r, uc_vec3 const& s)
    {
      return sqrt_det_g * (g * r).cross(g * s);
    }

    double acos_deg(double x) { return scitbx::rad_as_deg(std::acos(x)); }

    af::double6
    parameters_from_metrical_matrix(const double* metrical_matrix)
    {
      af::double6 params;
      for(std::size_t i=0;i<3;i++) {
        if (metrical_matrix[i] <= 0.) throw_corrupt_metrical_matrix();
        params[i] = std::sqrt(metrical_matrix[i]);
      }
      params[3] = acos_deg(metrical_matrix[5] / params[1] / params[2]);
      params[4] = acos_deg(metrical_matrix[4] / params[2] / params[0]);
      params[5] = acos_deg(metrical_matrix[3] / params[0] / params[1]);
      return params;
    }

    uc_sym_mat3
    construct_metrical_matrix(
      af::double6 const& params, uc_vec3 const& cos_ang)
    {
      return uc_sym_mat3(
       params[0] * params[0],
       params[1] * params[1],
       params[2] * params[2],
       params[0] * params[1] * cos_ang[2],
       params[0] * params[2] * cos_ang[1],
       params[1] * params[2] * cos_ang[0]);
    }

  } // namespace <anonymous>

  void unit_cell::init_volume()
  {
    /* V = a * b * c * sqrt(1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2
                              + 2 * cos(alpha) * cos(beta) * cos(gamma))
     */
    double d = 1.;
    for(std::size_t i=0;i<3;i++) d -= cos_ang_[i] * cos_ang_[i];
    d += 2. * cos_ang_[0] * cos_ang_[1] * cos_ang_[2];
    if (d < 0.) throw error("Square of unit cell volume is negative.");
        volume_ = params_[0] * params_[1] * params_[2] * std::sqrt(d);
    if (volume_ <= 0.) throw error("Unit cell volume is zero or negative.");
  }

  void unit_cell::init_reciprocal()
  {
    // Transformation Lattice Constants -> Reciprocal Lattice Constants
    // after Kleber, W., 17. Aufl., Verlag Technik GmbH Berlin 1990, P.352
    static const char*
      error_msg = "Error computing reciprocal unit cell angles.";
    for(std::size_t i=0;i<3;i++) r_params_[i] = params_[(i + 1) % 3]
                                              * params_[(i + 2) % 3]
                                              * sin_ang_[i] / volume_;
    for(std::size_t i=0;i<3;i++) {
      double denom = sin_ang_[(i + 1) % 3] * sin_ang_[(i + 2) % 3];
      if (denom == 0) throw error (error_msg);
      r_cos_ang_[i] = (  cos_ang_[(i + 1) % 3]
                       * cos_ang_[(i + 2) % 3]
                       - cos_ang_[i])
                    / denom;
    }
    for(std::size_t i=0;i<3;i++) {
      if (r_cos_ang_[i] < -1 || r_cos_ang_[i] > 1) {
        throw error(error_msg);
      }
      double a_rad = std::acos(r_cos_ang_[i]);
      r_params_[i+3] = scitbx::rad_as_deg(a_rad);
      r_sin_ang_[i] = std::sin(a_rad);
      r_cos_ang_[i] = std::cos(a_rad);
    }
  }

  void unit_cell::init_orth_and_frac_matrices()
  {
    // Crystallographic Basis: D = {a,b,c}
    // Cartesian Basis:        C = {i,j,k}
    //
    // PDB convention:
    //   i || a
    //   j is in (a,b) plane
    //   k = i x j

    double s1rca2 = std::sqrt(1. - r_cos_ang_[0] * r_cos_ang_[0]);
    if (s1rca2 == 0.) {
      throw error(
       "Reciprocal unit cell alpha angle is zero or extremely close to zero.");
    }

    // fractional to cartesian
    orth_[0] =  params_[0];
    orth_[1] =  cos_ang_[2] * params_[1];
    orth_[2] =  cos_ang_[1] * params_[2];
    orth_[3] =  0.;
    orth_[4] =  sin_ang_[2] * params_[1];
    orth_[5] = -sin_ang_[1] * r_cos_ang_[0] * params_[2];
    orth_[6] =  0.;
    orth_[7] =  0.;
    orth_[8] =  sin_ang_[1] * params_[2] * s1rca2;

    // cartesian to fractional
    frac_[0] =  1. / params_[0];
    frac_[1] = -cos_ang_[2] / (sin_ang_[2] * params_[0]);
    frac_[2] = -(  cos_ang_[2] * sin_ang_[1] * r_cos_ang_[0]
                 + cos_ang_[1] * sin_ang_[2])
             / (sin_ang_[1] * s1rca2 * sin_ang_[2] * params_[0]);
    frac_[3] =  0.;
    frac_[4] =  1. / (sin_ang_[2] * params_[1]);
    frac_[5] =  r_cos_ang_[0] / (s1rca2 * sin_ang_[2] * params_[1]);
    frac_[6] =  0.;
    frac_[7] =  0.;
    frac_[8] =  1. / (sin_ang_[1] * s1rca2 * params_[2]);
  }

  void unit_cell::init_metrical_matrices()
  {
    metr_mx_ = construct_metrical_matrix(params_, cos_ang_);
    r_metr_mx_ = construct_metrical_matrix(r_params_, r_cos_ang_);

    using scitbx::constants::pi_180;
    af::versa<double, af::c_grid<2> > &f = d_metrical_matrix_d_params_;
    f.resize(af::c_grid<2>(6,6), 0.);
    f(0,0) = 2*params_[0];
    f(1,1) = 2*params_[1];
    f(2,2) = 2*params_[2];
    f(3,0) = params_[1]*cos_ang_[2];
    f(3,1) = params_[0]*cos_ang_[2];
    f(3,5) = -params_[0]*params_[1]*sin_ang_[2]*pi_180;
    f(4,0) = params_[2]*cos_ang_[1];
    f(4,2) = params_[0]*cos_ang_[1];
    f(4,4) = -params_[0]*params_[2]*sin_ang_[1]*pi_180;
    f(5,1) = params_[2]*cos_ang_[0];
    f(5,2) = params_[1]*cos_ang_[0];
    f(5,3) = -params_[1]*params_[2]*sin_ang_[0]*pi_180;
  }

  void unit_cell::init_tensor_rank_2_orth_and_frac_linear_maps() {
    uc_mat3 const &o = orthogonalization_matrix();
    af::double6   &f = u_star_to_u_iso_linear_form_;
    f[0] = o[0]*o[0];
    f[1] = o[1]*o[1] + o[4]*o[4];
    f[2] = o[2]*o[2] + o[5]*o[5] + o[8]*o[8];
    f[3] = 2*o[0]*o[1];
    f[4] = 2*o[0]*o[2];
    f[5] = 2*(o[1]*o[2] + o[4]*o[5]);
    f *= 1./3;

    u_star_to_u_cart_linear_map_.resize(
      af::c_grid<2>(6,6), static_cast<double>(0));
    af::versa<double, af::c_grid<2> > &L = u_star_to_u_cart_linear_map_;
    L(0,0) = o[0]*o[0];
    L(0,1) = o[1]*o[1];
    L(0,2) = o[2]*o[2];
    L(0,3) = 2*o[0]*o[1];
    L(0,4) = 2*o[0]*o[2];
    L(0,5) = 2*o[1]*o[2];
    L(1,1) = o[4]*o[4];
    L(1,2) = o[5]*o[5];
    L(1,5) = 2*o[4]*o[5];
    L(2,2) = o[8]*o[8];
    L(3,1) = o[1]*o[4];
    L(3,2) = o[2]*o[5];
    L(3,3) = o[0]*o[4];
    L(3,4) = o[0]*o[5];
    L(3,5) = o[2]*o[4] + o[1]*o[5];
    L(4,2) = o[2]*o[8];
    L(4,4) = o[0]*o[8];
    L(4,5) = o[1]*o[8];
    L(5,2) = o[5]*o[8];
    L(5,5) = o[4]*o[8];

    af::double6  &c = u_star_to_u_cif_linear_map_;
    c[0] = 1./(r_params_[0] * r_params_[0]);
    c[1] = 1./(r_params_[1] * r_params_[1]);
    c[2] = 1./(r_params_[2] * r_params_[2]);
    c[3] = 1./(r_params_[0] * r_params_[1]);
    c[4] = 1./(r_params_[0] * r_params_[2]);
    c[5] = 1./(r_params_[1] * r_params_[2]);
  }

  void unit_cell::initialize()
  {
    std::size_t i;
    for(i=0;i<6;i++) {
      if (params_[i] <= 0.) {
        throw error("Unit cell parameter is zero or negative.");
      }
    }
    for(i=3;i<6;i++) {
      double a_deg = params_[i];
      if (a_deg >= 180.) {
        throw error(
          "Unit cell angle is greater than or equal to 180 degrees.");
      }
      double a_rad = scitbx::deg_as_rad(a_deg);
      cos_ang_[i-3] = std::cos(a_rad);
      sin_ang_[i-3] = std::sin(a_rad);
      if (sin_ang_[i-3] == 0.) {
        throw error("Unit cell angle is zero or or extremely close to zero.");
      }
    }
    init_volume();
    init_reciprocal();
    init_metrical_matrices();
    init_orth_and_frac_matrices();
    init_tensor_rank_2_orth_and_frac_linear_maps();
    longest_vector_sq_ = -1.;
    shortest_vector_sq_ = -1.;
  }

  unit_cell::unit_cell(af::small<double, 6> const& parameters)
  : params_(1,1,1,90,90,90)
  {
    std::copy(parameters.begin(), parameters.end(), params_.begin());
    initialize();
  }

  unit_cell::unit_cell(af::double6 const& parameters)
  : params_(parameters)
  {
    initialize();
  }

  unit_cell::unit_cell(uc_sym_mat3 const& metrical_matrix)
  : params_(parameters_from_metrical_matrix(metrical_matrix.begin()))
  {
    try {
      initialize();
    }
    catch (error const&) {
      throw_corrupt_metrical_matrix();
    }
  }

  unit_cell::unit_cell(uc_mat3 const& orthogonalization_matrix)
  : params_(parameters_from_metrical_matrix(
      orthogonalization_matrix.self_transpose_times_self().begin()))
  {
    try {
      initialize();
    }
    catch (error const&) {
      throw error("Corrupt orthogonalization matrix.");
    }
  }

  // used by reciprocal()
  unit_cell::unit_cell(
    af::double6 const& params,
    af::double3 const& sin_ang,
    af::double3 const& cos_ang,
    double volume,
    uc_sym_mat3 const& metr_mx,
    af::double6 const& r_params,
    af::double3 const& r_sin_ang,
    af::double3 const& r_cos_ang,
    uc_sym_mat3 const& r_metr_mx)
  :
    params_(params),
    sin_ang_(sin_ang),
    cos_ang_(cos_ang),
    volume_(volume),
    metr_mx_(metr_mx),
    r_params_(r_params),
    r_sin_ang_(r_sin_ang),
    r_cos_ang_(r_cos_ang),
    r_metr_mx_(r_metr_mx),
    longest_vector_sq_(-1.),
    shortest_vector_sq_(-1.)
  {
    init_orth_and_frac_matrices();
  }

  unit_cell
  unit_cell::reciprocal() const
  {
    return unit_cell(
      r_params_,
      r_sin_ang_,
      r_cos_ang_,
      1. / volume_,
      r_metr_mx_,
      params_,
      sin_ang_,
      cos_ang_,
      metr_mx_);
  }

  double
  unit_cell::longest_vector_sq() const
  {
    if (longest_vector_sq_ < 0.) {
      longest_vector_sq_ = 0.;
      int corner[3];
      for (corner[0] = 0; corner[0] <= 1; corner[0]++)
      for (corner[1] = 0; corner[1] <= 1; corner[1]++)
      for (corner[2] = 0; corner[2] <= 1; corner[2]++) {
        fractional<> xf;
        for(std::size_t i=0;i<3;i++) xf[i] = corner[i];
        scitbx::math::update_max(longest_vector_sq_, length_sq(xf));
      }
    }
    return longest_vector_sq_;
  }

  double
  unit_cell::shortest_vector_sq() const
  {
    if (shortest_vector_sq_ < 0.) {
      af::double6 gruber_params = fast_minimum_reduction<>(*this)
        .as_gruber_matrix();
      shortest_vector_sq_ = gruber_params[0];
      for(std::size_t i=1;i<3;i++) {
        scitbx::math::update_min(shortest_vector_sq_, gruber_params[i]);
      }
    }
    return shortest_vector_sq_;
  }

  bool
  unit_cell::is_degenerate(double min_min_length_over_max_length,
                           double min_volume_over_min_length)
  {
    if (volume_ == 0) return true;
    double min_length = std::min(std::min(params_[0], params_[1]), params_[2]);
    if (volume_ < min_length * min_volume_over_min_length) return true;
    double max_length = std::max(std::max(params_[0], params_[1]), params_[2]);
    if (min_length < max_length * min_min_length_over_max_length) return true;
    return false;
  }

  bool
  unit_cell::is_similar_to(unit_cell const& other,
                           double relative_length_tolerance,
                           double absolute_angle_tolerance) const
  {
    using scitbx::fn::absolute;
    const double* l1 = params_.begin();
    const double* l2 = other.params_.begin();
    for(std::size_t i=0;i<3;i++) {
      if (absolute(std::min(l1[i], l2[i]) / std::max(l1[i], l2[i]) - 1)
          > relative_length_tolerance) {
        return false;
      }
    }
    const double* a1 = l1 + 3;
    const double* a2 = l2 + 3;
    for(std::size_t i=0;i<3;i++) {
      if (absolute(a1[i] - a2[i]) > absolute_angle_tolerance) {
        return false;
      }
    }
    return true;
  }

  af::shared<scitbx::mat3<int> >
  unit_cell::similarity_transformations(
    unit_cell const& other,
    double relative_length_tolerance,
    double absolute_angle_tolerance,
    int unimodular_generator_range) const
  {
    af::shared<scitbx::mat3<int> > result;
    scitbx::mat3<int> identity_matrix(1);
    if (is_similar_to(
          other, relative_length_tolerance, absolute_angle_tolerance)) {
      result.push_back(identity_matrix);
    }
    scitbx::math::unimodular_generator<int>
      unimodular_generator(unimodular_generator_range);
    while (!unimodular_generator.at_end()) {
      scitbx::mat3<int> c_inv_r = unimodular_generator.next();
      unit_cell other_cb = other.change_basis(uc_mat3(c_inv_r));
      if (is_similar_to(
            other_cb, relative_length_tolerance, absolute_angle_tolerance)
          && c_inv_r != identity_matrix) {
        result.push_back(c_inv_r);
      }
    }
    return result;
  }

  unit_cell
  unit_cell::change_basis(uc_mat3 const& c_inv_r, double r_den) const
  {
    if (r_den == 0) {
      return unit_cell(metr_mx_.tensor_transpose_transform(c_inv_r));
    }
    return unit_cell(metr_mx_.tensor_transpose_transform(c_inv_r/r_den));
  }

  unit_cell
  unit_cell::change_basis(sgtbx::rot_mx const& c_inv_r) const
  {
    return change_basis(c_inv_r.as_double(), 0);
  }

  unit_cell
  unit_cell::change_basis(sgtbx::change_of_basis_op const& cb_op) const
  {
    return change_basis(cb_op.c_inv().r().as_double(), 0);
  }

  miller::index<>
  unit_cell::max_miller_indices(double d_min, double tolerance) const
  {
    CCTBX_ASSERT(d_min > 0);
    CCTBX_ASSERT(tolerance >= 0);
    miller::index<> max_h;
    for(std::size_t i=0;i<3;i++) {
      uc_vec3 u(0,0,0);
      uc_vec3 v(0,0,0);
      u[(i + 1) % 3] = 1.;
      v[(i + 2) % 3] = 1.;
      // length of uxv is not used => sqrt(det(metr_mx)) is simply set to 1
      uc_vec3 uxv = cross_g(1., r_metr_mx_, u, v);
      double uxv2 = dot_g(uxv, r_metr_mx_, uxv);
      CCTBX_ASSERT(uxv2 != 0);
      double uxv_abs = std::sqrt(uxv2);
      CCTBX_ASSERT(uxv_abs != 0);
      max_h[i] = static_cast<int>(uxv[i] / uxv_abs / d_min + tolerance);
    }
    return max_h;
  }

  int
  unit_cell::compare_orthorhombic(
    const unit_cell& other) const
  {
    af::double6 const& lhs = params_;
    af::double6 const& rhs = other.params_;
    for(std::size_t i=0;i<3;i++) {
      if (lhs[i] < rhs[i]) return -1;
      if (lhs[i] > rhs[i]) return  1;
    }
    return 0;
  }

  int
  unit_cell::compare_monoclinic(
    const unit_cell& other,
    unsigned unique_axis,
    double angular_tolerance) const
  {
    CCTBX_ASSERT(unique_axis < 3);
    af::double6 const& lhs = params_;
    af::double6 const& rhs = other.params_;
    double lhs_ang = lhs[unique_axis+3];
    double rhs_ang = rhs[unique_axis+3];
    using scitbx::fn::absolute;
    if (absolute(lhs_ang - rhs_ang) < angular_tolerance) {
      return compare_orthorhombic(other);
    }
    double lhs_ang_d90 = absolute(lhs_ang - 90);
    double rhs_ang_d90 = absolute(rhs_ang - 90);
    if (absolute(lhs_ang_d90 - rhs_ang_d90) > angular_tolerance) {
      if (lhs_ang_d90 < rhs_ang_d90) return -1;
      if (lhs_ang_d90 > rhs_ang_d90) return  1;
    }
    else {
      if (lhs_ang > 90 && rhs_ang < 90) return -1;
      if (lhs_ang < 90 && rhs_ang > 90) return  1;
    }
    if (lhs_ang > rhs_ang) return -1;
    if (lhs_ang < rhs_ang) return  1;
    return 0;
  }

  sgtbx::change_of_basis_op const&
  unit_cell::change_of_basis_op_for_best_monoclinic_beta() const
  {
    static const sgtbx::change_of_basis_op cb_op_identity;
    static const sgtbx::change_of_basis_op cb_op_acbc("a+c,b,c");
    unit_cell alt = change_basis(cb_op_acbc);
    double beta = params_[4];
    double beta_alt = alt.params_[4];
    if (beta_alt >= 90 && beta_alt < beta) return cb_op_acbc;
    return cb_op_identity;
  }

}} // namespace cctbx::uctbx
