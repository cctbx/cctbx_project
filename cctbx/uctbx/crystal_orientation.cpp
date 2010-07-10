#include <cctbx/crystal_orientation.h>
#include <scitbx/math/unimodular_generator.h>

namespace cctbx {

cctbx::crystal_orientation::crystal_orientation(
  oc_mat3 const& matrix, bool const& reciprocal_flag):
  Astar_(reciprocal_flag?matrix:matrix.inverse()){}

cctbx::uctbx::unit_cell
cctbx::crystal_orientation::unit_cell() const {
  oc_mat3 direct = direct_matrix();
  oc_vec3 p[3];
  for (int i=0; i<3; ++i){
     p[i]  = direct.get_row(i);
  }
  return cctbx::uctbx::unit_cell(scitbx::sym_mat3<double>(
    p[0]*p[0], p[1]*p[1], p[2]*p[2], p[0]*p[1], p[2]*p[0], p[1]*p[2]));
}

cctbx::uctbx::unit_cell
cctbx::crystal_orientation::unit_cell_inverse() const {
  return unit_cell().reciprocal();
}

cctbx::oc_mat3
cctbx::crystal_orientation::direct_matrix() const{
  return Astar_.inverse();
}

cctbx::oc_mat3
cctbx::crystal_orientation::reciprocal_matrix() const{
  return Astar_;
}

crystal_orientation
cctbx::crystal_orientation::change_basis(
  cctbx::sgtbx::change_of_basis_op const& cb_op) const {
  return crystal_orientation(
    Astar_ * (cb_op.c().r().as_double()), reciprocal );
}

crystal_orientation
cctbx::crystal_orientation::change_basis(oc_mat3 const& rot) const {
  return crystal_orientation( Astar_ * rot, reciprocal );
}

cctbx::oc_mat3
cctbx::crystal_orientation::best_similarity_transformation(
  crystal_orientation const& other,
  double const& fractional_length_tolerance,
  int unimodular_generator_range) const{

  scitbx::mat3<double> orientation_similarity(1); //initially the identity
  double minimum_orientation_bases_zsc = difference_Z_score(other);

  scitbx::math::unimodular_generator<int>
      unimodular_generator(unimodular_generator_range);
  while (!unimodular_generator.at_end()) {
      scitbx::mat3<int> m = unimodular_generator.next();
      oc_mat3 c_inv_r(m(0,0), m(0,1), m(0,2),
                      m(1,0), m(1,1), m(1,2),
                      m(2,0), m(2,1), m(2,2));
      crystal_orientation mod_copy = other.change_basis(c_inv_r.inverse());
      double R = difference_Z_score(mod_copy);
      if (R < minimum_orientation_bases_zsc){
        minimum_orientation_bases_zsc = R;
        orientation_similarity = c_inv_r;
      }
  }
  SCITBX_ASSERT(difference_Z_score(
    other.change_basis(orientation_similarity.inverse())) < fractional_length_tolerance);
  return orientation_similarity;
}

} // namespace cctbx
