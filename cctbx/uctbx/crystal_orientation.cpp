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

void
cctbx::crystal_orientation::change_basis(
  cctbx::sgtbx::change_of_basis_op const& cb_op){
  Astar_ = Astar_ * (cb_op.c().r().as_double());
}

void
cctbx::crystal_orientation::change_basis(oc_mat3 const& rot){
  Astar_ = Astar_ * rot;
}

cctbx::oc_mat3
cctbx::crystal_orientation::best_similarity_transformation(
  crystal_orientation const& other, int unimodular_generator_range) const{

  scitbx::mat3<double> orientation_similarity(1); //initially the identity
  double minimum_orientation_bases_msd = direct_mean_square_difference(other);

  scitbx::math::unimodular_generator<double>
      unimodular_generator(unimodular_generator_range);
  while (!unimodular_generator.at_end()) {
      oc_mat3 c_inv_r = unimodular_generator.next();
      crystal_orientation ref_copy = other;
      ref_copy.change_basis(c_inv_r.inverse());
      double R = direct_mean_square_difference(ref_copy);
      if (R < minimum_orientation_bases_msd){
        minimum_orientation_bases_msd = R;
        orientation_similarity = c_inv_r;
      }
  }
  return orientation_similarity;
}

} // namespace cctbx
