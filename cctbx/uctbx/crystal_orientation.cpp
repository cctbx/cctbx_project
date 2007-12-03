#include <cctbx/crystal_orientation.h>

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

} // namespace cctbx
