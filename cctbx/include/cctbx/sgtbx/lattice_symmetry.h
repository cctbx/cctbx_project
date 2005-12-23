#ifndef CCTBX_SGTBX_LATTICE_SYMMETRY_H
#define CCTBX_SGTBX_LATTICE_SYMMETRY_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/math/utils.h>
#include <scitbx/constants.h>

namespace cctbx { namespace sgtbx {

//! Determination of lattice symmetry (Bravais type).
namespace lattice_symmetry {

  //! Determines maximum Le Page (1982) delta of all two-fold axes.
  /*! See also: group()
   */
  double
  find_max_delta(
    uctbx::unit_cell const& reduced_cell,
    sgtbx::space_group const& space_group);

  struct reduced_cell_two_fold_info
  {
    sg_mat3 r;
    sg_vec3 u;
    sg_vec3 h;
  };

  //! Reformatted output of cctbx/cctbx/examples/reduced_cell_two_folds.py
  static const reduced_cell_two_fold_info reduced_cell_two_folds[] = {
    sg_mat3(-1,-1,-1,0,0,1,0,1,0),sg_vec3(-1,1,1),sg_vec3(0,1,1),
    sg_mat3(-1,-1,0,0,1,0,0,-1,-1),sg_vec3(1,-2,1),sg_vec3(0,1,0),
    sg_mat3(-1,-1,0,0,1,0,0,0,-1),sg_vec3(-1,2,0),sg_vec3(0,1,0),
    sg_mat3(-1,-1,0,0,1,0,0,1,-1),sg_vec3(-1,2,1),sg_vec3(0,1,0),
    sg_mat3(-1,-1,1,0,0,-1,0,-1,0),sg_vec3(1,-1,1),sg_vec3(0,-1,1),
    sg_mat3(-1,0,-1,0,-1,-1,0,0,1),sg_vec3(-1,-1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,-1,0,-1,0,0,0,1),sg_vec3(-1,0,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,-1,0,-1,1,0,0,1),sg_vec3(-1,1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,0,-1,0,-1,1,-1,0),sg_vec3(0,-1,1),sg_vec3(1,-1,1),
    sg_mat3(-1,0,0,-1,0,1,-1,1,0),sg_vec3(0,1,1),sg_vec3(-1,1,1),
    sg_mat3(-1,0,0,-1,1,-1,0,0,-1),sg_vec3(0,1,0),sg_vec3(1,-2,1),
    sg_mat3(-1,0,0,-1,1,0,0,0,-1),sg_vec3(0,1,0),sg_vec3(-1,2,0),
    sg_mat3(-1,0,0,-1,1,1,0,0,-1),sg_vec3(0,1,0),sg_vec3(-1,2,1),
    sg_mat3(-1,0,0,0,-1,-1,0,0,1),sg_vec3(0,-1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,0,0,-1,0,-1,-1,1),sg_vec3(0,0,1),sg_vec3(-1,-1,2),
    sg_mat3(-1,0,0,0,-1,0,-1,0,1),sg_vec3(0,0,1),sg_vec3(-1,0,2),
    sg_mat3(-1,0,0,0,-1,0,-1,1,1),sg_vec3(0,0,1),sg_vec3(-1,1,2),
    sg_mat3(-1,0,0,0,-1,0,0,-1,1),sg_vec3(0,0,1),sg_vec3(0,-1,2),
    sg_mat3(-1,0,0,0,-1,0,0,0,1),sg_vec3(0,0,1),sg_vec3(0,0,1),
    sg_mat3(-1,0,0,0,-1,0,0,1,1),sg_vec3(0,0,1),sg_vec3(0,1,2),
    sg_mat3(-1,0,0,0,-1,0,1,-1,1),sg_vec3(0,0,1),sg_vec3(1,-1,2),
    sg_mat3(-1,0,0,0,-1,0,1,0,1),sg_vec3(0,0,1),sg_vec3(1,0,2),
    sg_mat3(-1,0,0,0,-1,0,1,1,1),sg_vec3(0,0,1),sg_vec3(1,1,2),
    sg_mat3(-1,0,0,0,-1,1,0,0,1),sg_vec3(0,1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,0,0,0,-1,0,-1,0),sg_vec3(0,-1,1),sg_vec3(0,-1,1),
    sg_mat3(-1,0,0,0,0,1,0,1,0),sg_vec3(0,1,1),sg_vec3(0,1,1),
    sg_mat3(-1,0,0,0,1,-1,0,0,-1),sg_vec3(0,1,0),sg_vec3(0,-2,1),
    sg_mat3(-1,0,0,0,1,0,0,-1,-1),sg_vec3(0,-2,1),sg_vec3(0,1,0),
    sg_mat3(-1,0,0,0,1,0,0,0,-1),sg_vec3(0,1,0),sg_vec3(0,1,0),
    sg_mat3(-1,0,0,0,1,0,0,1,-1),sg_vec3(0,2,1),sg_vec3(0,1,0),
    sg_mat3(-1,0,0,0,1,1,0,0,-1),sg_vec3(0,1,0),sg_vec3(0,2,1),
    sg_mat3(-1,0,0,1,0,-1,-1,-1,0),sg_vec3(0,-1,1),sg_vec3(-1,-1,1),
    sg_mat3(-1,0,0,1,0,1,1,1,0),sg_vec3(0,1,1),sg_vec3(1,1,1),
    sg_mat3(-1,0,0,1,1,-1,0,0,-1),sg_vec3(0,1,0),sg_vec3(-1,-2,1),
    sg_mat3(-1,0,0,1,1,0,0,0,-1),sg_vec3(0,1,0),sg_vec3(1,2,0),
    sg_mat3(-1,0,0,1,1,1,0,0,-1),sg_vec3(0,1,0),sg_vec3(1,2,1),
    sg_mat3(-1,0,1,0,-1,-1,0,0,1),sg_vec3(1,-1,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,1,0,-1,0,0,0,1),sg_vec3(1,0,2),sg_vec3(0,0,1),
    sg_mat3(-1,0,1,0,-1,1,0,0,1),sg_vec3(1,1,2),sg_vec3(0,0,1),
    sg_mat3(-1,1,-1,0,0,-1,0,-1,0),sg_vec3(-1,-1,1),sg_vec3(0,-1,1),
    sg_mat3(-1,1,0,0,1,0,0,-1,-1),sg_vec3(-1,-2,1),sg_vec3(0,1,0),
    sg_mat3(-1,1,0,0,1,0,0,0,-1),sg_vec3(1,2,0),sg_vec3(0,1,0),
    sg_mat3(-1,1,0,0,1,0,0,1,-1),sg_vec3(1,2,1),sg_vec3(0,1,0),
    sg_mat3(-1,1,1,0,0,1,0,1,0),sg_vec3(1,1,1),sg_vec3(0,1,1),
    sg_mat3(0,-1,-1,-1,0,1,0,0,-1),sg_vec3(-1,1,0),sg_vec3(-1,1,1),
    sg_mat3(0,-1,-1,0,-1,0,-1,1,0),sg_vec3(-1,0,1),sg_vec3(-1,1,1),
    sg_mat3(0,-1,0,-1,0,0,-1,1,-1),sg_vec3(-1,1,1),sg_vec3(-1,1,0),
    sg_mat3(0,-1,0,-1,0,0,0,0,-1),sg_vec3(-1,1,0),sg_vec3(-1,1,0),
    sg_mat3(0,-1,0,-1,0,0,1,-1,-1),sg_vec3(1,-1,1),sg_vec3(-1,1,0),
    sg_mat3(0,-1,1,-1,0,-1,0,0,-1),sg_vec3(-1,1,0),sg_vec3(1,-1,1),
    sg_mat3(0,-1,1,0,-1,0,1,-1,0),sg_vec3(1,0,1),sg_vec3(1,-1,1),
    sg_mat3(0,0,-1,-1,-1,1,-1,0,0),sg_vec3(-1,1,1),sg_vec3(-1,0,1),
    sg_mat3(0,0,-1,0,-1,0,-1,0,0),sg_vec3(-1,0,1),sg_vec3(-1,0,1),
    sg_mat3(0,0,-1,1,-1,-1,-1,0,0),sg_vec3(-1,-1,1),sg_vec3(-1,0,1),
    sg_mat3(0,0,1,-1,-1,-1,1,0,0),sg_vec3(1,-1,1),sg_vec3(1,0,1),
    sg_mat3(0,0,1,0,-1,0,1,0,0),sg_vec3(1,0,1),sg_vec3(1,0,1),
    sg_mat3(0,0,1,1,-1,1,1,0,0),sg_vec3(1,1,1),sg_vec3(1,0,1),
    sg_mat3(0,1,-1,0,-1,0,-1,-1,0),sg_vec3(-1,0,1),sg_vec3(-1,-1,1),
    sg_mat3(0,1,-1,1,0,-1,0,0,-1),sg_vec3(1,1,0),sg_vec3(-1,-1,1),
    sg_mat3(0,1,0,1,0,0,-1,-1,-1),sg_vec3(-1,-1,1),sg_vec3(1,1,0),
    sg_mat3(0,1,0,1,0,0,0,0,-1),sg_vec3(1,1,0),sg_vec3(1,1,0),
    sg_mat3(0,1,0,1,0,0,1,1,-1),sg_vec3(1,1,1),sg_vec3(1,1,0),
    sg_mat3(0,1,1,0,-1,0,1,1,0),sg_vec3(1,0,1),sg_vec3(1,1,1),
    sg_mat3(0,1,1,1,0,1,0,0,-1),sg_vec3(1,1,0),sg_vec3(1,1,1),
    sg_mat3(1,-1,-1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,1,1),
    sg_mat3(1,-1,0,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,1,0),
    sg_mat3(1,-1,1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,-1,1),
    sg_mat3(1,0,-1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,0,1),
    sg_mat3(1,0,0,-1,-1,0,-1,0,-1),sg_vec3(-2,1,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,-1,-1,0,0,0,-1),sg_vec3(-2,1,0),sg_vec3(1,0,0),
    sg_mat3(1,0,0,-1,-1,0,1,0,-1),sg_vec3(2,-1,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,0,-1,0,-1,0,-1),sg_vec3(-2,0,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(1,0,0),
    sg_mat3(1,0,0,0,-1,0,1,0,-1),sg_vec3(2,0,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,1,-1,0,-1,0,-1),sg_vec3(-2,-1,1),sg_vec3(1,0,0),
    sg_mat3(1,0,0,1,-1,0,0,0,-1),sg_vec3(2,1,0),sg_vec3(1,0,0),
    sg_mat3(1,0,0,1,-1,0,1,0,-1),sg_vec3(2,1,1),sg_vec3(1,0,0),
    sg_mat3(1,0,1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,0,1),
    sg_mat3(1,1,-1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,-1,1),
    sg_mat3(1,1,0,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,1,0),
    sg_mat3(1,1,1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,1,1)
  };

  /*! References:
        Y. Le Page
        The derivation of the axes of the conventional unit cell from the
        dimensions of the Buerger-reduced cell
        J. Appl. Cryst. (1982). 15, 255-259

        Andrey A. Lebedev, Alexei A. Vagin & Garib N. Murshudov
        Acta Cryst. (2006). D62, 83-95.
        http://journals.iucr.org/d/issues/2006/01/00/ba5089/index.html
        Appendix A1. Algorithms used in the determination of twinning
          operators and their type of merohedry
  */
  space_group
  group(
    uctbx::unit_cell const& reduced_cell,
    double max_delta=3.,
    bool enforce_max_delta_for_generated_two_folds=true);

}}} // namespace cctbx::sgtbx::lattice_symmetry

#endif // CCTBX_SGTBX_LATTICE_SYMMETRY_H
