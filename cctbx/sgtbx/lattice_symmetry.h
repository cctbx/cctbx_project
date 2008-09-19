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

  //! Reformatted output of cctbx/examples/reduced_cell_two_folds.py
  static const reduced_cell_two_fold_info reduced_cell_two_folds[] = {
    {sg_mat3(-1,-1,-1,0,0,1,0,1,0),sg_vec3(-1,1,1),sg_vec3(0,1,1)},
    {sg_mat3(-1,-1,0,0,1,0,0,-1,-1),sg_vec3(1,-2,1),sg_vec3(0,1,0)},
    {sg_mat3(-1,-1,0,0,1,0,0,0,-1),sg_vec3(-1,2,0),sg_vec3(0,1,0)},
    {sg_mat3(-1,-1,0,0,1,0,0,1,-1),sg_vec3(-1,2,1),sg_vec3(0,1,0)},
    {sg_mat3(-1,-1,1,0,0,-1,0,-1,0),sg_vec3(1,-1,1),sg_vec3(0,-1,1)},
    {sg_mat3(-1,0,-1,0,-1,-1,0,0,1),sg_vec3(-1,-1,2),sg_vec3(0,0,1)},
    {sg_mat3(-1,0,-1,0,-1,0,0,0,1),sg_vec3(-1,0,2),sg_vec3(0,0,1)},
    {sg_mat3(-1,0,-1,0,-1,1,0,0,1),sg_vec3(-1,1,2),sg_vec3(0,0,1)},
    {sg_mat3(-1,0,0,-1,0,-1,1,-1,0),sg_vec3(0,-1,1),sg_vec3(1,-1,1)},
    {sg_mat3(-1,0,0,-1,0,1,-1,1,0),sg_vec3(0,1,1),sg_vec3(-1,1,1)},
    {sg_mat3(-1,0,0,-1,1,-1,0,0,-1),sg_vec3(0,1,0),sg_vec3(1,-2,1)},
    {sg_mat3(-1,0,0,-1,1,0,0,0,-1),sg_vec3(0,1,0),sg_vec3(-1,2,0)},
    {sg_mat3(-1,0,0,-1,1,1,0,0,-1),sg_vec3(0,1,0),sg_vec3(-1,2,1)},
    {sg_mat3(-1,0,0,0,-1,-1,0,0,1),sg_vec3(0,-1,2),sg_vec3(0,0,1)},
    {sg_mat3(-1,0,0,0,-1,0,-1,-1,1),sg_vec3(0,0,1),sg_vec3(-1,-1,2)},
    {sg_mat3(-1,0,0,0,-1,0,-1,0,1),sg_vec3(0,0,1),sg_vec3(-1,0,2)},
    {sg_mat3(-1,0,0,0,-1,0,-1,1,1),sg_vec3(0,0,1),sg_vec3(-1,1,2)},
    {sg_mat3(-1,0,0,0,-1,0,0,-1,1),sg_vec3(0,0,1),sg_vec3(0,-1,2)},
    {sg_mat3(-1,0,0,0,-1,0,0,0,1),sg_vec3(0,0,1),sg_vec3(0,0,1)},
    {sg_mat3(-1,0,0,0,-1,0,0,1,1),sg_vec3(0,0,1),sg_vec3(0,1,2)},
    {sg_mat3(-1,0,0,0,-1,0,1,-1,1),sg_vec3(0,0,1),sg_vec3(1,-1,2)},
    {sg_mat3(-1,0,0,0,-1,0,1,0,1),sg_vec3(0,0,1),sg_vec3(1,0,2)},
    {sg_mat3(-1,0,0,0,-1,0,1,1,1),sg_vec3(0,0,1),sg_vec3(1,1,2)},
    {sg_mat3(-1,0,0,0,-1,1,0,0,1),sg_vec3(0,1,2),sg_vec3(0,0,1)},
    {sg_mat3(-1,0,0,0,0,-1,0,-1,0),sg_vec3(0,-1,1),sg_vec3(0,-1,1)},
    {sg_mat3(-1,0,0,0,0,1,0,1,0),sg_vec3(0,1,1),sg_vec3(0,1,1)},
    {sg_mat3(-1,0,0,0,1,-1,0,0,-1),sg_vec3(0,1,0),sg_vec3(0,-2,1)},
    {sg_mat3(-1,0,0,0,1,0,0,-1,-1),sg_vec3(0,-2,1),sg_vec3(0,1,0)},
    {sg_mat3(-1,0,0,0,1,0,0,0,-1),sg_vec3(0,1,0),sg_vec3(0,1,0)},
    {sg_mat3(-1,0,0,0,1,0,0,1,-1),sg_vec3(0,2,1),sg_vec3(0,1,0)},
    {sg_mat3(-1,0,0,0,1,1,0,0,-1),sg_vec3(0,1,0),sg_vec3(0,2,1)},
    {sg_mat3(-1,0,0,1,0,-1,-1,-1,0),sg_vec3(0,-1,1),sg_vec3(-1,-1,1)},
    {sg_mat3(-1,0,0,1,0,1,1,1,0),sg_vec3(0,1,1),sg_vec3(1,1,1)},
    {sg_mat3(-1,0,0,1,1,-1,0,0,-1),sg_vec3(0,1,0),sg_vec3(-1,-2,1)},
    {sg_mat3(-1,0,0,1,1,0,0,0,-1),sg_vec3(0,1,0),sg_vec3(1,2,0)},
    {sg_mat3(-1,0,0,1,1,1,0,0,-1),sg_vec3(0,1,0),sg_vec3(1,2,1)},
    {sg_mat3(-1,0,1,0,-1,-1,0,0,1),sg_vec3(1,-1,2),sg_vec3(0,0,1)},
    {sg_mat3(-1,0,1,0,-1,0,0,0,1),sg_vec3(1,0,2),sg_vec3(0,0,1)},
    {sg_mat3(-1,0,1,0,-1,1,0,0,1),sg_vec3(1,1,2),sg_vec3(0,0,1)},
    {sg_mat3(-1,1,-1,0,0,-1,0,-1,0),sg_vec3(-1,-1,1),sg_vec3(0,-1,1)},
    {sg_mat3(-1,1,0,0,1,0,0,-1,-1),sg_vec3(-1,-2,1),sg_vec3(0,1,0)},
    {sg_mat3(-1,1,0,0,1,0,0,0,-1),sg_vec3(1,2,0),sg_vec3(0,1,0)},
    {sg_mat3(-1,1,0,0,1,0,0,1,-1),sg_vec3(1,2,1),sg_vec3(0,1,0)},
    {sg_mat3(-1,1,1,0,0,1,0,1,0),sg_vec3(1,1,1),sg_vec3(0,1,1)},
    {sg_mat3(0,-1,-1,-1,0,1,0,0,-1),sg_vec3(-1,1,0),sg_vec3(-1,1,1)},
    {sg_mat3(0,-1,-1,0,-1,0,-1,1,0),sg_vec3(-1,0,1),sg_vec3(-1,1,1)},
    {sg_mat3(0,-1,0,-1,0,0,-1,1,-1),sg_vec3(-1,1,1),sg_vec3(-1,1,0)},
    {sg_mat3(0,-1,0,-1,0,0,0,0,-1),sg_vec3(-1,1,0),sg_vec3(-1,1,0)},
    {sg_mat3(0,-1,0,-1,0,0,1,-1,-1),sg_vec3(1,-1,1),sg_vec3(-1,1,0)},
    {sg_mat3(0,-1,1,-1,0,-1,0,0,-1),sg_vec3(-1,1,0),sg_vec3(1,-1,1)},
    {sg_mat3(0,-1,1,0,-1,0,1,-1,0),sg_vec3(1,0,1),sg_vec3(1,-1,1)},
    {sg_mat3(0,0,-1,-1,-1,1,-1,0,0),sg_vec3(-1,1,1),sg_vec3(-1,0,1)},
    {sg_mat3(0,0,-1,0,-1,0,-1,0,0),sg_vec3(-1,0,1),sg_vec3(-1,0,1)},
    {sg_mat3(0,0,-1,1,-1,-1,-1,0,0),sg_vec3(-1,-1,1),sg_vec3(-1,0,1)},
    {sg_mat3(0,0,1,-1,-1,-1,1,0,0),sg_vec3(1,-1,1),sg_vec3(1,0,1)},
    {sg_mat3(0,0,1,0,-1,0,1,0,0),sg_vec3(1,0,1),sg_vec3(1,0,1)},
    {sg_mat3(0,0,1,1,-1,1,1,0,0),sg_vec3(1,1,1),sg_vec3(1,0,1)},
    {sg_mat3(0,1,-1,0,-1,0,-1,-1,0),sg_vec3(-1,0,1),sg_vec3(-1,-1,1)},
    {sg_mat3(0,1,-1,1,0,-1,0,0,-1),sg_vec3(1,1,0),sg_vec3(-1,-1,1)},
    {sg_mat3(0,1,0,1,0,0,-1,-1,-1),sg_vec3(-1,-1,1),sg_vec3(1,1,0)},
    {sg_mat3(0,1,0,1,0,0,0,0,-1),sg_vec3(1,1,0),sg_vec3(1,1,0)},
    {sg_mat3(0,1,0,1,0,0,1,1,-1),sg_vec3(1,1,1),sg_vec3(1,1,0)},
    {sg_mat3(0,1,1,0,-1,0,1,1,0),sg_vec3(1,0,1),sg_vec3(1,1,1)},
    {sg_mat3(0,1,1,1,0,1,0,0,-1),sg_vec3(1,1,0),sg_vec3(1,1,1)},
    {sg_mat3(1,-1,-1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,1,1)},
    {sg_mat3(1,-1,0,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,1,0)},
    {sg_mat3(1,-1,1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,-1,1)},
    {sg_mat3(1,0,-1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,0,1)},
    {sg_mat3(1,0,0,-1,-1,0,-1,0,-1),sg_vec3(-2,1,1),sg_vec3(1,0,0)},
    {sg_mat3(1,0,0,-1,-1,0,0,0,-1),sg_vec3(-2,1,0),sg_vec3(1,0,0)},
    {sg_mat3(1,0,0,-1,-1,0,1,0,-1),sg_vec3(2,-1,1),sg_vec3(1,0,0)},
    {sg_mat3(1,0,0,0,-1,0,-1,0,-1),sg_vec3(-2,0,1),sg_vec3(1,0,0)},
    {sg_mat3(1,0,0,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(1,0,0)},
    {sg_mat3(1,0,0,0,-1,0,1,0,-1),sg_vec3(2,0,1),sg_vec3(1,0,0)},
    {sg_mat3(1,0,0,1,-1,0,-1,0,-1),sg_vec3(-2,-1,1),sg_vec3(1,0,0)},
    {sg_mat3(1,0,0,1,-1,0,0,0,-1),sg_vec3(2,1,0),sg_vec3(1,0,0)},
    {sg_mat3(1,0,0,1,-1,0,1,0,-1),sg_vec3(2,1,1),sg_vec3(1,0,0)},
    {sg_mat3(1,0,1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,0,1)},
    {sg_mat3(1,1,-1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(-2,-1,1)},
    {sg_mat3(1,1,0,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,1,0)},
    {sg_mat3(1,1,1,0,-1,0,0,0,-1),sg_vec3(1,0,0),sg_vec3(2,1,1)}
  };

  /*! References:
        http://www.ccp4.ac.uk/newsletters/newsletter44/articles/explore_metric_symmetry.html

        From section "2.1. Determination of the lattice symmetry"
        in the CCP4 newsletter 44:

          The determination of the lattice symmetry is based on
          ideas by Le Page (1982) and Lebedev et al. (2006). Given
          a reduced cell (e.g. Grosse-Kunstleve et al. 2004a), it is
          sufficient to search for two-fold axes to determine the full
          symmetry. Subjecting the two-folds to group multiplication
          produces the higher-order symmetry elements, if present.

          Le Page (1982) searches for the two-folds by computing
          angles between certain vectors in direct space and reciprocal
          space. This search is relatively expensive. Recently Lebedev
          et al. (2006) introduced the idea of simply enumerating all
          3x3 matrices with elements {-1,0,1} and determinant one. As
          an additional requirement group multiplication based on each
          matrix individually has to produce matrices exclusively with
          elements {-1,0,1}. There are only 480 matrices that conform
          to all requirements. Lebedev et al. (2006) argue that this
          set covers all possible symmetry operations for reduced
          cells. We were able to confirm this intuitive argument
          empirically via simple brute-force tests.

          Only 81 of the 480 selected matrices correspond to
          two-folds. These are easily detected by establishing
          which of the matrices produce the identity matrix when
          multiplied with themselves (and are not the identity matrix
          to start out with). To replace the expensive search for
          two-folds in the original Le Page (1982) algorithm, the
          81 two-fold matrices are tabulated along with the axis
          directions in direct space and reciprocal space. The axis
          direction in direct space is determined as described by
          Grosse-Kunstleve (1999). The axis direction in reciprocal
          space is determined with the same algorithm, but using
          the transpose of the matrix. The complete implementation
          of the algorithm for generating the table (essentially
          just six lines of Python code) can be found in the file
          cctbx/examples/reduced_cell_two_folds.py in the cctbx
          distributions.

          The search for two-folds computes the Le Page (1982) delta
          for each of the 81 tabulated pairs of axis directions. The
          corresponding symmetry matrix for each pair is immediately
          available from the table. In contrast, the original Le Page
          algorithm requires the evaluation of 2391 pairs of axis
          directions, and the computation of the symmetry matrices
          involves expensive trigonometric functions (sin, cos)
          and change-of-basis calculations.

          The matrices with a Le Page delta smaller than a given
          threshold are sorted, smallest delta first. Successive
          group multiplication as described in Grosse-Kunstleve
          et al. (2004b) and Sauter et al. (2006) yields
          the final highest lattice symmetry. The complete
          search algorithm is implemented in the file
          cctbx/sgtbx/lattice_symmetry.cpp.

        Grosse-Kunstleve, R.W. (1999). Acta Cryst. A55, 383-395.

        Grosse-Kunstleve, R.W., Sauter, N.K. & Adams, P.D. (2004a).
        Acta Cryst. A60, 1-6.

        Grosse-Kunstleve, R.W., Sauter, N.K. & Adams, P.D. (2004b).
        Newsletter of the IUCr Commission on Crystallographic Computing 3,
        22-31.

        Y. Le Page
        The derivation of the axes of the conventional unit cell from the
        dimensions of the Buerger-reduced cell
        J. Appl. Cryst. (1982). 15, 255-259

        Andrey A. Lebedev, Alexei A. Vagin & Garib N. Murshudov
        Acta Cryst. (2006). D62, 83-95.
        http://journals.iucr.org/d/issues/2006/01/00/ba5089/index.html
        Appendix A1. Algorithms used in the determination of twinning
          operators and their type of merohedry

        Sauter, N.K., Grosse-Kunstleve, R.W. & Adams, P.D. (2006).
        J. Appl. Cryst.  39, 158-168.
   */
  space_group
  group(
    uctbx::unit_cell const& reduced_cell,
    double max_delta=3.,
    bool enforce_max_delta_for_generated_two_folds=true);

}}} // namespace cctbx::sgtbx::lattice_symmetry

#endif // CCTBX_SGTBX_LATTICE_SYMMETRY_H
