#pragma once

#include <scitbx/array_family/tiny.h>
#include <scitbx/vec3.h>

namespace scitbx { namespace math {

/// Tetrahedron
/**
 * T is the type of scalars used for the coordinates of the vertices.
 */
template<typename T>
class tetrahedron
{
public:
  /// Type of scalars
  typedef T scalar_t;

  /// Type of a vertex
  typedef vec3<scalar_t> vertex_t;

  /// Type of the sequence of vertices
  typedef af::tiny<vertex_t, 4> vertices_t;

  /// Type of each gradient wrt a vertex
  typedef vec3<scalar_t> gradient_t;

  /// Type of the sequence of gradients
  typedef af::tiny<gradient_t, 4> gradients_t;

  /// Construct a tetrahedron with the given vertices
  tetrahedron(vertices_t const &v) : v(v) {
    determinant = 1./6*(v[1]-v[0])*(v[2]-v[0]).cross(v[3]-v[0]);
  }

  /// Vertices
  vertices_t const &vertices() const { return v; }

  /// Volume
  scalar_t volume() { return std::abs(determinant); }

  /// Gradients of volume
  /**
   *  The i-th component of the result is the gradient of volume()
   *  with respect to the i-th vertex, for i=0,1,2,3
   *
   *  C.f. Mathematica notebook tetrahedron_volume.nb in the same directory
   *  as this file for a check of the formula.
   *
   * Those formulae are one possibility among the set of possible formulae
   * we wrote about in [1], eqn. (107). Unfortunately, that equation is actually
   * incorrect as there is a missing minus sign depending on the indices.
   * We should issue a correction.
   *
   * [1] Luc J. Bourhis, Oleg V. Dolomanov, Richard J. Gildea,
   * Judith A. K. Howard, and Horst Puschmann,
   * The anatomy of a comprehensive constrained, restrained refinement program
   * for the modern computing environment -- Olex2 dissected,
   * Acta Crystallographica Section A 71 (2015), 59--75
   */
  gradients_t gradients() const {
    scalar_t f = determinant > 0 ? +1 : -1;
    f *= 1./6;
    return gradients_t(f*(v[3] - v[1]).cross(v[2] - v[1]),
                       f*(v[2] - v[0]).cross(v[3] - v[0]),
                       f*(v[3] - v[0]).cross(v[1] - v[0]),
                       f*(v[1] - v[0]).cross(v[2] - v[0]));
  }

private:
  /// Vertices
  vertices_t v;

  scalar_t determinant;
};

}}
