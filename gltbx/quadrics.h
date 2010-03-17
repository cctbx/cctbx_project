#ifndef GLTBX_QUADRICS_H
#define GLTBX_QUADRICS_H

#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/accessors/f_grid.h>
#include <scitbx/matrix/eigensystem.h>

// these includes last to avoid Visual C++ 7.1, 8.0 failures
#include <gltbx/util.h>
#include <gltbx/error.h>

namespace gltbx { namespace quadrics {

namespace af = scitbx::af;

template <class Heir>
class prototype
{
  public:
    prototype(GLenum draw_style=GLU_FILL,
              GLenum orientation=GLU_OUTSIDE,
              GLenum normals=GLU_SMOOTH)
      : q(gluNewQuadric())
    {
      gluQuadricDrawStyle(q, draw_style);
      gluQuadricOrientation(q, orientation);
      gluQuadricNormals(q, normals);
      index = glGenLists(1);
      if (index == 0) handle_error();
    }

    ~prototype() {
      gluDeleteQuadric(q);
      glDeleteLists(index, 1);
    }

  protected:
    GLUquadricObj *q;
    GLuint index;
};


/// A prototype of a cylinder stored in a diplay list
class proto_cylinder : public prototype<proto_cylinder>
{
  public:
    typedef scitbx::vec3<GLdouble> vec3_t;

    proto_cylinder(GLdouble top_to_base_radius_ratio,
                   GLint slices, GLint stacks,
                   GLenum draw_style=GLU_FILL,
                   GLenum orientation=GLU_OUTSIDE,
                   GLenum normals=GLU_SMOOTH)
      : top_to_base_radius_ratio_(top_to_base_radius_ratio),
        slices_(slices), stacks_(stacks),
        prototype<proto_cylinder>(draw_style, orientation, normals)
    {
      draw_prototype();
    }

    proto_cylinder(GLint slices,
                   GLenum draw_style=GLU_FILL,
                   GLenum orientation=GLU_OUTSIDE,
                   GLenum normals=GLU_SMOOTH)
      : top_to_base_radius_ratio_(1),
        slices_(slices), stacks_(1),
        prototype<proto_cylinder>(draw_style, orientation, normals)
    {
      draw_prototype();
    }

    /// Draw a cylinder from start to end with the given radius at start.
    void draw(vec3_t const &start, vec3_t const &end, GLdouble base_radius)
    {
      GLTBX_SCOPE_PUSH_MATRIX;
      vec3_t u = end - start;
      GLdouble height = u.length();
      if (height == 0) return;
      glTranslated(start[0], start[1], start[2]);
      u = u/height;
      vec3_t e_z (0,0,1);
      vec3_t v = e_z.cross(u);
      GLdouble l = v.length();
      if (l < 1e-8) {
        if (e_z*u < 0) height = -height;
      }
      else {
        GLdouble theta = scitbx::rad_as_deg(std::asin(v.length()));
        glRotated(theta, v[0], v[1], v[2]);
      }
      glScaled(base_radius, base_radius, height);
      glCallList(index);
    }

    /// Draw the prototype as per specified by the arguments passed to the
    /// constructor.
    void draw_prototype() {
      GLTBX_SCOPE_PUSH_MATRIX;
      glLoadIdentity();
      glNewList(index, GL_COMPILE);
      gluCylinder(q, 1., top_to_base_radius_ratio_, 1., slices_, stacks_);
      glEndList();
    }

  private:
    GLdouble top_to_base_radius_ratio_;
    GLint slices_, stacks_;
};


/// The change of coordinates frame from the frame where an ellipsoid is
/// represented by a given symmetric matrix to the frame where it is a sphere.
/* Those two frames shall be orthonormal.
   That change of frame is decomposed into:

    - first a translation of the coordinates frame to the centre of
      the ellipsoid,
    - then a rotation of the axes so that they end up along the principal axes
      of the ellipsoid,
    - finally a scaling of each axis unit length to match the 1/2 lengths of
      the ellipsoid;

   Mathematically, the ellispoid is the locus of points whose coordinates in the
   starting frame are \f$x\f$ and which satisfy
   \f$x^T U^{-1} x = cst\f$
   where \f$U\f$ is the metrics passed to the constructor.
   Since \f$U = R \Delta R^T\f$ for some orthogonal matrix \f$R\f$
   and diagonal matrix \f$\Delta\f$, the coordinates \f$y=\Delta^{-1/2} R^T x\f$
   satisfy \f$y^T y = cst\f$, i.e. this is a sphere in those coordinates.
   Thus the operator to transform the starting frame to the frame where
   coordinates are \f$y\f$, is \f$M = R \Delta^{1/2}\f$. Moreover, \f$R\f$ is
   the matrix whose columns are the eigenvectors of \f$U\f$ and \f$\Delta\f$
   has the corresponding eigenvalues. Hence we can use
   \c scitbx::matrix::eigensystem::real_symmetric to get those matrices.
*/
class ellipsoid_to_sphere_transform
{
  private:
    GLdouble m[16];

    bool npd;

  public:
    /// Do-nothing constructor
    ellipsoid_to_sphere_transform() {}

    /// Construct the change of frame from the frame where the ellipsoid
    /// is represented by the given metrics at the given centre
    ellipsoid_to_sphere_transform(scitbx::vec3<GLdouble> const &centre,
                                  scitbx::sym_mat3<GLdouble> const &metrics)
      : npd(false)
    {
      scitbx::matrix::eigensystem::real_symmetric<GLdouble> es(metrics);
      af::const_ref<GLdouble, af::c_grid<2> > e = es.vectors().const_ref();
      scitbx::vec3<GLdouble> e0(e[0], e[1], e[2]),
                             e1(e[3], e[4], e[5]),
                             e2=e0.cross(e1);
      af::const_ref<GLdouble> eigenval = es.values().const_ref();
      if (eigenval[0] <= 0 || eigenval[1] <= 0 || eigenval[2] <= 0) {
        npd = true;
        return;
      }
      e0 *= std::sqrt(eigenval[0]);
      e1 *= std::sqrt(eigenval[1]);
      e2 *= std::sqrt(eigenval[2]);

      m[0] = e0[0]; m[4] = e1[0]; m[ 8] = e2[0]; m[12] = centre[0];
      m[1] = e0[1]; m[5] = e1[1]; m[ 9] = e2[1]; m[13] = centre[1];
      m[2] = e0[2]; m[6] = e1[2]; m[10] = e2[2]; m[14] = centre[2];
      m[3] =   0  ; m[7] =   0  ; m[11] =   0  ; m[15] =     1    ;
    };

    /// Whether the metric was non-positive definite.
    /** If this is so, then the other members can't be relied upon */
    bool non_positive_definite() const { return npd; }

    /// The rotation-scaling part of the change of basis
    scitbx::mat3<GLdouble> linear_part() const {
      GLTBX_ASSERT(!non_positive_definite());
      return scitbx::mat3<GLdouble>(m[0], m[4], m[ 8],
                                    m[1], m[5], m[ 9],
                                    m[2], m[6], m[10]);
    }

    /// The translation part of the change of basis
    scitbx::vec3<GLdouble> translation_part() const {
      GLTBX_ASSERT(!non_positive_definite());
      return scitbx::vec3<GLdouble>(&m[12]);
    }

    GLdouble const *matrix() const {
      GLTBX_ASSERT(!non_positive_definite());
      return m;
    }
};


/// A prototype of an ellispoid stored in a display list
class proto_ellipsoid : public prototype<proto_ellipsoid>
{
  public:
    typedef scitbx::vec3<GLdouble> vec3_t;
    typedef scitbx::sym_mat3<GLdouble> metrics_t;

    proto_ellipsoid(GLint slices, GLint stacks,
                    GLenum draw_style=GLU_FILL,
                    GLenum orientation=GLU_OUTSIDE,
                    GLenum normals=GLU_SMOOTH)
      : slices_(slices), stacks_(stacks),
        prototype<proto_ellipsoid>(draw_style, orientation, normals)
    {
      draw_prototype();
    }

    /// Draw an ellipsoid at the given centre with the current metrics
    /** We use \c ellipsoid_to_sphere_transform to do the math

      We don't use glTranslate, glRotate and glScale to break the change
      of frame as it is  built by \c ellipsoid_to_sphere_transformso because
      glRotate would require us to compute the angle and vector of the rotation
      which would then be used to compute the rotation matrix on the graphic card:
      much more efficient to pass that matrix directly to OpenGL. But then
      we can as well fill the translation part and do the scaling part, hence
      doing it all in one glMultMatrix.
    */
    void draw(vec3_t const &centre, metrics_t const &metrics) {
      ellipsoid_to_sphere_transform m(centre, metrics);
      draw(m);
    }

    /// An overloaded version which opens the possibility of saving the
    /// expensive computation of the change of basis for each ellipsoid
    void draw(ellipsoid_to_sphere_transform const &m) {
      GLTBX_SCOPE_PUSH_MATRIX;
      glMultMatrixd(m.matrix());
      glCallList(index);
    }

    /// Draw the prototype as per specified by the arguments passed to the
    /// constructor.
    void draw_prototype() {
      GLTBX_SCOPE_PUSH_MATRIX;
      glLoadIdentity();
      glNewList(index, GL_COMPILE);
      gluSphere(q, 1., slices_, stacks_);
      glEndList();
    }

  private:
    GLint slices_, stacks_;
};


}} // gltbx::quadrics

#endif // GUARD
