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
   and diagonal matrix \f$\Delta\f$, the coordinates \f$y=\Delta^{-1/2} R^T\f$
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

    bool ndp;

  public:
    typedef scitbx::af::ref<GLdouble, scitbx::af::f_grid<2> >
            mat_ref_type;
    typedef scitbx::af::const_ref<GLdouble, scitbx::af::f_grid<2> >
            mat_const_ref_type;
    typedef scitbx::af::const_ref<GLdouble> vec_const_ref_type;

    /// Construct the change of frame from the frame where the ellipsoid
    /// is represented by the given metrics at the given centre
    ellipsoid_to_sphere_transform(scitbx::vec3<GLdouble> const &centre,
                                  scitbx::sym_mat3<GLdouble> const &metrics)
      : ndp(false)
    {
      mat_ref_type trs(m, 4, 4);
      scitbx::matrix::eigensystem::real_symmetric<GLdouble> es(metrics);
      mat_const_ref_type eigenvec(es.vectors().begin(), 3, 3);
      vec_const_ref_type eigenval(es.values().begin(), 3);
      for(int j=0; j<3; ++j) {
        if (eigenval[j] <= 0) {
          ndp = true;
          return;
        }
        GLdouble sqrt_eigenval = std::sqrt(eigenval[j]);
        for(int i=0; i<3; ++i) {
          trs(i,j) = sqrt_eigenval*eigenvec(i,j);
        }
      }
      for (int i=0; i<3; ++i) trs(i, 3) = centre[i];
      for (int j=0; j<3; ++j) trs(3, j) = 0;
      trs(3, 3) = 1;
    };

    /// Whether the metric was non-positive definite.
    /** If this is so, then the other members can't be relied upon */
    bool non_positive_definite() const { return ndp; }

    /// The rotation-scaling part of the change of basis
    scitbx::mat3<GLdouble> linear_part() const {
      GLTBX_ASSERT(!non_positive_definite());
      mat_const_ref_type r(m, 4, 4);
      return scitbx::mat3<GLdouble>(r(0,0), r(0,1), r(0,2),
                                    r(1,0), r(1,1), r(1,2),
                                    r(2,0), r(2,1), r(2,2));
    }

    /// The translation part of the change of basis
    scitbx::vec3<GLdouble> translation_part() const {
      GLTBX_ASSERT(!non_positive_definite());
      mat_const_ref_type r(m, 4, 4);
      return scitbx::vec3<GLdouble>(&r(0, 3));
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
