#ifndef GLTBX_QUADRICS_H
#define GLTBX_QUADRICS_H

#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/f_grid.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/matrix/eigensystem.h>

// these includes last to avoid Visual C++ 7.1, 8.0 failures
#include <gltbx/util.h>
#include <gltbx/error.h>

namespace gltbx { namespace quadrics {

namespace af = scitbx::af;

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
class proto_cylinder : public prototype
{
  public:
    typedef scitbx::vec3<GLdouble> vec3_t;

    proto_cylinder(GLdouble top_to_base_radius_ratio,
                   GLint slices, GLint stacks=1,
                   GLenum draw_style=GLU_FILL,
                   GLenum orientation=GLU_OUTSIDE,
                   GLenum normals=GLU_SMOOTH)
      : top_to_base_radius_ratio_(top_to_base_radius_ratio),
        slices(slices), stacks(stacks),
        prototype(draw_style, orientation, normals)
    {
      draw_prototype();
    }

    proto_cylinder(GLint slices,
                   GLenum draw_style=GLU_FILL,
                   GLenum orientation=GLU_OUTSIDE,
                   GLenum normals=GLU_SMOOTH)
      : top_to_base_radius_ratio_(1),
        slices(slices), stacks(1),
        prototype(draw_style, orientation, normals)
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
      GLdouble s = v.length(), c=e_z*u;
      GLdouble theta = scitbx::rad_as_deg(std::atan2(s,c));
      glRotated(theta, v[0], v[1], v[2]);
      glScaled(base_radius, base_radius, height);
      glCallList(index);
    }

    /// Draw the prototype as per specified by the arguments passed to the
    /// constructor.
    void draw_prototype() {
      GLTBX_SCOPE_PUSH_MATRIX;
      glLoadIdentity();
      glNewList(index, GL_COMPILE);
      gluCylinder(q, 1., top_to_base_radius_ratio_, 1., slices, stacks);
      glEndList();
    }

  private:
    GLdouble top_to_base_radius_ratio_;
    GLint slices, stacks;
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
    /// is represented by the given metrics at the given centre.
    /** If the metrics is non-positive-definite, each non-positive eigenvalue
        is changed to a tiny positive value. As a result, the displayed
        shape will be an ellipsoid flattened in the direction of the
        eigenvectors corresponding to those non-positive eigenvalues.
     */
    ellipsoid_to_sphere_transform(scitbx::vec3<GLdouble> const &centre,
                                  scitbx::sym_mat3<GLdouble> const &metrics)
      : npd(false)
    {
      scitbx::matrix::eigensystem::real_symmetric<GLdouble> es(metrics);
      af::const_ref<GLdouble, af::c_grid<2> > e = es.vectors().const_ref();
      scitbx::vec3<GLdouble> e0(e[0], e[1], e[2]),
                             e1(e[3], e[4], e[5]),
                             e2=e0.cross(e1);
      af::ref<GLdouble> eigenval = es.values().ref();
      GLdouble const tiny = +0.0005;
      for (int i=0; i<3; ++i) {
        if (eigenval[i] <= 0) {
          eigenval[i] = tiny;
          npd = true;
        }
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
    bool non_positive_definite() const { return npd; }

    /// The rotation-scaling part of the change of basis
    scitbx::mat3<GLdouble> linear_part() const {
      return scitbx::mat3<GLdouble>(m[0], m[4], m[ 8],
                                    m[1], m[5], m[ 9],
                                    m[2], m[6], m[10]);
    }

    /// The translation part of the change of basis
    scitbx::vec3<GLdouble> translation_part() const {
      return scitbx::vec3<GLdouble>(&m[12]);
    }

    GLdouble const *matrix() const {
      return m;
    }
};


class ellipsoid_principal_sections_texture
{
public:

  ellipsoid_principal_sections_texture(GLdouble darkening, int n_s, int n_t)
  {
    GLTBX_ASSERT(0 <= darkening && darkening <= 1)(darkening);
    GLubyte lumin = darkening*255;
    typedef af::f_grid<2> grid_t;
    af::versa<GLubyte, grid_t> texture_image(grid_t(n_s, n_t), 255);
    af::ref<GLubyte, grid_t> img = texture_image.ref();
    for (int s=0; s<n_s; ++s) {
      img(s, n_t/2) = lumin;
    }
    for (int t=0; t<n_t; ++t) {
      img(0, t) = img(n_s/4, t)
      = img(n_s/2, t) = img(3*n_s/4, t) = lumin;
    }

    glPushAttrib(GL_CLIENT_PIXEL_STORE_BIT);
    glGenTextures(1, &texture_name);
    glBindTexture(GL_TEXTURE_2D, texture_name);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_COMPRESSED_LUMINANCE, n_s, n_t, 0,
                 GL_LUMINANCE, GL_UNSIGNED_BYTE, img.begin());
    glBindTexture(GL_TEXTURE_2D, 0);
    glPopAttrib();
  }

  void bind() {
    glBindTexture(GL_TEXTURE_2D, texture_name);
  }

  void unbind() {
    glBindTexture(GL_TEXTURE_2D, 0);
  }

private:
  GLuint texture_name;
};



/// A prototype of an ellispoid stored in a display list
class proto_ellipsoid : public prototype
{
  public:
    typedef scitbx::vec3<GLdouble> vec3_t;
    typedef scitbx::sym_mat3<GLdouble> metrics_t;

    proto_ellipsoid(GLint slices, GLint stacks,
                    GLenum draw_style=GLU_FILL,
                    GLenum orientation=GLU_OUTSIDE,
                    GLenum normals=GLU_SMOOTH)
      : slices(slices), stacks(stacks),
        prototype(draw_style, orientation, normals)
    {
      gluQuadricTexture(q, GL_TRUE);
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
      glRotatef(90.f, 0.f, 1.f, 0.f); // bring z along the eigenvector
                                      // corresponding to the highest eigenvalue
      gluSphere(q, 1., slices, stacks);
      glEndList();
    }

  private:
    GLint slices, stacks;
};


}} // gltbx::quadrics

#endif // GUARD
