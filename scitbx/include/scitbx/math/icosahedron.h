#ifndef SCITBX_MATH_ICOSAHEDRON_H
#define SCITBX_MATH_ICOSAHEDRON_H

#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <cmath>

namespace scitbx { namespace math {

  //! Cartesian coordinates of points on a sphere.
  /*! This code was contributed by Erik McKee, Texas A&M University.
   */
  template <typename FloatType=double>
  class icosahedron
  {
    public:
      unsigned level;
      af::shared<scitbx::vec3<FloatType> > sites;

    private:
      static
      void
      make_icosahedron(scitbx::vec3<FloatType>* v)
      {
        static const FloatType phi_ratio = (1.0+std::sqrt(5.0))/2.0;
        static const FloatType b_ = 1.0 / (std::sqrt(1.0+phi_ratio*phi_ratio));
        static const FloatType a_ = phi_ratio * b_;
        FloatType a = a_;
        FloatType b = b_;
        for (unsigned i = 0; i < 2; i++) {
          a = -a;
          for (unsigned j = 0; j < 2; j++) {
            b = -b;
            (*v)[0] = 0.0;
            (*v)[1] = a;
            (*v)[2] = b;
            v++;
            (*v)[0] = b;
            (*v)[1] = 0.0;
            (*v)[2] = a;
            v++;
            (*v)[0] = a;
            (*v)[1] = b;
            (*v)[2] = 0.0;
            v++;
          }
        }
      }

      void
      sub_triangle(
        scitbx::vec3<FloatType> const& a,
        scitbx::vec3<FloatType> const& b,
        scitbx::vec3<FloatType> const& c,
        unsigned lvl)
      {
        if (lvl == 0) {
          sites.push_back((a+b+c).normalize());
        }
        else {
          /* compute midpoints */
          scitbx::vec3<FloatType> d = (a+b).normalize();
          scitbx::vec3<FloatType> e = (b+c).normalize();
          scitbx::vec3<FloatType> f = (a+c).normalize();
          sub_triangle(a, d, f, lvl-1);
          sub_triangle(d, b, e, lvl-1);
          sub_triangle(d, e, f, lvl-1);
          sub_triangle(f, e, c, lvl-1);
        }
      }

    public:
      //! Default constructor. Some data members are not initialized!
      icosahedron() {}

      /*! Builds points on a sphere.
       */
      /*! level=0 builds the vertices of the base icosahedron. The
          sites are uniformly distributed.

          level>0 recursively sub-divides the 20 triangles of the base
          icosahedron level times. The number of sites is therefore
          80 * 4**(level-1). The sites are approximately uniformly
          distributed.
       */
      icosahedron(unsigned level_)
      :
        level(level_)
      {
        if (level == 0) {
          sites.resize(12);
          make_icosahedron(sites.begin());
        }
        else {
          af::tiny<scitbx::vec3<FloatType>, 12> v;
          make_icosahedron(v.begin());
          std::size_t four_pow_level_minus_one = 1;
          for (unsigned i=0;i<level-1;i++) four_pow_level_minus_one *= 4;
          sites.reserve(80 * four_pow_level_minus_one);
          for (unsigned i = 0; i < 10; i++) {
            for (unsigned j = i + 1; j < 11; j++) {
              if ((v[i]-v[j]).length_sq() < 1.2) {
                for (unsigned k = j + 1; k < 12; k++) {
                  if (   (v[i]-v[k]).length_sq() < 1.2
                      && (v[j]-v[k]).length_sq() < 1.2)
                    sub_triangle(v[i], v[j], v[k], level);
                }
              }
            }
          }
          SCITBX_ASSERT(sites.size() == 80 * four_pow_level_minus_one);
        }
      }

      //! Maximum distance to the next neighbors.
      /*! For level=0 this is the distance to the five next neighbors
          of each icosahedron vertex.

          For level>0 this is the maximum distance of any site
          to the three next neighbors.

          Distances are tabulated for levels 0 through 7.
          An exception is thrown if level>7.
       */
      double
      next_neighbors_distance() const
      {
        static const af::tiny<double, 8> known_distances(
          1.05146222424,
          0.353098248494,
          0.185386249656,
          0.0947464326266,
          0.0476510500603,
          0.0238609877705,
          0.011934950279,
          0.00596803292972);
        if (level < known_distances.size()) return known_distances[level];
        throw std::runtime_error("next_neighbors_distance not known.");
      }
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_ICOSAHEDRON_H
