#ifndef SCITBX_MATH_ICOSAHEDRON_H
#define SCITBX_MATH_ICOSAHEDRON_H

#include <scitbx/array_family/shared.h>
#include <math.h>

namespace scitbx { namespace math {
  template <typename FloatType = double>
    class icosahedron
    {
      public:
      int level;
      af::shared<scitbx::vec3<FloatType> > sites;
      private:
      af::shared<scitbx::vec3<FloatType> > v;
      FloatType distsq(scitbx::vec3<FloatType> s1, scitbx::vec3<FloatType> s2) {
        return (s1-s2).length();
      }

      scitbx::vec3<FloatType> midpt(scitbx::vec3<FloatType> s1, scitbx::vec3<FloatType> s2) {
        return (s1+s2)/2.0;
      }

      void make_icosahedron() {
        scitbx::vec3<FloatType> temp;
        static const FloatType mag_val = atan((1.0+sqrt(5.0))/2.0);
        static const FloatType s = sin(mag_val);
        static const FloatType c = cos(mag_val);
        FloatType a = s;
        FloatType b = c;

        for (int i = 0; i < 2; i++) {
          a = -a;
          for (int j = 0; j < 2; j++) {
            b = -b;

            temp[0] = 0.0;
            temp[1] = a;
            temp[2] = b;
            v.push_back(temp);

            temp[0] = b;
            temp[1] = 0.0;
            temp[2] = a;
            v.push_back(temp);

            temp[0] = a;
            temp[1] = b;
            temp[2] = 0.0;
            v.push_back(temp);
          }
        }
      }

      /*-
       * Recursively subdivides triangle ABC levl times, storing the result in
       * the probe array, PROBE->uco.
       *
       * STRATEGY: if not the bottom level, find the midpoints of each of the
       * sides of this sub_triangle and recurse, calling each of the four sub
       * -sub_triangles in turn.  If it's the bottom level, find the center, and
       * add that point to the output probe object.
       -*/
      void sub_triangle(scitbx::vec3<FloatType> a, scitbx::vec3<FloatType> b, scitbx::vec3<FloatType> c, int lvl) {
        scitbx::vec3<FloatType> d, e, f;

        if (lvl > 0) {
          /* find new midpoints */
          d = midpt(a, b);
          e = midpt(b, c);
          f = midpt(a, c);

          sub_triangle(a, d, f, lvl - 1);
          sub_triangle(d, b, e, lvl - 1);
          sub_triangle(d, e, f, lvl - 1);
          sub_triangle(f, e, c, lvl - 1);
        } else {
          /* normalize the vectors */
          d = a.normalize();
          e = b.normalize();
          f = c.normalize();
          /* find center */
          af::shared<scitbx::vec3<FloatType> > tmp;
          tmp.push_back(d);
          tmp.push_back(e);
          tmp.push_back(f);
          scitbx::vec3<FloatType> g;
          for(int i = 0; i < 3; i++) {
            g[i] = (d[i]+e[i]+f[i])/3.0;
          }
          sites.push_back(g.normalize());
        }
      }

      public:
      icosahedron() { }

      /*-
        builds a probe sphere with uniformly distributed points, by recursively
        sub-dividing an icosahedrom level times.
        -*/
      icosahedron(int level_) {
        int i, j, k;

        SCITBX_ASSERT(level_>=1);
        level=level_;
        sites.clear();

        make_icosahedron();
        for (i = 0; i < 10; i++) {
          for (j = i + 1; j < 11; j++) {
            if (distsq(v[i], v[j]) < 1.21) {
              for (k = j + 1; k < 12; k++) {
                if ((distsq(v[i], v[k]) < 1.21) &&
                    (distsq(v[j], v[k]) < 1.21))
                  sub_triangle(v[i], v[j], v[k],level);
              } /* k */
            } /* endif */
          } /* j */
        } /* i */
      }
    };
}}

#endif // SCITBX_MATH_ICOSAHEDRON_H
