#include "sensor.h"

rstbx::detector_model::sensor::sensor(
            const scitbx::vec3<double>& _origin,
            const scitbx::vec3<double>& _dir1,
            const scitbx::vec3<double>& _dir2,
            const scitbx::vec2<double>& _lim1,
            const scitbx::vec2<double>& _lim2):
    origin(_origin),
    dir1(_dir1.normalize()),
    dir2(_dir2.normalize()),
    lim1(_lim1),
    lim2(_lim2),
    normal(),
    D()
{
    update();
}

// other getters are inline
scitbx::mat3<double> rstbx::detector_model::sensor::get_d() const
{
    return D.inverse();
}

// setters
void rstbx::detector_model::sensor::set_origin(
                 const scitbx::vec3<double>& _origin)
{
    origin = _origin;
    update();
}

void rstbx::detector_model::sensor::set_frame(
                 const scitbx::vec3<double>& _origin,
                 const scitbx::vec3<double>& _dir1,
                 const scitbx::vec3<double>& _dir2)
{

  /* add some tests that the input directions _dir1, _dir2 are not colinear */

  /* assert(fabs(_dir1.angle(_dir2, deg = true) % 180.0) > 1.0); */

    origin = _origin;
    dir1 = _dir1.normalize();
    dir2 = _dir2.normalize();
    update();
}

void rstbx::detector_model::sensor::update()
{
    // ensure dir1, dir2 are orthonormal
    normal = (dir1.cross(dir2)).normalize();
    dir2 = normal.cross(dir1);

    distance = origin * normal;

    D.set_column(0, dir1);
    D.set_column(1, dir2);
    D.set_column(2, origin);
    D = D.inverse();
}
