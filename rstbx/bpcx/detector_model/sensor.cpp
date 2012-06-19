#include <scitbx/error.h>
#include "sensor.h"

rstbx::detector_model::sensor::sensor(
            const scitbx::vec3<double>& _origin,
            const scitbx::vec3<double>& _dir1,
            const scitbx::vec3<double>& _dir2,
            const scitbx::vec2<double>& _lim1,
            const scitbx::vec2<double>& _lim2):
    lim1(_lim1),
    lim2(_lim2),
    normal(),
    d(),
    D(),
    d_is_invertible(false)
{
    set_frame(_origin, _dir1, _dir2);
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

    // test that the input directions are not collinear
    // FIXME this assert apparently ineffective. Optimised out?
    assert(fabs(_dir1.angle(_dir2, deg = true) % 180.0) > 1.0);

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

    d.set_column(0, dir1);
    d.set_column(1, dir2);
    d.set_column(2, origin);
    try {
        D = d.inverse();
        d_is_invertible = true;
    }
    catch (scitbx::error const&) {
        D = scitbx::mat3<double>(.0,.0,.0,.0,.0,.0,.0,.0,.0);
        d_is_invertible = false;
    }
}
