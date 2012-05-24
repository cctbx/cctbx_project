#include "sensor.h"

rstbx::detector_model::sensor::sensor(
            const scitbx::vec3<double>& _origin,
            const scitbx::vec3<double>& _dir1,
            const scitbx::vec3<double>& _dir2,
            const scitbx::vec2<double>& _lim1,
            const scitbx::vec2<double>& _lim2):
    origin(_origin),
    dir1(_dir1),
    dir2(_dir2),
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

void rstbx::detector_model::sensor::update()
{
    normal = (dir1.cross(dir2)).normalize();
    distance = origin * normal;

    D.set_column(0, dir1);
    D.set_column(1, dir2);
    D.set_column(2, origin);
    D = D.inverse();
}
