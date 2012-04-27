#include "sensor.h"
#include <iostream>

using namespace std;

//trivial main function for sensor
int main()
{
    scitbx::vec3<double> origin(0, 0, 100);
    scitbx::vec3<double> dir1(1, 0, 0);
    scitbx::vec3<double> dir2(0, 1, 0);
    scitbx::vec2<double> lim(0, 50);

    rstbx::detector_model::sensor s(origin, dir1, dir2, lim, lim);

    cout << s.get_distance() << endl;

    return 0;
}
