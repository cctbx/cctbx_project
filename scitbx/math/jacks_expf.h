#ifndef SCITBX_MATH_JACKS_EXPF_H
#define SCITBX_MATH_JACKS_EXPF_H

#include <stdexcept>

namespace scitbx { namespace math {

//! http://forums.devshed.com/software-design-43/fast-exponential-function-algorithm-60039.html
/*! The only problem with the original code was the missing test for h != 0.
    Code further modified by applying obvious optimizations.
 */
inline
float
jacks_expf(
  float x)
{
  static const float binary[] = {
    1.0f,
    1.258925411794167f,
    1.584893192461113f,
    1.995262314968879f,
    2.511886431509580f,
    3.162277660168379f,
    3.981071705534972f,
    5.011872336272723f,
    6.309573444801932f,
    7.943282347242815f};
  static const float pow10tab[] = {
    1e00f,1e01f,1e02f,1e03f,1e04f,1e05f,1e06f,1e07f,1e08f,1e09f,
    1e10f,1e11f,1e12f,1e13f,1e14f,1e15f,1e16f,1e17f,1e18f,1e19f,
    1e20f,1e21f,1e22f,1e23f,1e24f,1e25f,1e26f,1e27f,1e28f,1e29f,
    1e30f,1e31f,1e32f,1e33f,1e34f,1e35f,1e36f,1e37f};
  if (   x > -2.980232061133847e-8f
      && x <  5.960463766996327e-8f) {
    return 1.0f;
  }
  float w = x * 0.434294481903251f;
  if (w < 0) w = -w;
  int g = (w < 38 ? int(w) : 38);
  if (g > 37) {
    if (x < 0) return 0;
    throw std::runtime_error("jacks_expf(): function argument out of range.");
  }
  int z = int((w-g)*10);
  float h = w - (float(g) + float(z)/10);
  float result = pow10tab[g] * binary[z];
  if (h != 0) {
    result *= (
      (10.423067565678043f / (
        5.211533782839021f - h - (
          9.430584850580696f / (
            1.886116970116139f / h + h)))) - 1);
  }
  if (x < 0) return 1 / result;
  return result;
}

}} // namespace scitbx::math

#endif // GUARD
