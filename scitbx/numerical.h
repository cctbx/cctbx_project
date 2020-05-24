#pragma once
#include <scitbx/array_family/shared.h>
#include <limits>
#include <cmath>

namespace scitbx { namespace math { namespace numerical {
  /* numerical differentiation utilities
  http://en.wikipedia.org/wiki/Numerical_differentiation
  evaluator is supposed to provide the folloving interface:
  double calculate(size_t index) - calculates the function value when index
  variable varies
  */
  template <typename float_t>
  struct differential {
    typedef scitbx::af::shared<float_t> out_array_t;

    static float_t get_delta() {
      static const float_t delta = sqrt(std::numeric_limits<float_t>::epsilon());
      return delta;
    }

    // symmetric 2-point differential
    template <typename in_array_t, typename evaluator_t>
    static out_array_t diff_2(in_array_t &vars, evaluator_t const &e) {
      const float_t delta = get_delta();
      out_array_t rv(vars.size());
      for (size_t i = 0; i < vars.size(); i++) {
        vars[i] += delta;
        float_t vl = e.calculate(i);
        vars[i] -= 2*delta;
        float_t vr = e.calculate(i);
        rv[i] = (vl -vr) / delta;
        vars[i] += delta;
      }
      return rv;
    }

    // symmetric 4-point differential
    template <typename in_array_t, typename evaluator_t>
    static out_array_t diff_4(in_array_t &vars, evaluator_t const &e) {
      const float_t delta = get_delta();
      out_array_t rv(vars.size());
      for (size_t i = 0; i < vars.size(); i++) {
        vars[i] += 2 * delta;
        const float_t v1 = e.calculate(i);
        vars[i] -= delta;
        const float_t v2 = e.calculate(i);
        vars[i] -= 2 * delta;
        const float_t v3 = e.calculate(i);
        vars[i] -= delta;
        const float_t v4 = e.calculate(i);
        vars[i] += 2 * delta;
        rv[i] = (-v1 + 8 * v2 - 8 * v3 + v4) / (12 * delta);
      }
      return rv;
    }
};

}}}
