from __future__ import absolute_import, division, print_function
from boost_adaptbx.deprecation_tools import deprecate_method
import warnings


def exercise():
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")

    # Choose a Boost Python class for which to deprecate a method. Here we
    # use rational.int.
    from boost_adaptbx.boost import rational

    original_value = rational.int().numerator()
    deprecate_method(rational.int, "numerator")
    new_value = rational.int().numerator()

  assert original_value == new_value
  assert len(w) == 1
  assert issubclass(w[-1].category, DeprecationWarning)
  assert "deprecated" in str(w[-1].message)


def run():
  exercise()
  print("OK")


if __name__ == "__main__":
  run()
