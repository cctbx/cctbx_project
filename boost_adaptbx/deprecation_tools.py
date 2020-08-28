from __future__ import absolute_import, division, print_function
import warnings


def deprecate_method(boost_object, method_name):
  original_method = getattr(boost_object, method_name)

  def deprecation_helper(*args, **kwargs):
    warnings.warn(
      "The method {method_name} is deprecated and will be removed shortly".format(
        method_name=method_name
      ),
      DeprecationWarning,
      stacklevel=2,
    )
    return original_method(*args, **kwargs)

  setattr(boost_object, method_name, deprecation_helper)
