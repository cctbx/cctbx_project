from __future__ import absolute_import, division, print_function

import logging
import pkg_resources


class ProfileModelFactory(object):
    """
    A factory to create a profile model
    """

    @staticmethod
    def from_dict(obj):
        """
        Given a dictionary, convert to a profile model
        """
        if obj is None:
            return None
        for entry_point in pkg_resources.iter_entry_points("dxtbx.profile_model"):
            if entry_point.name == obj["__id__"]:
                return entry_point.load().from_dict(obj)
        logging.getLogger("dxtbx.model.profile").warn(
            "No profile class %s registered" % obj["__id__"]
        )
        print(
            "dxtbx.model.profile: WARNING: No profile class %s registered"
            % obj["__id__"]
        )
        return None
