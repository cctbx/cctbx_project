from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatPYunspecified import FormatPYunspecified


class FormatPYCXI41(FormatPYunspecified):

    """PREFERENCE FILE.
     Treats any Pickle-format file lacking a DETECTOR_FORMAT_VERSION key
     automatically as format CXI 4.1.
     Rename this file to FormatPYCXI41.py and put in dxtbx/format path.
  """

    def _start(self):

        FormatPYunspecified.start_helper(
            self, version_token="distl.detector_format_version=CXI 4.1"
        )


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatPYCXI41.understand(arg))
