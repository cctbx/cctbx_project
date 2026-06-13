'''
iotbx.file_io: fast file-type detection (get_file_type) and a parse-once read
facade (read_file) for DataManager datatypes. The modern entry point that
supersedes direct use of iotbx.file_reader.any_file for detection/content.
'''

from iotbx.file_io._maps import (  # noqa: F401
  any_file_type, data_manager_type, extension_to_datatype, CIF_SENTINEL)
from iotbx.file_io.detection import get_file_type  # noqa: F401
from iotbx.file_io.reader import read_file, FileIOResult  # noqa: F401
