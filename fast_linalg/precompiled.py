import platform
import libtbx.utils
import zipfile
import sys

openblas_version = '0.2.15-p1'
release_url_root = 'https://github.com/luc-j-bourhis/OpenBLAS/releases/download'

message = """
*** Configuration of module fast_linalg on this platform requires a
pre-compiled distribution of OpenBLAS but none is available. ***
Are you in the process of installing OpenBLAS with fast_linalg.build_openblas 
(y/n)? """

error_message = """
Then either the network is down; or no pre-compiled OpenBLAS is available for
your platform. In the latter case, your only option is to do:
libtbx.configure --exclude fast_linalg
and then to drop a message on the cctbx mailing list. Alternatively you can
try to install OpenBLAS in the cctbx build directory yourself with
fast_linalg.build_openblas.
"""

warning_message = """I assume you know what you are doing then and 
complete the configuration: libtbx.env.has_module("fast_linalg") will
return true although the fast_linalg module will be non-fonctional until
you install OpenBLAS with fast_linalg.build_openblas.
"""

def install_precompiled_openblas(platform_moniker):
  arch, _ = platform.architecture()
  if arch not in ('32bit', '64bit'):
    raise libtbx.utils.Sorry(
      'Module fast_linalg does not support %s architecture on this platform' %
      arch)
  url = '%s/v%s/openblas-%s-%s.zip' % (
    release_url_root, openblas_version, platform_moniker, arch)
  filename = 'openblas-%s.zip' % openblas_version
  try:
    libtbx.utils.retrieve_unless_exists(url, filename)
  except IOError:
    ans = None
    print message,
    sys.stdout.flush()
    while ans not in ('y', 'n'):
      ans = sys.stdin.readline().strip()
    print "\n"
    if ans == 'n':
      raise libtbx.utils.Sorry(error_message)
    else:
      print warning_message
      return
  archive = zipfile.ZipFile(filename)
  archive.extractall(abs(libtbx.env.build_path))
