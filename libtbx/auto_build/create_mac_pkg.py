
"""
Create a graphical Mac installer for a complete installation of CCTBX or
Phenix (or any other derived software).
"""

from __future__ import absolute_import, division, print_function

import os
import os.path as op
import shutil
import sys
import time
from optparse import SUPPRESS_HELP, OptionParser

from .installer_utils import *

# XXX HACK
libtbx_path = op.abspath(op.dirname(op.dirname(__file__)))
if (not libtbx_path in sys.path):
  sys.path.append(libtbx_path)

def run(args, out=sys.stdout):
  if (sys.platform != "darwin"):
    print("ERROR: this program can only be run on Macintosh systems.")
    return False
  # XXX this prevents tar on OS X from including resource fork files, which
  # break the object relocation.  thanks to Francis Reyes for pointing this out.
  os.environ["COPYFILE_DISABLE"] = "true"
  datestamp = time.strftime("%Y-%m-%d", time.localtime())
  parser = OptionParser()
  parser.add_option("--tmp_dir", dest="tmp_dir", action="store",
    help="Temporary (staging) directory", default=os.getcwd())
  parser.add_option("--dist-dir", dest="dist_dir", action="store",
    help="Distribution directory", default="dist")
  parser.add_option("--package_name", dest="package_name", action="store",
    help="Package name", default="CCTBX")
  parser.add_option("--version", dest="version", action="store",
    help="Software version", default=datestamp)
  parser.add_option("--license_file", dest="license_file", action="store",
    help="License file", default=None)
  parser.add_option("--background", dest="background", action="store",
    help="Background image", default=None)
  parser.add_option("--organization", dest="organization", action="store",
    help="Organization (domain name)", default="org.phenix-online")
  parser.add_option("--machine_type", dest="machine_type", action="store",
    help="Machine type (OS/architecture)", default=machine_type())
  parser.add_option("--overwrite", dest="overwrite", action="store_true",
    help="Overwrite temporary files")
  parser.add_option("--no_compression", dest="no_compression",
    action="store_true", help=SUPPRESS_HELP) # legacy parameter
    # .pkg files are already compressed. Double-compressing is pointless.
  options, args = parser.parse_args(args)
  arch_type = os.uname()[-1]
  system_type = "64-bit Intel Macs"
  if (arch_type != "x86_64"):
    system_type = "Intel Macs"
  system_type += " (OS %s or later)" % get_os_version()
  if (len(args) != 1):
    print("Usage: create_mac_pkg [options] PROGRAM_DIR")
    return False
  program_dir = args[0]
  if (not program_dir.startswith("/")):
    print("ERROR: absolute path required")
    return False
  if (not op.isdir(program_dir)):
    print("ERROR: '%s' is not a directory" % program_dir)
    return False
  dest_name = op.dirname(program_dir)
  if (dest_name == "/Applications"):
    dest_name = "the Applications folder"
  pkg_root = "pkg_root"
  if (op.exists(pkg_root)):
    if (not options.overwrite):
      print("ERROR: pkg_root already exists - run with --overwrite to ignore")
      return False
    shutil.rmtree(pkg_root)
  os.mkdir(pkg_root)

  base_name = "%s-%s" % (options.package_name.lower(), options.version)
  base_pkg = "%s.pkg"%base_name
  pkg_name = os.path.join(options.dist_dir, "%s-%s.pkg" %(base_name, options.machine_type))
  pkg_id = "%s.%s" % (options.organization, base_name)
  os.mkdir("resources")
  welcome = open("resources/welcome.txt", "w")
  welcome.write("""\
This package contains %(package)s version %(version)s compiled for %(arch)s. \
It will install the full command-line suite and graphical launcher(s) \
to %(dest_name)s.  If you need to install %(package)s to a different \
location, you must use the command-line installer (also available from \
our website).""" % { "package" : options.package_name,
                     "version" : options.version,
                     "arch"    : system_type,
                     "dest_name" : dest_name })
  # identify .app bundles in top-level directory
  app_bundle_xml = []
  for file_name in os.listdir(program_dir):
    if (file_name.endswith(".app")):
      full_path = op.join(program_dir, file_name)
      if (not op.isdir(full_path)) : continue
      app_bundle_xml.append("""\
  <dict>
    <key>BundleHasStrictIdentifier</key>
    <true/>
    <key>BundleIsRelocatable</key>
    <true/>
    <key>BundleIsVersionChecked</key>
    <true/>
    <key>BundleOverwriteAction</key>
    <string>upgrade</string>
    <key>RootRelativeBundlePath</key>
    <string>%s</string>
  </dict>
""" % full_path[1:])
  plist_file = "%s.plist" % options.package_name
  plist = open(plist_file, "w")
  plist.write("""\
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<array>
%s
</array>
</plist>
""" % "\n".join(app_bundle_xml))
  plist.close()
  # move directory
  # XXX os.path.relpath behavior is not consistent between Python versions -
  # the behavior we need here was not introduced until 2.7.  Since the program
  # directory is assumed to be absolute, we just chop off the leading '/'.
  rel_path = op.dirname(program_dir)[1:]
  pkg_path = op.join("pkg_root", rel_path)
  os.mkdir(pkg_path)
  shutil.move(program_dir, pkg_path)
  # write out distribution.xml
  misc_files = []
  if (options.background is not None):
    if (options.background.endswith(".jpg")):
      shutil.copyfile(options.background, "resources/background.jpg")
      misc_files.append(
        """<background file="background.jpg" mime-type="image/jpeg" """+
        """alignment="topleft" scaling="proportional"/>""")
    else :
      assert options.background.endswith(".png")
      shutil.copyfile(options.background, "resources/background.png")
      misc_files.append(
        """<background file="background.png" mime-type="image/png" """+
        """alignment="topleft" scaling="proportional"/>""")
  if (options.license_file is not None):
    shutil.copyfile(options.license_file, "resources/license.txt")
    misc_files.append(
      """<license    file="license.txt"    mime-type="text/plain" />""")
  distfile = open("distribution.xml", "w")
  distfile.write("""\
<?xml version="1.0" encoding="utf-8" standalone="no"?>
<installer-gui-script minSpecVersion="1">
    <title>%(package)s %(version)s</title>
    <organization>%(org)s</organization>
    <domains enable_localSystem="true"/>
    <options customize="never" require-scripts="true" rootVolumeOnly="true" />
    <!-- Define documents displayed at various steps -->
    <welcome    file="welcome.txt"    mime-type="text/plain" />
    %(misc_files)s
    <pkg-ref id="%(pkg_id)s"
             version="0"
             auth="root">%(base_pkg)s</pkg-ref>
    <choices-outline>
        <line choice="%(pkg_id)s"/>
    </choices-outline>
    <choice
        id="%(pkg_id)s"
        visible="false"
        title="%(package)s"
        description="%(package)s"
        start_selected="true">
      <pkg-ref id="%(pkg_id)s"/>
    </choice>
</installer-gui-script>
""" % { "package"  : options.package_name,
        "version"  : options.version,
        "pkg_id"   : pkg_id,
        "base_pkg" : base_pkg,
        "org"      : options.organization,
        "misc_files" : "\n".join(misc_files) })
  distfile.close()

  print("Fixing package permissions:", pkg_root, file=out)
  call(['chmod','-R','0755',pkg_root])

  # Run packaging commands
  pkg_args = [
    "pkgbuild",
    "--ownership", "recommended",
    "--root", "pkg_root",
    "--identifier", pkg_id,
    "--component-plist", plist_file,
    base_pkg,
  ]
  print("Calling pkgbuild:", pkg_args, file=out)
  call(pkg_args, out)
  product_args = [
    "productbuild",
    "--distribution", "distribution.xml",
    "--resources", "resources",
    "--package-path", ".",
    "--version", options.version,
    pkg_name,
  ]
  print("Calling productbuild:", product_args, file=out)
  call(product_args, out)
  assert op.exists(pkg_name)
  return True

if (__name__ == "__main__"):
  if (not run(sys.argv[1:])):
    sys.exit(1)
