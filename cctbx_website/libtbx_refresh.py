from __future__ import absolute_import, division, print_function
from cctbx_website.command_line.extract_script_from_html_cctbx_doc import run
import libtbx.load_env

parent_dir = libtbx.env.dist_path("cctbx_website", default=None)
if parent_dir is not None:
  run(parent_dir)
