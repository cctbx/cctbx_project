"""Preserve backwards compatibility"""
from __future__ import absolute_import, division, print_function
from mmtbx.f_model.f_model_info import *

'''
This file is specifically to preserve backwards compatibility with earlier
versions of the Phenix GUI where the class for storing results in a pickle,
phenix.runtime.result_summary, has an import statement for this file.
Pickling an instance of result_summary stored the location of this file in the
cctbx hierarchy. The import statement in the result_summary class has been
moved outside of the class definition to prevent future issues. Moved to f_model.py
'''
