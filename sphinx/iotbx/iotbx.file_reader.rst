======================================
iotbx.file_reader - generic file input
======================================

The ``iotbx.file_reader`` module is intended to provide a single entry point
for reading most common crystallographic file formats.  This allows the
programmer to use the underlying input functions without needing to know the
specific APIs in detail, although the resulting objects will still be
format-specific.  It is also designed to support automatic file type
determination, first by guessing the format based on the file extension, then
by trying a succession of input methods until one finishes without an error.
This facility is used both for processing command-line arguments (especially
via the :py:mod:`iotbx.phil` extensions), and for handling file input in the
Phenix GUI.

In the simplest case, reading a file requires only a single method::

  >>> from iotbx.file_reader import any_file
  >>> input_file = any_file(sys.argv[1:])
  >>> file_data = input_file.file_object

Note that if the extension does not imply a particular format to try first, or
if parsing using the appropriate module fails due to corrupted file data, this
may be more inefficient than explicitly specifying the file
type, and should be used only when the format is not known in advance.  You
can alternately specify which input module to use::

  >>> pdb_in = any_file("model.pdb", force_type="pdb")
  >>> pdb_in.assert_file_type("pdb")
  >>> hierarchy = pdb_in.file_object.construct_hierarchy()

  >>> mtz_in = any_file("data.mtz", force_type="hkl")
  >>> miller_arrays = mtz_in.file_server.miller_arrays

This will skip the automatic format detection and only try the specified
input method.  Several options are available for error handling; the default
behavior when ``force_type`` is set is to pass through any exceptions
encountered when calling the underlying input method::

  >>> from iotbx.file_reader import any_file
  >>> f = any_file("model.pdb", force_type="hkl")
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "/Users/nat/src/cctbx_project/iotbx/file_reader.py", line 197, in any_file
      raise_sorry_if_not_expected_format=raise_sorry_if_not_expected_format)
    File "/Users/nat/src/cctbx_project/iotbx/file_reader.py", line 251, in __init__
      read_method()
    File "/Users/nat/src/cctbx_project/iotbx/file_reader.py", line 341, in _try_as_hkl
      assert (hkl_file.file_type() is not None), "Not a valid reflections file."
  AssertionError: Not a valid reflections file.

This is fine for internal use when an unexpected file parsing error is likely
to be a bug in the code, but less suitable when processing user input.
Alternately, a :py:class:`libtbx.utils.Sorry` exception may be raised instead::

  >>> f = any_file("model.pdb", force_type="hkl", raise_sorry_if_errors=True)
  Sorry: Couldn't read 'model.pdb' as file type 'hkl': Not a valid reflections file.

For PDB, MTZ, and CIF files (the most commonly used formats in macromolecular
crystallograph), it is also possible to get similar behavior by treating the
file extension as an implicit replacement for ``force_type``::

  >>> from iotbx.file_reader import any_file
  >>> f = any_file("model.sca")
  >>> print f.file_type
  pdb
  >>> f = any_file("model.sca", raise_sorry_if_not_expected_format=True)
  Sorry: File format error:
  Not a valid reflections file.

The allowed file types are specified in the module::

  standard_file_descriptions = {
    'pdb'  : "Model",
    'hkl'  : "Reflections",
    'cif'  : "Restraints",
    'seq'  : "Sequence",
    'xplor_map'  : "XPLOR map",
    'ccp4_map' : "CCP4 map",
    'phil' : "Parameters",
    'xml'  : "XML",
    'pkl'  : "Python pickle",
    'txt'  : "Text",
    'mtz'  : "Reflections (MTZ)",
    'aln'  : "Sequence alignment",
    'hhr'  : "HHpred alignment",
    'img'  : "Detector image",
  }

However, in most cases only a subset of these will be tried automatically.

API documentation
=================

.. automodule:: iotbx.file_reader
    :members:
    :undoc-members:
    :show-inheritance:
