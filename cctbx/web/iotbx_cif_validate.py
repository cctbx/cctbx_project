from __future__ import absolute_import, division, print_function
from cctbx.web import cgi_utils
import iotbx.cif
from iotbx.cif import validation
import iotbx.pdb.mmcif

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("cif_file", None),
     ("cif_text", None),
     ("cif_dic", None),
     ("extract_miller_arrays", False),
     ("extract_crystal_structures", False),
     ("extract_pdb_hierarchy", False),
     ("diffraction_index_equivalent", None)))

  inp.extract_miller_arrays = inp.extract_miller_arrays == "True"
  inp.extract_crystal_structures = inp.extract_crystal_structures == "True"
  inp.extract_pdb_hierarchy = inp.extract_pdb_hierarchy == "True"
  return inp


def run_implementation(server_info, inp, status):
  cif_text = inp.cif_text
  if cif_text is None:
    cif_text = inp.cif_file

  reader = iotbx.cif.reader(input_string=cif_text, raise_if_errors=False)
  if reader.error_count():
    print("Errors encountered during parsing")
    return
  else:
    print("No parsing errors.")

  if inp.cif_dic is not None:
    cif_dic = validation.smart_load_dictionary(name=inp.cif_dic)
    print()
    print("Validating CIF against %s:" %inp.cif_dic)
    cif_model = reader.model()
    error_handler = cif_model.validate(cif_dic, show_warnings=True)
    if len(error_handler.warnings) + len(error_handler.errors) == 0:
      print("No validation errors found.")
    print()

  if inp.extract_miller_arrays:
    print("Extracting Miller arrays:")
    try:
      miller_arrays = reader.as_miller_arrays()
    except iotbx.cif.builders.CifBuilderError as e:
      print("CifBuilderError: %s" %str(e))
    else:
      for ma in miller_arrays:
        print(ma)
        ma.show_comprehensive_summary()
      if len(miller_arrays) == 0:
        print("No Miller arrays found.")
      print()

  if inp.extract_crystal_structures:
    print("Extracting crystal structures:")
    try:
      crystal_structures = reader.build_crystal_structures()
    except iotbx.cif.builders.CifBuilderError as e:
      print("CifBuilderError: %s" %str(e))
    else:
      for xs in crystal_structures.values():
        xs.show_summary().show_scatterers()
        print()
      if len(crystal_structures) == 0:
        print("No crystal structures found.")
    print()

  if inp.extract_pdb_hierarchy:
    print("Extracting pdb.hierarchy:")
    # TODO am I a dictionary ?- if so change me accordingly..
    hierarchy = iotbx.pdb.mmcif.pdb_hierarchy_builder(
      cif_model.blocks.values()[0]).hierarchy
    hierarchy.show()
    print()


def run(server_info, inp, status):
  print("<pre>")
  run_implementation(server_info=server_info, inp=inp, status=status)
  print("</pre>")
