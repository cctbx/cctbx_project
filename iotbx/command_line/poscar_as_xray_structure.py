def run(args):
  assert len(args) > 0
  import iotbx.poscar
  for file_name in args:
    print file_name
    poscar = iotbx.poscar.reader(
      lines=open(file_name).read().splitlines(),
      source_info=file_name)
    poscar.make_up_types_if_necessary()
    poscar.xray_structure(u_iso=0.05).show_summary().show_scatterers()
    print

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
