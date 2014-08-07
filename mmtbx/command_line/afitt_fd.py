from mmtbx.geometry_restraints.afitt import finite_difference_test


def run(pdb_file, cif_file, ligand_names,atom,scale,verbose):
  finite_difference_test(pdb_file, cif_file, ligand_names,atom,scale,verbose)

if (__name__ == "__main__"):
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("pdb_file", help="pdb file")
  parser.add_argument("cif_file", help="cif file", default=0)
  parser.add_argument("ligand_names", help="3-letter ligand names separated by commas")
  parser.add_argument("atom", help="i_seq of atom whose x-coord will be "
                                   "modified (to test AFITT this should "
                                   "be a ligand atom)")
  parser.add_argument( "-s", "--scale", default=1, help="weight applied" \
                      "to AFITT gradient/traget", type=float)
  parser.add_argument('-v', dest='verbose', action='store_true')
  parser.epilog='Example: phenix.python afitt_fd.py vAla3.pdb vAla3.cif NME 37'
  args = parser.parse_args()
  ligand_names=args.ligand_names.split(',')
  print args.verbose
  run(args.pdb_file, args.cif_file, ligand_names, int(args.atom),
      args.scale, args.verbose)

