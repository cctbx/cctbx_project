from __future__ import absolute_import, division, print_function
import os, sys
op = os.path

class build_hkl_cif(object):
  def __init__(self, cod_ids=None, ext=None):
    if ext is not None:
      assert ext in ('cif','hkl')
    envar = "COD_SVN_WORKING_COPY"
    cod_svn = os.environ.get(envar)
    if (cod_svn is None):
      msg = [
        "Environment variable %s not defined:" % envar,
        "  Usage:",
        "    mkdir /some/path",
        "    cd /some/path",
        "    svn checkout svn://www.crystallography.net/cod",
        "    export %s=/some/path/cod" % envar]
      raise RuntimeError("\n".join(msg))
    cif_dir = op.join(cod_svn, "cif")
    hkl_dir = op.join(cod_svn, "hkl")
    def cod_path(cod_dir, cod_id, ext):
      return op.join(cod_dir, cod_id[0], cod_id[1:3], cod_id[3:5], cod_id+ext)
    self.hkl_cif_pairs = {}
    self.hkl = {}
    self.cif = {}
    if (cod_ids is None or len(cod_ids)==0):
      if ext is None or ext == "hkl":
        for root, dirs, files in os.walk(hkl_dir):
          if '.svn' in dirs: dirs.remove('.svn')
          for file_name in sorted(files):
            if (file_name.startswith(".")): continue
            if (not file_name.endswith(".hkl")): continue
            cod_id = file_name[:-4]
            hkl_path = op.join(root, file_name)
            self.hkl.setdefault(cod_id, hkl_path)
      if ext is None or ext == "cif":
        for root, dirs, files in os.walk(cif_dir):
          if '.svn' in dirs: dirs.remove('.svn')
          for file_name in sorted(files):
            if (file_name.startswith(".")): continue
            if (not file_name.endswith(".cif")): continue
            cod_id = file_name[:-4]
            cif_path = op.join(root, file_name)
            self.cif.setdefault(cod_id, cif_path)
            if (cod_id in self.hkl):
              self.hkl_cif_pairs.setdefault(
                cod_id, (self.hkl[cod_id], cif_path))
    else:
      n_missing_all = 0
      for cod_id in cod_ids:
        hkl_path = cod_path(hkl_dir, cod_id, ".hkl")
        cif_path = cod_path(cif_dir, cod_id, ".cif")
        n_missing = 0
        if (op.isfile(cif_path)):
          self.cif.setdefault(cod_id, cif_path)
        else:
          print("Missing COD cif file:", cif_path)
          n_missing += 1
        if (op.isfile(hkl_path)):
          self.hkl.setdefault(cod_id, hkl_path)
        else:
          print("Missing COD hkl file:", hkl_path)
          n_missing += 1
        if (n_missing == 0):
          self.hkl_cif_pairs.setdefault(cod_id, (hkl_path, cif_path))
        else:
          n_missing_all += n_missing
      if (n_missing_all != 0):
        raise RuntimeError("Number of missing COD files: %d" % n_missing_all)

  def show_summary(self, out=None):
    if out is None: out = sys.stdout
    print("Number of hkl without cif:", len(
      set(self.hkl_cif_pairs.keys())-set(self.hkl.keys())), file=out)
    print("Number of cif: ", len(self.cif), file=out)
    print("Number of hkl: ", len(self.hkl), file=out)
    print("Number of hkl+cif:", len(self.hkl_cif_pairs), file=out)
