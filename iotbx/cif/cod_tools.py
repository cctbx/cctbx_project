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
    self.hkl_cif_pairs = {}
    self.hkl = {}
    self.cif = {}
    if (cod_ids is None or len(cod_ids)==0):
      if ext is None or ext == "hkl":
        for sub_dir in sorted(os.listdir(hkl_dir)):
          if (sub_dir.startswith(".")): continue
          hkl_sub_dir = op.join(hkl_dir, sub_dir)
          if (os.path.isfile(hkl_sub_dir)): continue
          for node in sorted(os.listdir(hkl_sub_dir)):
            if (node.startswith(".")): continue
            if (not node.endswith(".hkl")): continue
            cod_id = node[:-4]
            hkl_path = op.join(hkl_sub_dir, node)
            self.hkl.setdefault(cod_id, hkl_path)
      if ext is None or ext == "cif":
        for sub_dir in sorted(os.listdir(cif_dir)):
          if (sub_dir.startswith(".")): continue
          cif_sub_dir = op.join(cif_dir, sub_dir)
          for node in sorted(os.listdir(cif_sub_dir)):
            if (node.startswith(".")): continue
            if (not node.endswith(".cif")): continue
            cod_id = node[:-4]
            cif_path = op.join(cif_dir, sub_dir, cod_id+".cif")
            self.cif.setdefault(cod_id, cif_path)
            if (cod_id in self.hkl):
              self.hkl_cif_pairs.setdefault(cod_id, (self.hkl[cod_id], cif_path))
    else:
      n_missing_all = 0
      for cod_id in cod_ids:
        hkl_path = op.join(hkl_dir, cod_id[0], cod_id+".hkl")
        cif_path = op.join(cif_dir, cod_id[0], cod_id+".cif")
        n_missing = 0
        if (op.isfile(cif_path)):
          self.cif.setdefault(cod_id, cif_path)
        else:
          print "Missing COD cif file:", cif_path
          n_missing += 1
        if (op.isfile(hkl_path)):
          self.hkl.setdefault(cod_id, hkl_path)
        else:
          print "Missing COD hkl file:", hkl_path
          n_missing += 1
        if (n_missing == 0):
          self.hkl_cif_pairs.setdefault(cod_id, (hkl_path, cif_path))
        else:
          n_missing_all += n_missing
      if (n_missing_all != 0):
        raise RuntimeError("Number of missing COD files: %d" % n_missing_all)

  def show_summary(self, out=None):
    if out is None: out = sys.stdout
    print >> out, "Number of hkl without cif:", len(
      set(self.hkl_cif_pairs.keys())-set(self.hkl.keys()))
    print >> out, "Number of cif: ", len(self.cif)
    print >> out, "Number of hkl: ", len(self.hkl)
    print >> out, "Number of hkl+cif:", len(self.hkl_cif_pairs)
