
from libtbx import adopt_init_args

# XXX prototype for container for residue-specific validation information.
# not currently used anywhere, but eventually this can be used to display
# multi-criterion overviews of a structure (much more cleanly than in the
# current PHENIX GUI code).
class residue_info (object) :
  def __init__ (self,
                resname,
                resid,
                chain_id,
                altloc,
                xyz,
                sec_str=None,
                ramalyze_score=None,
                rotalyze_score=None,
                cbetadev_score=None,
                clash_max=None,
                bonds_max=None,
                angles_max=None,
                dihedral_max=None,
                chiral_max=None,
                plane_max=None,
                rscc=None,
                b_mean=None,
                mean_2fofc=None,
                mean_fc=None) :
    adopt_init_args(self, locals())

  def is_ramachandran_outlier (self) :
    if (self.ramalyze_score is None) : return None
    else : return (self.ramalyze_score < 0.002)

  def is_ramachandran_favored (self) :
    if (self.ramalyze_score is None) : return None
    else : return (self.ramalyze_score > 0.02)

  def is_rotalyze_outlier (self) :
    if (self.rotalyze_score is None) : return None
    else : return (self.rotalyze_score < 1.0)

  def is_clash_outlier (self) :
    return (self.clash_max >= 0.4)

  def is_cbeta_outlier (self) :
    return (self.cbetadev_score > 0.25)

  def get_atom_selection (self) :
    if (self.altloc == "") :
      return "(chain '%s' and resname %s and resid %s)" % (self.chain_id,
        self.resname, self.resid)
    else :
      return "(chain '%s' and resname %s and resid %s and altloc %s)" % \
        (self.chain_id, self.resname, self.resid, self.altloc)
