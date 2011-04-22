def str_ev(ev):
  return "[%d,%d,%d]" % ev

class literal_description:

  def __init__(self,symop):
    self.symop = symop # an instance of rt_mx
    from cctbx import sgtbx
    self.r_info = sgtbx.rot_mx_info(self.symop.r())
    self.t_info = sgtbx.translation_part_info(self.symop)

  def long_form(self):

    if str(self.t_info.origin_shift()) == "0,0,0":
      origin_message = ""
    else:
      origin_message = "origin shift: %11s" % (str(self.t_info.origin_shift()))

    if origin_message == "" and str(self.t_info.intrinsic_part()) == "0,0,0":
      intrinsic_message = ""
    else:
      intrinsic_message = "intrinsic part: %11s" % (str(self.t_info.intrinsic_part()))

    if (self.r_info.type() == 1):
      return "     no rotation             %s %s" % (
        intrinsic_message,
        origin_message)
    elif (self.r_info.type() == -1):
      return "     rotation type: %d %s %s" % (
        self.r_info.type(),
        intrinsic_message,
        origin_message)
    elif (abs(self.r_info.type()) == 2):
      return "     %d-fold about %10s %s %s" % (
        self.r_info.type(),
        str_ev(self.r_info.ev()),
        intrinsic_message,
        origin_message)
    else:
      return "%s %d-fold about %10s %s %s" % (
         {-1:"rev.",1:"fwd."}[self.r_info.sense()],
         self.r_info.type(),
         str_ev(self.r_info.ev()),
         intrinsic_message,
         origin_message)

  def cosets_form(self):
    return "%20s  %20s   Rotation: %4s ; direction: %10s ; screw/glide: %10s"%(
          self.symop,
          self.symop.r().as_hkl(),
          self.symop.r().info().type() ,
          self.symop.r().info().ev(),
          "("+self.symop.t().as_string()+")" )

  def labelit_check_pdb_symmetry_form(self):
    return "%22s  %18s  %s"%(
          self.symop,
          self.symop.r().as_hkl(),
          self.long_form())

  def labelit_check_pdb_symmetry_short(self):
    return "%18s  %s"%(
          self.symop.r().as_hkl(),
          self.long_form())

  def select(self,format):
    return apply(literal_description.__dict__[format],[self])
