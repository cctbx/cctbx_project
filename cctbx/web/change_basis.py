def number_from_string(s):
  flds = s.split("/")
  if (len(flds) == 2):
    from boost.rational import int as rint
    return rint(int(flds[0]), int(flds[1]))
  return (eval(s)+0)*1

def nums_from_string(string):
  for c in ",;|{}()[]":
    string = string.replace(c, " ")
  flds = string.split()
  nums = []
  for fld in flds:
    try:
      num = number_from_string(fld)
    except Exception:
      return None
    nums.append(num)
  return nums

def rt_from_string(string, default_r_identity=False, r_den=12**2, t_den=12**3):
  from cctbx import sgtbx
  from scitbx import matrix
  result = None
  s = string.strip()
  if (len(s) == 0):
    s = "a,b,c"
  try:
    cb_op = sgtbx.change_of_basis_op(
      symbol=s, stop_chars="", r_den=r_den, t_den=t_den)
  except RuntimeError, e:
    pass
  else:
    return cb_op.c_inv().as_rational()
  nums = nums_from_string(s)
  if (nums is not None):
    if (len(nums) == 12):
      if (string.find("|") < 0):
        result = matrix.rt((nums[:9], nums[9:]))
      else:
        result = matrix.rt(
          (nums[:3]+nums[4:7]+nums[8:11],
          (nums[3], nums[7], nums[11])))
    elif (len(nums) == 9):
      result = matrix.rt((nums[:9], (0,0,0)))
    elif (len(nums) == 3):
      if (default_r_identity):
        r = (1,0,0,0,1,0,0,0,1)
      else:
        r = (0,0,0,0,0,0,0,0,0)
      result = matrix.rt((r, nums))
  return result

def raise_uninterpretable(what, expression):
  from libtbx.utils import Sorry
  raise Sorry(
    'Uninterpretable expression for %s: "%s"' % (
      what, expression.replace("\n", " ").replace("\r", " ")))

def p_from_string(string):
  p = rt_from_string(string, default_r_identity=True)
  if (p is None):
    raise_uninterpretable("change-of-basis matrix", string)
  return p

def w_from_string(string):
  w = rt_from_string(string)
  if (w is None):
    raise_uninterpretable("symmetry matrix", string)
  return w

def xyz_from_string(string):
  nums = nums_from_string(string)
  if (nums is None or len(nums) != 3):
    raise_uninterpretable("coordinates", string)
  return tuple(nums)

def display_r(r):
  for ir in xrange(3):
    print "  (%8s %8s %8s)" % (r(ir,0), r(ir,1), r(ir,2))

def display_rt(rt):
  r, t = rt.r, rt.t.elems
  for ir in xrange(3):
    print "  (%8s %8s %8s | %8s)" % (r(ir,0), r(ir,1), r(ir,2), t[ir])

def interpret_form_data(form):
  from cctbx.web import cgi_utils
  inp = cgi_utils.inp_from_form(form,
    (("p_or_q", "P"),
     ("p_transpose", "off"),
     ("cb_expr", ""),
     ("obj_type", "xyz"),
     ("obj_expr", "")))
  return inp

def run(server_info, inp, status):
  print "<pre>"
  from scitbx import matrix
  p = p_from_string(string=inp.cb_expr)
  assert inp.p_or_q in ["P", "Q"]
  if (inp.p_or_q == "Q"):
    p = p.inverse()
  assert inp.p_transpose in ["off", "on"]
  if (inp.p_transpose == "on"):
    p = matrix.rt((p.r.transpose(), p.t))
  print "P:"
  display_rt(p)
  print
  q = p.inverse()
  print "Q:"
  display_rt(q)
  print
  if (len(inp.obj_expr.strip()) != 0):
    if (inp.obj_type in ["xyz", "hkl"]):
      triple = xyz_from_string(string=inp.obj_expr)
      if (inp.obj_type == "xyz"):
        print "Transformation law: (Q,q) xyz"
        print
        print "  xyz:", triple
        print
        print "  xyz':", (
          q.r.as_float() * matrix.col(triple) + q.t.as_float()).elems
        print
      else:
        print "Transformation law: hkl P"
        print
        print "  hkl:", triple
        print
        print "  hkl':", (matrix.row(triple) * p.r.as_float()).elems
        print
    elif (inp.obj_type == "unit_cell"):
      from cctbx import uctbx
      uc = uctbx.unit_cell(inp.obj_expr)
      print "Transformation law: Pt G P"
      print
      print "unit cell:", uc
      print
      g = matrix.sym(sym_mat3=uc.metrical_matrix())
      print "metrical matrix:"
      display_r(g)
      print
      gp = p.r.transpose().as_float() * g * p.r.as_float()
      print "metrical matrix':"
      print
      display_r(gp)
      print
      ucp = uctbx.unit_cell(metrical_matrix=gp.as_sym_mat3())
      print "unit cell':", ucp
      print
    elif (inp.obj_type == "Ww"):
      w = w_from_string(string=inp.obj_expr)
      print "Transformation law: (Q,q) (W,w) (P,p)"
      print
      print "(W, w):"
      display_rt(w)
      print
      wp = q.as_float() * w.as_float() * p.as_float()
      print "(W, w)':"
      display_rt(wp)
      print
    else:
      raise RuntimeError("Unknown obj_type: %s" % inp.obj_type)
  print "</pre>"
