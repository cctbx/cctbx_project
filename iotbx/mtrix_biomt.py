"""
Routines to read and write PDB-style BIOMTR records
"""

from __future__ import absolute_import, division, print_function
from scitbx import matrix
import iotbx.pdb
from libtbx.utils import Sorry
from six import string_types
from six.moves import range, zip
from cctbx.array_family import flex

class container(object):

  def __init__(self):
    self.r=[]
    self.t=[]
    self.coordinates_present=[]
    self.serial_number=[]

  @classmethod
  def from_lists(cls, rl, tl, snl, cpl):
    assert len(rl) == len(tl) == len(snl) == len(cpl)
    result = cls()
    for r,t,sn,cp in zip(rl,tl,snl,cpl):
      ignore_transform = r.is_r3_identity_matrix() and t.is_col_zero()
      result.add(
        r=r, t=t,
        coordinates_present=(cp or ignore_transform),
        serial_number=sn)
    return result

  def add(self, r, t, serial_number, coordinates_present=False):
    self.r.append(r)
    self.t.append(t)
    self.coordinates_present.append(coordinates_present)
    self.serial_number.append(serial_number)

  def validate(self, eps=1e-4):
    raise_sorry = False
    for i, r in enumerate(self.r):
      if(not r.is_r3_rotation_matrix(rms_tolerance=eps)):
        raise_sorry = True
        print ('  ERROR: matrix with serial number %s is not proper' % self.serial_number[i])
    if raise_sorry:
      raise Sorry("One or more rotation matrices is not proper. See the log for details.")
    present = False
    for (r,t,n,cp) in zip(self.r,self.t,self.serial_number,
                            self.coordinates_present):
      if(not (r.is_r3_identity_matrix() and t.is_col_zero())):
        present = present or cp
    return present

  def as_pdb_string(self):
    return format_MTRIX_pdb_string(
      rotation_matrices=self.r,
      translation_vectors=self.t,
      serial_numbers=self.serial_number,
      coordinates_present_flags=self.coordinates_present)

  def format_BIOMT_pdb_string(self):
    '''
    BIOMT data sample
    REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
    '''
    lines = []
    fmt1="REMARK 350   BIOMT1  %2d%10.6f%10.6f%10.6f     %10.5f"
    fmt2="REMARK 350   BIOMT2  %2d%10.6f%10.6f%10.6f     %10.5f"
    fmt3="REMARK 350   BIOMT3  %2d%10.6f%10.6f%10.6f     %10.5f"
    for sn_, r_, t_ in zip(self.serial_number, self.r, self.t):
      lines.append(fmt1%(int(sn_), r_[0],r_[1],r_[2], t_[0]))
      lines.append(fmt2%(int(sn_), r_[3],r_[4],r_[5], t_[1]))
      lines.append(fmt3%(int(sn_), r_[6],r_[7],r_[8], t_[2]))
    return "\n".join(lines)

  def is_empty(self):
    return len(self.r) == 0

def parse_MTRIX_BIOMT_records_cif(cif_block, recs='mtrix'):
  # this is temporarily work-around. Whole thing for matrices need to be
  # rewritten to be able to carry all the info from mmCIF, see
  # https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies#Anchor-Biol

  assert recs in ['mtrix', 'biomt']
  block_name = '_struct_ncs_oper' if recs=='mtrix' else '_pdbx_struct_oper_list'
  rots = []
  trans = []
  serial_number = []
  coordinates_present = []
  ncs_oper = cif_block.get('%s.id' % block_name)
  if ncs_oper is not None:
    for i,sn in enumerate(ncs_oper):
      # filter everything for X0 and P here because they represent some other translations
      # not related to whole molecule reproduction in BIOMT:
      # http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_pdbx_struct_oper_list.id.html
      if recs=='biomt' and (sn == "P" or sn == "X0" or sn == "H"):
        continue
      serial_number.append((sn,i))
      if recs == 'mtrix':
        code_loop_or_item = cif_block.get('%s.code' % block_name)
        if isinstance(code_loop_or_item, flex.std_string):
          coordinates_present.append(code_loop_or_item[i] == 'given')
        else:
          coordinates_present.append(code_loop_or_item == 'given')
      else:
        coordinates_present.append(False) # no way to figure out
      r = [(cif_block.get('%s.matrix[%s][%s]' %(block_name,x,y)))
        for x,y in ('11', '12', '13', '21', '22', '23', '31','32', '33')]
      if not isinstance(r[0], string_types):
        r = [elem[i] for elem in r]
      try:
        rots.append(matrix.sqr([float(r_elem) for r_elem in r]) )
        t = [(cif_block.get('%s.vector[%s]' %(block_name, x)))
          for x in '123']
        if not isinstance(t[0], string_types):
          t = [elem[i] for elem in t]
        trans.append(matrix.col([float(t_elem) for t_elem in t]))
      except ValueError:
        raise Sorry("Error in %s information. Likely '?' instead of a number." % block_name)
  if recs=='mtrix':
    # sort records by serial number
    serial_number.sort()
    items_order = [i for (_,i) in serial_number]
    trans = [trans[i] for i in items_order]
    rots = [rots[i] for i in items_order]
    coordinates_present = [coordinates_present[i] for i in items_order]
    serial_number = [j for (j,_) in serial_number]
  return rots, trans, serial_number, coordinates_present

def process_BIOMT_records_cif(cif_block):
  rots, trans, serial_number, coordinates_present = parse_MTRIX_BIOMT_records_cif(cif_block, 'biomt')
  return container.from_lists(rots, trans, serial_number, coordinates_present)

def process_MTRIX_records_cif(cif_block):
  rots, trans, serial_number, coordinates_present = parse_MTRIX_BIOMT_records_cif(cif_block, 'mtrix')
  return container.from_lists(rots, trans, serial_number, coordinates_present)

def process_BIOMT_records_pdb(lines):
  '''(pdb_data,boolean,float) -> group of lists
  extract REMARK 350 BIOMT information, information that provides rotation matrices
  and translation  data, required for generating  a complete multimer from the asymmetric unit.

  BIOMT data sample:
  REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
  REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
  REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
  REMARK 350   BIOMT1   2  0.559048 -0.789435  0.253492        0.30000
  REMARK 350   BIOMT2   2  0.722264  0.313528 -0.616470        0.00100
  REMARK 350   BIOMT3   2  0.407186  0.527724  0.745457        0.05000

  The data from every 3 lines will be combined to a list
  [[x11,x12,x13,x21,x22,x23,x31,x32,x33],[x1,x2,x3]]
  the first component, with the 9 numbers are the rotation matrix,
  and the second component, with the 3 numbers is the translation information
  x is the serial number of the operation

  for x=2 the transformation_data will be
  [[0.559048,-0.789435,0.253492,0.722264,0.313528,-0.616470,0.407186,0.527724,0.745457][0.30000,0.00100,0.05000]]

  the result is a list of libtx group_arg constructs, each one contains
    pdb_inp.process_BIOMT_records()[1].values              # give a list containing float type numbers
    pdb_inp.process_BIOMT_records()[1].coordinates_present # True when transformatin included in pdb file
    pdb_inp.process_BIOMT_records()[1].serial_number       # is an integer

  @author: Youval Dar (2013)
  '''
  source_info = lines # XXX
  if not source_info:
    return container()                # check if any BIOMT info is available
    # collecting the data from the remarks. Checking that we are collecting only data
    # and not part of the remarks header by verifying that the 3rd component contains "BIOMT"
    # and that the length of that component is 6
  biomt_data = [[float(_x_elem) for _x_elem in  x.split()[3:]] for x in source_info if (
    x.split()[2].find('BIOMT') > -1) and (len(x.split()[2]) == 6)]
  # test that there is no missing data
  if len(biomt_data)%3 != 0:
    raise RuntimeError(
      "Improper or missing set of PDB BIOMAT records. Data length = %s" % \
      str(len(biomt_data)))
  # test that the length of the data match the serial number, that there are no missing records
  # temporary workaround, could be plain text over there instead of
  # expected number of records, see 5l93
  try:
    temp = int(source_info[-1].split()[3])
  except ValueError:
    temp = 0
  if len(biomt_data)/3.0 != temp:
    raise RuntimeError(
      "Missing record sets in PDB BIOMAT records \n" + \
      "Actual number of BIOMT matrices: {} \n".format(len(biomt_data)/3.0)+\
      "expected according to serial number: {} \n".format(temp))
  rots = []
  trans = []
  serial_number = []
  coordinates_present = []
  for i in range(len(biomt_data)//3):
    # i is the group number in biomt_data
    # Each group composed from 3 data sets.
    j = 3*i;
    rotation_data = \
      biomt_data[j][1:4] + biomt_data[j+1][1:4] + biomt_data[j+2][1:4]
    translation_data = \
      [biomt_data[j][-1], biomt_data[j+1][-1], biomt_data[j+2][-1]]
    rots.append(matrix.sqr(rotation_data))
    trans.append(matrix.col(translation_data))
    # For BIOMT the first transform is the identity matrix
    ignore_transform = rots[-1].is_r3_identity_matrix() and trans[-1].is_col_zero()
    coordinates_present.append(ignore_transform)
    serial_number.append(i+1)
  return container.from_lists(rots, trans, serial_number, coordinates_present)

def process_MTRIX_records_pdb(lines):
  """
  Read MTRIX records from a pdb file
  """
  storage = {}
  for line in lines:
    if (line.startswith("MTRIX") and line[5:6] in ["1", "2", "3"]):
      r = read_mtrix_record(line=line)
      stored = storage.get(r.serial_number)
      if (stored is None):
        values = [[None]*9,[None]*3]
        done = [None]*3
        present = [None]*3
        # Create a dictionary record for a serial number
        stored = storage[r.serial_number] = (values, done, present)
      else:
        values, done, present = stored
      for i_col,v in enumerate(r.r):
        values[0][(r.n-1)*3+i_col] = v
      values[1][r.n-1] = r.t
      done[r.n-1] = r.n
      present[r.n-1] = r.coordinates_present
  rots = []
  trans = []
  serial_number = []
  coordinates_present = []
  for sn in sorted(storage.keys()):
    values, done, present = storage[sn]
    serial_number.append(sn)
    rots.append(matrix.sqr(values[0]))
    trans.append(matrix.col(values[1]))
    if (sorted(done) != [1,2,3] or len(set(present)) != 1):
      raise RuntimeError("Improper set of PDB MTRIX records")
    ignore_transform = rots[-1].is_r3_identity_matrix() and trans[-1].is_col_zero()
    coordinates_present.append(present[0]==1 or ignore_transform)
  return container.from_lists(rots, trans, serial_number, coordinates_present)

def format_MTRIX_pdb_string(rotation_matrices, translation_vectors,
      serial_numbers=None, coordinates_present_flags=None):
  '''
  MTRIX data sample
  REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
  MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
  MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
  MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
  MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
  MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
  MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
  '''
  assert len(rotation_matrices) == len(translation_vectors)
  if(serial_numbers is None): serial_numbers = range(0, len(rotation_matrices))
  if(coordinates_present_flags is None):
    coordinates_present_flags = [False]*len(rotation_matrices)
  lines = []
  fmt1="MTRIX1  %2d%10.6f%10.6f%10.6f     %10.5f    %s"
  fmt2="MTRIX2  %2d%10.6f%10.6f%10.6f     %10.5f    %s"
  fmt3="MTRIX3  %2d%10.6f%10.6f%10.6f     %10.5f    %s"
  for sn_, r_, t_, p_ in zip(serial_numbers, rotation_matrices,
                             translation_vectors, coordinates_present_flags):
    flag = " "
    if p_: flag="1"
    lines.append(fmt1%(int(sn_), r_[0],r_[1],r_[2], t_[0], flag))
    lines.append(fmt2%(int(sn_), r_[3],r_[4],r_[5], t_[1], flag))
    lines.append(fmt3%(int(sn_), r_[6],r_[7],r_[8], t_[2], flag))
  return "\n".join(lines)

class read_mtrix_record(iotbx.pdb.read_scale_record):

  __slots__ = iotbx.pdb.read_scale_record.__slots__ + [
    "serial_number", "coordinates_present"]

  def __init__(O, line):
    iotbx.pdb.read_scale_record.__init__(O, line=line, source_info="")
    O.serial_number = line[7:10]
    O.coordinates_present = (len(line) >= 60 and line[59] != " ")
