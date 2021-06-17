from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
from cctbx import uctbx
from cctbx import crystal
from cctbx.array_family import flex
import cctbx.sgtbx.lattice_symmetry
from scitbx import matrix
import math, sys
import scitbx.math
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from boost_adaptbx.boost import rational
from libtbx.math_utils import ifloor
from six.moves import cStringIO as StringIO
from six.moves import range
from six.moves import zip

def divisor(n, pairs=False):
  """ find all divisors """
  i_n = int( n )
  d_n = float(n)
  w_n = int( math.floor( math.sqrt(n) ) )
  result = []
  for trial_div in range(1, w_n+1):
    if d_n%trial_div == 0:
      a = trial_div
      b = i_n//trial_div
      if not pairs:
        result.append( a )
        if a != b:
          result.append( b )
      if pairs:
        result.append( (a,b) )

  return result


def find_triples_slowly(n):
  """ finds all triples for which a*b*c==n"""
  # do this by finding all triplets for which
  # the above condition holds
  #
  # lets do this the stupid way, i.e. like this:
  divs = divisor(n)
  triples = []
  for fst in divs:
    for snd in divs:
      for trd in divs:
        tmp = fst*snd*trd
        if tmp == n:
          trip = (fst,snd,trd)
          if not trip in triples:
            triples.append( trip )
  return triples

def sort_triple( triple ):
  order = flex.sort_permutation( flex.int(triple),False )
  sorted_triple = ( triple[order[0]],
                    triple[order[1]],
                    triple[order[2]] )

  return sorted_triple


def make_permutations( triple ):
  # no nonsense hard coded combinations
  c1 = (triple[0], triple[1], triple[2])
  c2 = (triple[0], triple[2], triple[1])
  c3 = (triple[1], triple[0], triple[2])
  c4 = (triple[1], triple[2], triple[0])
  c5 = (triple[2], triple[0], triple[1])
  c6 = (triple[2], triple[1], triple[0])
  permuts = [c1,c2,c3,c4,c5,c6 ]
  result = []
  for item in permuts:
    if not item in result:
      result.append( item )

  return result


def make_triples(n):
  """ finds all triples for which a*b*c==n"""
  # do this by finding all triplets for which
  # the above condition holds
  #
  divs = divisor(n,True)
  triples = []
  for item in divs:
    div_a = divisor(item[0],True)
    div_b = divisor(item[1],True)
    for a in div_a:
      trial_trip = ( a[0], a[1], item[1] )
      trial_trip = sort_triple( trial_trip )
      if not trial_trip in triples:
        triples.append( trial_trip )
    for b in div_b:
      trial_trip = ( b[0], b[1], item[0] )
      trial_trip = sort_triple( trial_trip )
      if not trial_trip in triples:
        triples.append( trial_trip )

  return triples

def find_triples( n ):
  # fist make a list of ordered triples
  triples = make_triples( n )
  result = []
  # now we have get all combinations
  for trip in triples:
    result += make_permutations( trip )
  return result

def generate_matrix( order ):
  """ make all sub lattices generators of order order """
  # first get all triples
  triples = find_triples( order )
  matrices = []
  for triple in triples:
    # make a list of d and e values
    a = triple[0]
    b = triple[1]
    c = triple[2]
    d_and_e = []
    f_list = []

    if a%2 == 0:
      # a is even
      tmp = -a//2+1
      while tmp <= a//2:
        d_and_e.append( tmp )
        tmp +=1
    if a%2 != 0 :
      # a is odd
      tmp = -(a-1)//2
      while tmp <= (a-1)//2:
        d_and_e.append( tmp )
        tmp += 1

    if b%2 == 0:
      # b is even
      tmp = -b//2+1
      while tmp <= b//2:
        f_list.append( tmp )
        tmp +=1

    if b%2 != 0:
      # b is odd
      tmp = -(b-1)//2
      while tmp <= (b-1)//2:
        f_list.append( tmp )
        tmp += 1

    for d in d_and_e :
      for e in d_and_e :
        for f in f_list :
          mat = [rational.int(a),
                 rational.int(d),
                 rational.int(e),
                 rational.int(0),
                 rational.int(b),
                 rational.int(f),
                 rational.int(0),
                 rational.int(0),
                 rational.int(c)]
          matrices.append(  matrix.sqr(mat)  )
  return matrices


def generate_matrix_up_to_order( end_order,start_order=1 ):
  cum = 0
  all_matrices = []
  for n in range(int(start_order),int(end_order)+1):
    all_matrices += generate_matrix( n )
  return all_matrices

def generate_cb_op_up_to_order( end_order, start_order=1):
  def string_it( row ):
    result = ""
    abc = ['a','b','c']
    part = 0
    for n,j in zip(row,abc):
      nn=""
      if n != 0:
        part += 1
        if n >0:
          if part==1:
            if n==1:
              nn = j
            else:
              nn = str(n)+j
          if part > 1:
            if n==1:
              nn = "+"+j
            else:
              nn = "+"+str(n)+j
        if n < 0 :
          if part==1:
            if n==-1:
              nn = "-"+j
            else:
              nn = str(n)+j
          else:
            nn = str(n)+j
      result += nn
    return result
  mats = generate_matrix_up_to_order(end_order,start_order)
  cb_ops = []
  for mat in mats:
    new_a = [mat[0], mat[3], mat[6] ] # mat[0:3]
    new_b = [mat[1], mat[4], mat[7] ] # mat[3:6]
    new_c = [mat[2], mat[5], mat[8] ]
    tmp = string_it(new_a)+","+string_it(new_b)+","+string_it(new_c)
    cb_ops.append( tmp )
  return cb_ops


def rt_mx_as_rational(rot_mat):
  # make sure one provides an integer matrix!
  tmp_mat = rot_mat.r().as_double()
  rational_list = []
  for ij in tmp_mat:
    rational_list.append( rational.int( ifloor(ij) ) )
  return matrix.sqr( rational_list )

def generate_inverse_matrix_up_to_order( end_order,start_order=1 ):
  mat_list = generate_matrix_up_to_order( end_order, start_order )
  inv_mat_list = []
  for mat in mat_list:
    inv_mat  = mat.inverse()
    inv_mat_list.append( inv_mat )
  return( inv_mat_list )

class make_list_of_target_xs_up_to_order(object):
  def __init__(self,
               basic_xs,
               order=1,
               max_delta=2):

    self.basic_xs = basic_xs
    self.basic_xs_n = self.basic_xs.change_basis( self.basic_xs.change_of_basis_op_to_niggli_cell() )
    self.basic_to_niggli_cb_op = self.basic_xs.change_of_basis_op_to_niggli_cell()


    self.basis = matrix.sqr( self.basic_xs_n.unit_cell().orthogonalization_matrix() )
    self.order = order
    self.max_delta = max_delta

    self.matrices = generate_matrix_up_to_order( self.order )
    self.cb_ops = generate_cb_op_up_to_order( self.order )

    self.xs_list = []
    self.sg_list = []

    self.extra_cb_op = []
    for mat, cb_op in zip(self.matrices, self.cb_ops):
      self.make_new_xs( mat, cb_op )

  def make_new_xs(self,
                  mat,
                  cb_op,
                  to_reference=True):
    # make new lattice
    new_basis = self.basis*mat.as_float()
    new_uc = uctbx.unit_cell( orthogonalization_matrix = new_basis )

    tmp_xs = crystal.symmetry( unit_cell=new_uc,
                               space_group=sgtbx.lattice_symmetry.group(new_uc,self.max_delta),
                               assert_is_compatible_unit_cell=False,
                             )

    extra_cb_op = tmp_xs.change_of_basis_op_to_reference_setting()
    self.extra_cb_op.append( extra_cb_op )
    #
    new_sg = None
    try:
      new_sg = sgtbx.space_group_info( group=self.basic_xs_n.space_group()).change_basis( cb_op )
    except Exception: pass

    if to_reference:
      tmp_xs = tmp_xs.change_basis( extra_cb_op )
      try:
        new_sg = new_sg.change_basis( extra_cb_op )
      except Exception: pass

    self.xs_list.append( tmp_xs )
    self.sg_list.append( new_sg )

class compare_lattice(object):
  def __init__(self,
               xs_a,
               xs_b,
               max_delta=2.0,
               out=None,
               relative_length_tolerance=0.05,
               absolute_angle_tolerance=10.0,
               order=1):

    self.relative_length_tolerance = relative_length_tolerance
    self.absolute_angle_tolerance = absolute_angle_tolerance
    self.order=order

    self.max_delta=max_delta
    self.out = out
    if self.out is None:
      self.out = sys.stdout

    self.xs_a = xs_a
    self.xs_b = xs_b
    # Go to niggli cell
    self.xs_a_n = self.xs_a.niggli_cell()
    self.xs_b_n = self.xs_b.niggli_cell()

    if ( self.xs_a.niggli_cell().unit_cell().volume() > self.xs_b.niggli_cell().unit_cell().volume() ):
      self.xs_a = xs_b
      self.xs_b = xs_a
      # go to niggli cell please
      self.xs_a_n = self.xs_a.niggli_cell()
      self.xs_b_n = self.xs_b.niggli_cell()


    volume_ratio   = self.xs_b_n.unit_cell().volume() / self.xs_a_n.unit_cell().volume()
    n_volume_ratio = math.floor( volume_ratio + 1 )

    self.show_basic_info()


    self.basis_a = matrix.sqr( self.xs_a_n.unit_cell().orthogonalization_matrix() )
    print("Cartesian basis (column) vectors of lego cell:", file=self.out)
    tmp_bas = self.basis_a.as_list_of_lists()
    print("  / %5.1f %5.1f %5.1f \\  " %(tmp_bas[0][0], tmp_bas[0][1], tmp_bas[0][2]), file=self.out)
    print("  | %5.1f %5.1f %5.1f |  " %(tmp_bas[1][0], tmp_bas[1][1], tmp_bas[1][2]), file=self.out)
    print("  \\ %5.1f %5.1f %5.1f /  " %(tmp_bas[2][0], tmp_bas[2][1], tmp_bas[2][2]), file=self.out)
    print(file=self.out)
    self.basis_b = matrix.sqr( self.xs_b_n.unit_cell().orthogonalization_matrix() )
    print("Cartesian basis (column) vectors of target cell:", file=self.out)
    tmp_bas = self.basis_b.as_list_of_lists()
    print("  / %5.1f %5.1f %5.1f \\  " %(tmp_bas[0][0], tmp_bas[0][1], tmp_bas[0][2]), file=self.out)
    print("  | %5.1f %5.1f %5.1f |  " %(tmp_bas[1][0], tmp_bas[1][1], tmp_bas[1][2]), file=self.out)
    print("  \\ %5.1f %5.1f %5.1f /  " %(tmp_bas[2][0], tmp_bas[2][1], tmp_bas[2][2]), file=self.out)
    print(file=self.out)

    self.lattice_symm_a = sgtbx.lattice_symmetry.group( self.xs_a_n.unit_cell(),self.max_delta )
    self.lattice_symm_b = sgtbx.lattice_symmetry.group( self.xs_b_n.unit_cell(),self.max_delta )

    self.xs_ref = crystal.symmetry( unit_cell=self.xs_b_n.unit_cell(),
                                    space_group=self.lattice_symm_b,
                                    assert_is_compatible_unit_cell=False )

    # make a list of matrices
    self.matrix_list = generate_matrix_up_to_order( n_volume_ratio, max(n_volume_ratio-1,1) )
    self.uc_and_symm_list=[]
    self.possible_solutions = []
    print("A total of %4i matrices in the hermite normal form have been generated."%( len( self.matrix_list ) ), file=self.out)
    print("The volume changes they cause lie between %4i and %4i."%(n_volume_ratio, max(n_volume_ratio-1,1)), file=self.out)

    count =0

    print(file=self.out)
    print("Trying all matrices", file=self.out)
    print(file=self.out)
    print("   1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0", file=self.out)
    print("  ", end=' ', file=self.out)
    for mat in self.matrix_list:
      count+=1
      found_it=False
      tmp_xs = self.make_new_cell_and_symmetry( mat )

      self.uc_and_symm_list.append( tmp_xs )

      tmp_gen = self.xs_b_n.unit_cell().similarity_transformations(
        tmp_xs[2].unit_cell(),
        self.relative_length_tolerance,
        self.absolute_angle_tolerance,
        self.order)

      if tmp_gen.size()>0:
        found_it = True
        cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx( sgtbx.rot_mx(tmp_gen[0]))).inverse()
        tmp_sol = (mat,cb_op,tmp_xs[2].change_basis( cb_op ) , tmp_xs[3]) # matrix, corresponding cb_op, cell + lattice sym
        self.possible_solutions.append( tmp_sol )
        #self.show_solution(tmp_sol)
      if found_it:
        print("*", end=' ', file=self.out)
        self.out.flush()
      else:
        if count%10==0:
          print("|", end=' ', file=self.out)
          self.out.flush()
        else:
          print(".", end=' ', file=self.out)
          self.out.flush()
      if count%20==0:
        print(file=self.out)
        print("  ", end=' ', file=self.out)
        self.out.flush()

    print(file=self.out)
    print(" Listing all possible solutions", file=self.out)
    count=0
    if len(self.possible_solutions)==0:
      print(file=out)
      print("No relations found for this particular sublattice of the specified target cell.", file=out)
      print("   ", file=out)
      print(file=out)

    for tmp_sol in self.possible_solutions:
      count+=1
      print(file=self.out)
      print("Solution %4i"%(count), file=self.out)
      self.show_solution( tmp_sol )


  def show_basic_info(self):
    print(file=self.out)
    print("Crystal symmetries in supplied setting", file=self.out)
    print(file=self.out)
    print("Target crystal symmetry:", file=self.out)
    self.xs_b.show_summary(f=self.out,
                           prefix="    ")
    print("Building block crystal symmetry: ", file=self.out)
    self.xs_a.show_summary(f=self.out,
                           prefix="    ")
    print(file=self.out)

    print("Crystal symmetries in Niggli setting", file=self.out)
    print(file=self.out)
    print("Target crystal symmetry:", file=self.out)
    self.xs_b_n.show_summary(f=self.out,
                           prefix="    ")
    print("Building block (lego cell) crystal symmetry: ", file=self.out)
    self.xs_a_n.show_summary(f=self.out,
                           prefix="    ")
    print(file=self.out)
    print("Volume ratio between target and lego cell: %5.2f"%( self.xs_b_n.unit_cell().volume() / self.xs_a_n.unit_cell().volume() ), file=self.out)
    print(file=self.out)


  def show_solution(self, sol_entry):
    mat = sol_entry[0].as_list_of_lists()
    print("--------------------------------------------------------------", file=self.out)
    print("Target unit cell :     %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f"%(self.xs_b_n.unit_cell().parameters()[0],
                                                                            self.xs_b_n.unit_cell().parameters()[1],
                                                                            self.xs_b_n.unit_cell().parameters()[2],
                                                                            self.xs_b_n.unit_cell().parameters()[3],
                                                                            self.xs_b_n.unit_cell().parameters()[4],
                                                                            self.xs_b_n.unit_cell().parameters()[5]), file=self.out)
    print("Lego cell :            %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f"%(self.xs_a_n.unit_cell().parameters()[0],
                                                                            self.xs_a_n.unit_cell().parameters()[1],
                                                                            self.xs_a_n.unit_cell().parameters()[2],
                                                                            self.xs_a_n.unit_cell().parameters()[3],
                                                                            self.xs_a_n.unit_cell().parameters()[4],
                                                                            self.xs_a_n.unit_cell().parameters()[5]), file=self.out)
    print(file=self.out)
    print("               /%4i %4i %4i  \\  "%(mat[0][0],mat[0][1],mat[0][2]), file=self.out)
    print("matrix :  M =  |%4i %4i %4i  |  "%(mat[1][0],mat[1][1],mat[1][2]), file=self.out)
    print("               \\%4i %4i %4i  /  "%(mat[2][0],mat[2][1],mat[2][2]), file=self.out)
    print(file=self.out)
    print("Additional Niggli transform:     ", sol_entry[3].as_xyz(), file=self.out)
    print("Additional similarity transform: ", sol_entry[1].as_xyz(), file=self.out)
    print("Resulting unit cell :  %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f"%(
      sol_entry[2].unit_cell().parameters()[0],
      sol_entry[2].unit_cell().parameters()[1],
      sol_entry[2].unit_cell().parameters()[2],
      sol_entry[2].unit_cell().parameters()[3],
      sol_entry[2].unit_cell().parameters()[4],
      sol_entry[2].unit_cell().parameters()[5]), file=self.out)
    print("Deviations :           %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f"%(
      100*(self.xs_b_n.unit_cell().parameters()[0]-sol_entry[2].unit_cell().parameters()[0])/self.xs_b_n.unit_cell().parameters()[0],
      100*(self.xs_b_n.unit_cell().parameters()[1]-sol_entry[2].unit_cell().parameters()[1])/self.xs_b_n.unit_cell().parameters()[1],
      100*(self.xs_b_n.unit_cell().parameters()[2]-sol_entry[2].unit_cell().parameters()[2])/self.xs_b_n.unit_cell().parameters()[2],
      (self.xs_b_n.unit_cell().parameters()[3]-sol_entry[2].unit_cell().parameters()[3]),
      (self.xs_b_n.unit_cell().parameters()[4]-sol_entry[2].unit_cell().parameters()[4]),
      (self.xs_b_n.unit_cell().parameters()[5]-sol_entry[2].unit_cell().parameters()[5])  ), file=self.out)
    print("Deviations for unit cell lengths are listed in %.", file=self.out)
    print("Angular deviations are listed in degrees.", file=self.out)
    print(file=self.out)
    print(" ", file=self.out)
    print("--------------------------------------------------------------", file=self.out)

  def make_new_cell_and_symmetry(self,
                                 mat):
    # make new lattice
    new_basis = self.basis_a*mat.as_float()
    new_uc = uctbx.unit_cell( orthogonalization_matrix = new_basis )

    tmp_xs = crystal.symmetry( new_uc, "P1" ) # if latice symm is supplied, we loose track of cb op sequence!

    extra_cb_op = sgtbx.change_of_basis_op( tmp_xs.change_of_basis_op_to_niggli_cell().as_xyz() )

    # get the niggli cell please
    new_uc = new_uc.niggli_cell()
    # get the lattice symmetry please
    tmp_xs = tmp_xs.change_basis( extra_cb_op )
    lattice_group = tmp_xs.space_group_info()
    return( (new_uc,
             lattice_group,
             tmp_xs,
             extra_cb_op) )

def tst_compare():
  uc1=uctbx.unit_cell( '61.28,95.92,145.02,90,90,90' )
  xs1 = crystal.symmetry(uc1, "P212121")
  uc2=uctbx.unit_cell('115.5,149.0,115.60,90,115.3,90' )
  xs2 = crystal.symmetry(uc2, "P1211")
  out=StringIO()
  compare_object= compare_lattice(xs1,xs2,order=1,out=out)
  assert len( compare_object.possible_solutions )==1

def tst_sublattice():
  # compare results to table 2, Acta Cryst A36, 242-248, (1980) Billiet, Rolley Le-Coz
  i = [2, 3,  4,  5,  6,  7,  8,   9,   10 ]
  N = [7, 13, 35, 31, 91, 57, 155, 130, 217]
  for ii, iN in zip(i,N):
    tmp = generate_matrix( ii )
    assert ( len(tmp) == iN )

def tst_make_bigger_cell():
  uc1=uctbx.unit_cell( '61.28,95.92,145.02,90,90,90' )
  xs1 = crystal.symmetry(uc1, "P212121")
  tmp  = make_list_of_target_xs_up_to_order(xs1,order=5)
  for xs,mat,cb_op in zip(tmp.xs_list,
                          tmp.matrices,
                          tmp.extra_cb_op):
    ratio = xs.unit_cell().volume()/xs1.unit_cell().volume()
    det = float(mat.determinant())
    det_ref= float(cb_op.c().r().determinant())
    assert approx_equal( ratio/(det/det_ref),1.0, eps=0.001 )

def tst_centered_cells():
  out = StringIO()
  uc1 = uctbx.unit_cell( '13,19,23,90,97,90' )
  xs1 = crystal.symmetry( uc1, "C2" )
  cb_op = xs1.change_of_basis_op_to_primitive_setting()
  xs1p = xs1.change_basis( cb_op )
  xs2 = crystal.symmetry( uc1, "P1" )
  comp = compare_lattice( xs1, xs2, out=out )
  assert len( comp.possible_solutions )==1

  uc1 = uctbx.unit_cell( '13,13,23,90,90,120' )
  xs1 = crystal.symmetry( uc1, "R32:H" )
  cb_op = xs1.change_of_basis_op_to_primitive_setting()
  xs1p = xs1.change_basis( cb_op )
  xs2 = crystal.symmetry( uc1, "P1" )
  comp = compare_lattice( xs1, xs2, out=out )
  assert len( comp.possible_solutions )==1




def exercise():
  tst_sublattice()
  tst_compare()
  tst_make_bigger_cell()
  tst_centered_cells()
  print(format_cpu_times())

if (__name__ == "__main__"):
  exercise()
