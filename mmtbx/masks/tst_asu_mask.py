import sys, os, time
import iotbx.pdb
import iotbx.pdb.remark_2_interpretation
from scitbx.array_family import flex
import iotbx
import StringIO
from mmtbx import masks
import cctbx
from cctbx import miller
from cctbx.sgtbx import space_group_info
from cctbx.development import random_structure
from cctbx import crystal
from cctbx import xray
from libtbx.test_utils import approx_equal, is_below_limit
from libtbx.utils import format_cpu_times
from cctbx.development import debug_utils

# cStringIO does not work
#  import cStringIO
# cout = cStringIO.StringIO("\n")
cout = StringIO.StringIO()

# modified cctbx.sgtbx.space_group_info.any_compatible_unit_cell
# from cctbx/sgtbx/__init__.py
def random_compatible_unit_cell(self, volume=None, asu_volume=None):
  import random
  from cctbx import uctbx
  assert [volume, asu_volume].count(None) == 1
  if (volume is None):
    volume = asu_volume * self.group().order_z()
  sg_number = self.type().number()
  rnd = []
  for i in xrange(6):
    rnd.append( random.random() )
    if( i<3 ):
      rnd[i] *= 4.0
  if   (sg_number <   3):
    alpha = 5.0 + 150.0*rnd[3]
    beta = 5.0 + min(360.0-alpha-35.0, 165.0)*rnd[4]
    gamma_min = 5.0 + max(alpha,beta)-min(alpha,beta)
    gamma_max = min(min(alpha+beta,175.0),360.0-alpha-beta-6.0)
    assert gamma_max >= gamma_min
    gamma = gamma_min + (gamma_max-gamma_min)*rnd[5]
    assert alpha>=5.0 and beta>=5.0 and gamma>=5.0 \
        and (alpha+beta+gamma) <= 355.0
    params = (1.+rnd[0], 1.+rnd[1], 1.+rnd[2], alpha, beta, gamma)
  elif (sg_number <  16):
    params = (1.0+rnd[0], 1.0+rnd[1], 1.0+rnd[2], 90, 90.5+80.0*rnd[5], 90)
  elif (sg_number <  75):
    params = (1.+rnd[0], 1.+rnd[2], 1.+rnd[1], 90, 90, 90)
  elif (sg_number < 143):
    params = (1.+rnd[0], 1.+rnd[0], 1.+rnd[2], 90, 90, 90)
  elif (sg_number < 195):
    params = (1.+rnd[2], 1.+rnd[2], 1.+rnd[1], 90, 90, 120)
  else:
    params = (1., 1., 1., 90, 90, 90)
  unit_cell = uctbx.unit_cell(params).change_basis(
    cb_op=self.change_of_basis_op_to_reference_setting().inverse())
  f = (volume / unit_cell.volume())**(1/3.)
  params = list(unit_cell.parameters())
  for i in xrange(3): params[i] *= f
  return uctbx.unit_cell(params)

# overriding 'any cell' with 'random cell'
cctbx.sgtbx.space_group_info.any_compatible_unit_cell_original = \
    cctbx.sgtbx.space_group_info.any_compatible_unit_cell
cctbx.sgtbx.space_group_info.any_compatible_unit_cell = \
    random_compatible_unit_cell


def compare_fc(obs, other, tolerance = 1.0E-9):
  assert obs.is_complex_array()
  assert other.is_complex_array(), other.__class__
  matching = miller.match_indices(obs.indices(), other.indices())
  data0 = obs.select(matching.pairs().column(0)).data()
  data = other.select(matching.pairs().column(1)).data()
  assert data0.size() == data.size(), str(data0.size()) + " != " \
      + str(data.size())
  assert data.size() > 1, str(data.size())
  max_rel_dif = 0.0
  max_dif = 0.0
  max_mx = 0.0
  for i in xrange(data.size()):
    dif = abs(data[i]-data0[i])
    mx = max( abs(data[i]),abs(data0[i]) )
    if mx > tolerance*1.0E-2:
      rel_dif = dif / mx
    else:
      rel_dif = 0.0
    if rel_dif > max_rel_dif:
      max_rel_dif = rel_dif
      max_dif = dif
      max_mx = mx
  assert ((max_rel_dif <= tolerance) or (max_mx <= tolerance*1.0E-2)), \
    "max  rel_dif = "+ str(max_rel_dif)+ "   dif = "+str(max_dif)+"    mx =" \
    +str(max_mx)
  return data.size() # max_rel_dif


def get_radii(structure):
  from cctbx.eltbx import van_der_waals_radii
  unknown = []
  atom_radii = []
  for i_seq, scatterer in enumerate(structure.scatterers()):
    try:
      atom_radii.append(
           van_der_waals_radii.vdw.table[scatterer.element_symbol()])
    except:
      unknown.append(scatterer.element_symbol())
  return atom_radii


SpaceGroups = ("P 21 21 21", "P 21", "P 1 1 21", "P 21/n", "P1",
    "Fm3m", "R3", "P61", "I41", "I 1 m 1")
Elements = ("N", "C", "O", "H", "Ca", "C", "B", "Li", "Ru", "N", "H", "H",
    "Mg", "Se")
def make_atoms(n_atoms):
  assert n_atoms>0
  atoms = []
  for i in xrange(n_atoms):
    if i < len(Elements) :
      atoms.append( Elements[i] )
    else:
      atoms.append("C")
  return atoms

def build_struc(spgr_symbol, n, atom_volume):
  group = cctbx.sgtbx.space_group_info( spgr_symbol )
  cell = group.any_compatible_unit_cell_original( asu_volume = n * atom_volume )
  symmetry = crystal.symmetry(unit_cell=cell,
                              space_group_symbol=spgr_symbol)
  structure = xray.structure(crystal_symmetry=symmetry)
  for i in xrange(n):
    if i < len(Elements):
      element = Elements[i]
    else:
      element = "C"
    site = ( (i%(n//2))/float(n), (n-i%(n//3))/float(n), (i%(n//4))/float(n) )
    scatterer = xray.scatterer(
                   site = site,
                   u = 0.1,
                   occupancy = 1.0,
                   scattering_type = element)
    structure.add_scatterer(scatterer)
  return structure


def zero_test(asu_mask, fc, tolerance = 1.0E-9):
  radii = []
  sites = []
  assert len(radii) == len(sites)
  asu_mask.compute( sites, radii )
  fm_asu = asu_mask.structure_factors( fc.indices() )
  fm_asu = fc.set().array( data = fm_asu )
  max_zero = flex.max( flex.abs(fm_asu.data()) )
  assert isinstance(max_zero, float), max_zero.__class__
  assert max_zero < tolerance, "Maximum deviation from zero = "+str(max_zero)

def check_group(group):
  ops = group.smx()
  ltrs = group.ltr()
  print "NNN = ", ltrs.__class__
  cb = group.type().cb_op()
  print "CB r= ", cb.c().r().as_double()
  print "CB t= ", cb.c().t()
  ident = cctbx.sgtbx.rt_mx()
  for ltr in ltrs:
    print "LTR= ", ltr
  for op in ops:
    print "Op: ", op, " Order=", op.r().order(), "  r.den= ", op.r().den(), \
        " tr.den= ", op.t().den()


def compare_masks(struc, opts):
  tolerance = opts.tolerance
  resolution = opts.resolution
  solvent_radius = opts.solvent_radius
  shrink_radius = opts.shrink_radius
  verbose = opts.verbose
  cout.truncate(0)
  time_p1 = 0.0
  time_asu = 0.0
  time_orig = 0.0
  params = masks.mask_master_params.extract()
  params.ignore_hydrogens = False
  params.ignore_zero_occupancy_atoms = False
  params.solvent_radius = solvent_radius
  params.shrink_truncation_radius = shrink_radius
  fc = struc.structure_factors(d_min = resolution).f_calc()
  while fc.data().size() <= 3 :
    resolution /= 1.2
    assert resolution > 1.0E-3
    fc = struc.structure_factors( d_min = resolution).f_calc()
  print >>cout, "Resolution= ", resolution, "  solvent radius= ", \
      solvent_radius, "  shrink radius= ", shrink_radius,  "  Tolerance= ", \
      tolerance, "  Number of reflection= ", fc.data().size()
  struc.show_summary(cout)
  group = struc.space_group()
  print >>cout, "Cell volume= ", struc.unit_cell().volume(), \
    "  Group order= ", group.order_z(), " p= ", group.order_p()
  print >>cout, "Hall  symbol: ", group.type().hall_symbol()
  #check_group(group)

  tb = time.time()
  asu_mask = masks.atom_mask(
      unit_cell = struc.unit_cell(),
      group = struc.space_group(),
      resolution = fc.d_min(),
      grid_step_factor = params.grid_step_factor,
      solvent_radius = params.solvent_radius,
      shrink_truncation_radius = params.shrink_truncation_radius )
  te = time.time()
  time_asu += (te-tb)
  grid =  asu_mask.grid_size()
  print >>cout, "asu mask grid = ", grid
  zero_test(asu_mask, fc, tolerance = tolerance)
  radii = get_radii(struc)
  assert len(radii) == len(struc.sites_frac())
  tb = time.time()
  asu_mask.compute( struc.sites_frac(), radii )
  te = time.time()
  time_asu += (te-tb)
  print >>cout, "   n asu atoms= ", asu_mask.n_asu_atoms(), \
      "   has-enclosed= ", asu_mask.debug_has_enclosed_box
  tb = time.time()
  fm_asu = asu_mask.structure_factors( fc.indices() )
  fm_asu = fc.set().array( data = fm_asu )
  te = time.time()
  time_asu_sf = te-tb
  time_asu += (te-tb)
  # save files
  if not opts.save_files is None:
    tmp_file = open(opts.save_files + ".pdb", "w")
    print >>tmp_file, struc.as_pdb_file()
    tmp_file.close()
    asu_mask.xplor_write_map(opts.save_files + "_mask.map")
    asu_mask.xplor_write_map(opts.save_files + "_inverted_mask.map", 1, True)
    # also save structure factors
    import iotbx.mtz
    mtzo = iotbx.mtz.object()
    mtzo.set_title("mask test")
    mtzo.add_history(line="start")
    mtzo.set_space_group_info(fm_asu.space_group_info())
    mtzo.set_hkl_base(fm_asu.unit_cell())
    crystal = mtzo.add_crystal(
      name="mask_test_crystal",
      project_name="mask_test_project",
      unit_cell=fm_asu.unit_cell())
    dataset = crystal.add_dataset(
      name="mask_test_dataset",
      wavelength=1)
    assert dataset.add_miller_array(
      miller_array=fm_asu,
      column_root_label="F",
      #column_types=column_types
      ) is dataset
    mtzo.add_history(line="done")
    mtzo.write(opts.save_files + "_sf.mtz")
  #
  # ========= old mask =============
  #
  tb = time.time()
  struc_p1 = struc.expand_to_p1()
  te = time.time()
  time_p1_exp = (te-tb)
  time_p1 += (te-tb)
  fc_p1 = fc.deep_copy()
  fc_p1 = fc_p1.customized_copy(crystal_symmetry = struc_p1.crystal_symmetry())
  tb = time.time()
  blk_p1 = masks.bulk_solvent(
    xray_structure = struc_p1,
    gridding_n_real = grid,
    ignore_zero_occupancy_atoms = params.ignore_zero_occupancy_atoms,
    ignore_hydrogen_atoms = params.ignore_hydrogens,
    solvent_radius = params.solvent_radius,
    shrink_truncation_radius = params.shrink_truncation_radius)
  te = time.time()
  time_p1_msk = (te-tb)
  time_p1 += (te-tb)
  tb = time.time()
  fm_p1 = blk_p1.structure_factors( miller_set = fc_p1 )
  te = time.time()
  time_p1_sf = (te-tb)
  time_p1 += (te-tb)
  blk_p1.show_summary(cout)
  ### original mask
  tb = time.time()
  blk_o = masks.bulk_solvent(
    xray_structure = struc,
    gridding_n_real = grid,
    ignore_zero_occupancy_atoms = params.ignore_zero_occupancy_atoms,
    ignore_hydrogen_atoms = params.ignore_hydrogens,
    solvent_radius = params.solvent_radius,
    shrink_truncation_radius = params.shrink_truncation_radius)
  te = time.time()
  time_orig_msk = (te-tb)
  time_orig += (te-tb)
  tb = time.time()
  fm_o = blk_o.structure_factors( miller_set = fc )
  te = time.time()
  time_orig_sf = (te-tb)
  time_orig += (te-tb)
  print >>cout, "Number of reflections ::: Fm asu = ", fm_asu.data().size(), \
    "Fm P1 = ", fm_p1.data().size()
  print >>cout, "Time ( ms )    P1= ", time_p1*1000.0, "   orig= ", \
      time_orig*1000.0, "    asu= ", time_asu*1000.0
  print >>cout, "Times ( ms ) mask_asu= ", asu_mask.debug_mask_asu_time, \
      " atoms_to_asu= ", asu_mask.debug_atoms_to_asu_time, \
      " accessible= ", asu_mask.debug_accessible_time, \
      " contact= ", asu_mask.debug_contact_time, \
      " Fc= ", time_asu_sf*1000.0, \
      " fft= ", asu_mask.debug_fft_time
  print >>cout, "Times ( ms ) orig:  mask= ", time_orig_msk*1000.0, "  Fc=", \
      time_orig_sf*1000.0
  print >>cout, "Times ( ms ) p1 :  expand= ", time_p1_exp*1000.0, "  mask= ", \
      time_p1_msk*1000.0, "  Fc=", time_p1_sf*1000.0
  assert fm_asu.data().size() == fm_o.data().size()
  t_v1 = asu_mask.contact_surface_fraction
  t_v2 = blk_p1.contact_surface_fraction
  t_v3 = max( abs(t_v1), abs(t_v2) )
  if t_v3 > 1.0E-6:
    t_v4 = abs( t_v1 - t_v2 ) / t_v3
  else:
    t_v4 = 0.0
  if( t_v4>1.0E-6 ):
    if not opts.failed_file is None:
      tmp_file = open(opts.failed_file, "w")
      print >>tmp_file, struc.as_pdb_file()
      tmp_file.close()
    raise "Not equal solvent volume"

  assert approx_equal(
    asu_mask.contact_surface_fraction, blk_p1.contact_surface_fraction)
  assert approx_equal(
    asu_mask.accessible_surface_fraction, blk_p1.accessible_surface_fraction)
  assert is_below_limit(
    value=asu_mask.accessible_surface_fraction,
    limit=asu_mask.contact_surface_fraction)
  n_compared = compare_fc(fm_asu, fm_p1, tolerance = tolerance)
  assert n_compared == fm_asu.data().size(), \
    "N compared refls: "+str(n_compared) + " != " + str(fm_asu.data().size())
  assert n_compared >0
  if verbose:
    print cout.getvalue()
  # test that second calculation will produce the same results
  asu_mask.compute( struc.sites_frac(), radii )
  fm_asu2 = asu_mask.structure_factors( fc.indices() )
  fm_asu2 = fc.set().array( data = fm_asu2 )
  n_compared = compare_fc(fm_asu, fm_asu2, tolerance = tolerance)
  assert n_compared == fm_asu.data().size(), \
    "N compared refls: "+str(n_compared) + " != " + str(fm_asu.data().size())
  cout.truncate(0)


def standard_tests(groups, options):
  if options.verbose:
    print "Standard tests, n space groups = ", len(groups), "\n options="
    print options, "\n"
  solvent_radius = options.solvent_radius
  shrink_radius = options.shrink_radius
  for sg in groups:
    for islv in xrange(3):
      slv_rad = solvent_radius*islv
      for ishr in xrange(3):
        shr_rad = shrink_radius*0.5*ishr
        struc = build_struc(sg, options.n_atoms,  options.atom_volume)
        options.solvent_radius = slv_rad
        options.shrink_radius = shr_rad
        compare_masks(struc, options)


def random_tests(groups, opts):
  import random
  atoms = make_atoms(opts.n_atoms)
  resolution = opts.resolution
  print "Number of space groups: ", len(groups)
  print "Number of random tests per space group: ", opts.random, "\n"
  for sg in groups:
    print "Space group= ", sg, "  n tests= ", opts.random
    group = space_group_info(sg)
    print "       HM= ", group.type().universal_hermann_mauguin_symbol(), \
        "  LOOKUP= ", group.type().lookup_symbol(), "  HALL= ", \
        group.type().hall_symbol()
    for i in xrange(opts.random):
      if i==0 :
        slv_rad = 1.1
        shr_rad = 0.9
        res = resolution
      elif i==1 :
        slv_rad = 1.1
        shr_rad = 0.0
        res = resolution
      elif i==2 :
        slv_rad = 0.0
        shr_rad = 0.0
        res = resolution
      else:
        slv_rad = 3.0*random.random()
        shr_rad = 1.33333*slv_rad*random.random()
        res = resolution + random.random()
      struc = None
      try:
        # occationally this fails with groups: P6522, P3112, P4322, P4122, Pma2
        # and small number of atoms
        struc = random_structure.xray_structure(
            space_group_info = group,
            volume_per_atom = opts.atom_volume,
            general_positions_only = False, #True,
            elements = atoms
            )
      except:
        print "Failed to generate random structure:  atom_volume= ", \
          opts.atom_volume, " group= ", group,  "\n   atoms= ", atoms
      if not struc is None:
        opts.resolution = res
        opts.shrink_radius = shr_rad
        opts.solvent_radius = slv_rad
        compare_masks( struc, opts)

def get_resolution(pdb_input, default_resolution):
  strs = pdb_input.extract_remark_iii_records(2)
  res =  iotbx.pdb.remark_2_interpretation.extract_resolution(strs)
  if res is None:
    return default_resolution
  else:
    return res[0]

def cci_vetted_tests( options) :
  print "CCI tests,  \n options="
  print options, "\n"
  n_files = options.cci
  assert n_files > 0
  d = os.environ.get("CCI_REFINE_VETTED")
  assert not d is None, "Tests on CCI structures requested, but" \
    " CCI_REFINE_VETTED is not defined"
  assert os.path.isdir( d ), d
  resolution = options.resolution
  print "Testing files in ", d
  fls = os.listdir( d )
  n = 0
  for f in fls:
    ffull = os.path.join(d, f)
    freal = os.path.abspath(ffull)
    freal = os.path.realpath(freal)
    if os.path.isfile(freal):
      fbase =  os.path.basename(f).lower()
      if fbase.find("pdb") != -1 :
        n = n + 1
        print "Processing file: ", f
        pdb_inp = iotbx.pdb.input(source_info = None, file_name = ffull)
        struc = pdb_inp.xray_structure_simple()
        options.resolution = get_resolution( pdb_inp, resolution)
        compare_masks(struc, options)
        if( n>=n_files ):
          break
  print "Number of structures tested: ", n
  assert n>0, "No CCI files have been tested"

def generate_cb(grp, ncb):
  halls = masks.generate_groups(grp,ncb)
  for h in halls:
    print "Generated: ", h
  return halls

def run():
  import optparse
  parser = optparse.OptionParser()
  parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
      default=False, help="be verbose")
  parser.add_option("-g", "--space_group", action="store", type="string",
      dest="space_group", help="space group symbol or number")
  parser.add_option("-n", "--n_atoms", action="store", type="int",
      dest="n_atoms", default=4, help="number of atoms in the asymmetric unit")
  parser.add_option("-a", "--atom_volume", action="store", type="float",
      dest="atom_volume", default=50.0, help="volume of one atom in agstrom^3")
  parser.add_option("-z", "--random", action="store", type="int",
      dest="random", default=0, help="number of random structures per space group")
  parser.add_option("-c", "--cci", action="store", type="int",
      dest="cci", default=0, help="number of structures from CCI PDB database")
  parser.add_option("-f", "--pdb_file", action="store", type="string",
      dest="file", help="pdb file to test")
  parser.add_option("-r", "--resolution", action="store", type="float",
      dest="resolution", default=1.972, help="resolution")
  parser.add_option("-t", "--tolerance", action="store", type="float",
      dest="tolerance", default=1.0E-6, help="resolution")
  parser.add_option("--solvent_radius", action="store", type="float",
      dest="solvent_radius", default=1.1, help="solvent radius")
  parser.add_option("--shrink_radius", action="store", type="float",
      dest="shrink_radius", default=0.9, help="shrink truncation radius")
  parser.add_option("--save_failed", action="store", type="string",
      dest="failed_file", help="file name for the pdb of the failed structure")
  parser.add_option("--save_files", action="store", type="string",
      dest="save_files", help="base file name for pdb/mask/sf files to save")
  parser.add_option("--change_basis", action="store", type="int",
      dest="change_basis", default=0, help="number of basis")
  parser.add_option("--groups_file", action="store", type="string",
      dest="groups_file", help="file containing space group, one per line")

  (opts, args) = parser.parse_args()

  groups = []
  if (not ((opts.space_group is None) and (opts.groups_file is None))) \
      and (opts.random == 0):
    opts.random = 1
  if not opts.groups_file is None:
    tmp_file = open(opts.groups_file, "r")
    for line in tmp_file.readlines(): # newlines retained
      groups.append( line.strip() ) # removes whitespace in the begining and end
    tmp_file.close()
  if (opts.space_group is None) & (len(groups)==0) :
    groups.extend(SpaceGroups)
  elif opts.space_group == "all" :
    for isg in xrange(1,231):
      groups.append(str(isg))
  elif opts.space_group == "all530":
    it = cctbx.sgtbx.space_group_symbol_iterator()
    while( True ):
      symbol = it.next()
      # TODO: the following  does not work
      #if( symbol.number()==0 ):
      #  break
      groups.append(symbol.hermann_mauguin())
      if( symbol.number()==230 ):
        break
  elif not opts.space_group is None:
    groups.append(opts.space_group)

  if opts.change_basis != 0:
    cb_groups = []
    for grp in groups:
      halls = generate_cb(grp, opts.change_basis)
      for hall in halls:
        cb_groups.append( hall )
    tmp_file = open("generated_groups.txt", "w")
    for g in cb_groups:
      print >>tmp_file, g
    tmp_file.close()
    groups = cb_groups

  if opts.random > 0:
    random_tests(groups, opts)
  if opts.cci >0:
    cci_vetted_tests(opts)
  if not opts.file is None:
    pdb_inp = iotbx.pdb.input(source_info = None, file_name = opts.file)
    struc = pdb_inp.xray_structure_simple()
    opts.resolution = get_resolution(pdb_inp, opts.resolution)
    compare_masks(struc, opts)
  elif opts.cci == 0 and opts.random == 0 and opts.space_group is None:
    standard_tests(groups, opts)

  print format_cpu_times()

def exercise_mask_data_1(space_group_info, n_sites=100):
  from cctbx import maptbx
  from cctbx.masks import vdw_radii_from_xray_structure
  for d_min in [1, 1.5, 2.1]:
    for resolution_factor in [1./2, 1./3, 1./4, 1./5]:
      xrs = random_structure.xray_structure(
        space_group_info=space_group_info,
        elements=(("O","N","C")*(n_sites//3+1))[:n_sites],
        volume_per_atom=30,
        min_distance=1)
      atom_radii = vdw_radii_from_xray_structure(xray_structure = xrs)
      asu_mask = masks.atom_mask(
        unit_cell                = xrs.unit_cell(),
        group                    = xrs.space_group(),
        resolution               = d_min,
        grid_step_factor         = resolution_factor,
        solvent_radius           = 1.0,
        shrink_truncation_radius = 1.0)
      asu_mask.compute(xrs.sites_frac(), atom_radii)
      mask_data = asu_mask.mask_data_whole_uc()
      assert flex.min(mask_data) == 0.0
      # It's not just 0 and 1 ...
      assert flex.max(mask_data) == xrs.space_group().order_z()
      # In fact, it is a mixture ...
      if 0: # XXX this will rightfully crash
        mask_data_ = mask_data / xrs.space_group().order_z()
        s0 = mask_data_ < 0.5
        s1 = mask_data_ > 0.5
        if(mask_data_.size() != s0.count(True)+s1.count(True)):
          for d in mask_data_:
            if(d != 0 and d != 1): print d, xrs.space_group().order_z()
          assert mask_data_.size() == s0.count(True)+s1.count(True), [
            mask_data_.size()-(s0.count(True)+s1.count(True))]
      if(0): # XXX This would crash with the message: "... The grid is not ..."
        cr_gr = maptbx.crystal_gridding(
          unit_cell         = xrs.unit_cell(),
          d_min             = d_min,
          resolution_factor = resolution_factor)
        asu_mask = masks.atom_mask(
          unit_cell                = xrs.unit_cell(),
          space_group              = xrs.space_group(),
          gridding_n_real          = cr_gr.n_real(),
          solvent_radius           = 1.0,
          shrink_truncation_radius = 1.0)
        asu_mask.compute(xrs.sites_frac(), atom_radii)

def exercise_mask_data_2(space_group_info, n_sites=100, d_min=2.0,
                         resolution_factor=1./4):
  from cctbx import maptbx
  from cctbx.masks import vdw_radii_from_xray_structure
  for yn in [0,1]:
    xrs = random_structure.xray_structure(
      space_group_info=space_group_info,
      elements=(("O","N","C")*(n_sites//3+1))[:n_sites],
      volume_per_atom=50,
      min_distance=1.5)
    atom_radii = vdw_radii_from_xray_structure(xray_structure = xrs)
    asu_mask = masks.atom_mask(
      unit_cell                = xrs.unit_cell(),
      group                    = xrs.space_group(),
      resolution               = d_min,
      grid_step_factor         = resolution_factor,
      solvent_radius           = 1.0,
      shrink_truncation_radius = 1.0)
    asu_mask.compute(xrs.sites_frac(), atom_radii)
    mask_data = asu_mask.mask_data_whole_uc()
    mask_data = mask_data / xrs.space_group().order_z()
    if(yn == 1):
      mask_data = maptbx.copy(mask_data, flex.grid(mask_data.focus()))
    fc = xrs.structure_factors(d_min = d_min).f_calc()
    f_mask_1 = fc.set().array(data = asu_mask.structure_factors(fc.indices()))
    f_mask_2 = f_mask_1.structure_factors_from_map(map=mask_data,
      use_scale = True, anomalous_flag = False, use_sg = True)
    fm1 = abs(f_mask_1).data()
    fm2 = abs(f_mask_2).data()
    r = flex.sum( flex.abs( fm1 - fm2 ) ) / flex.sum( fm1 + fm2 )
    assert approx_equal(r, 0.0)

def run_call_back(flags, space_group_info):
  exercise_mask_data_1(space_group_info)
  exercise_mask_data_2(space_group_info)

if (__name__ == "__main__"):
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  try:
    run()
  except :
    log = cout.getvalue()
    if len(log) != 0:
      print "<<<<<<<< Start Log:"
      print cout.getvalue()
      print ">>>>>>>> End Log"
    raise
