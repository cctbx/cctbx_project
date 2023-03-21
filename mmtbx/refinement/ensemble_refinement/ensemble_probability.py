from __future__ import absolute_import, division, print_function
import sys
from iotbx.option_parser import iotbx_option_parser
from mmtbx import utils
from iotbx import file_reader
import iotbx.phil
import libtbx.phil
from libtbx.utils import Sorry
from cctbx import maptbx
from cctbx.array_family import flex
from libtbx import group_args
from cctbx import adptbx
from iotbx import pdb
from libtbx.utils import Sorry, multi_out
from six.moves import cStringIO as StringIO
from six.moves import range
from iotbx import extract_xtal_data

master_phil = iotbx.phil.parse("""\
ensemble_probability {
  verbose = True
    .type = bool
    .help = '''Verbose'''
  assign_sigma_from_map = False
    .type = bool
    .help = '''ensemble mtz file containing precalculated map coeffs'''
  ensemble_sigma_map_input = None
    .type = str
    .help = '''ensemble mtz file containing precalculated map coeffs'''
  output_model_and_model_ave_mtz = False
    .type = bool
    .help = '''Output individual and average Fcalc maps'''
  fcalc_high_resolution = 2.0
    .type = float
    .help = '''high resolution limit for Fcalc map'''
  fcalc_low_resolution = None
    .type = float
    .help = '''low resolution limit for Fcalc map'''
  residue_detail = False
    .type = bool
    .help = '''Average probability per residue'''
  ignore_hd = True
    .type = bool
    .help = '''Ignore H/D, only applicable when residue detail is used'''
  sort_ensemble_by_nll_score = False
    .type = bool
    .help = '''ordered each model in ensemble by negative log liklihood score and output'''
  fobs_vs_fcalc_post_nll = False
    .type = bool
    .help = '''Recalc R factors based removing models or atoms based on nll score'''
}
""")

class map_cc_funct(object):
  def __init__(self, map_1,
                     xray_structure,
                     fft_map,
                     atom_radius,
                     hydrogen_atom_radius,
                     model_i,
                     number_previous_scatters,
                     ignore_hd = False,
                     residue_detail = True,
                     selection = None,
                     pdb_hierarchy = None):
    self.xray_structure = xray_structure
    self.selection = selection
    self.pdb_hierarchy = pdb_hierarchy
    self.result = []
    self.map_1_size = map_1.size()
    self.map_1_stat = maptbx.statistics(map_1)
    self.atoms_with_labels = None
    self.residue_detail = residue_detail
    self.model_i = model_i
    if(pdb_hierarchy is not None):
      self.atoms_with_labels = list(pdb_hierarchy.atoms_with_labels())
    scatterers = self.xray_structure.scatterers()
    sigma_occ = flex.double()
    if(self.selection is None):
      self.selection = flex.bool(scatterers.size(), True)
    real_map_unpadded = fft_map.real_map_unpadded()
    sites_cart = self.xray_structure.sites_cart()

    if not self.residue_detail:
      self.gifes = [None,]*scatterers.size()
      self._result = [None,]*scatterers.size()
      #
      atom_radii = flex.double(scatterers.size(), atom_radius)
      for i_seq, sc in enumerate(scatterers):
        if(self.selection[i_seq]):
          if(sc.element_symbol().strip().lower() in ["h","d"]):
            atom_radii[i_seq] = hydrogen_atom_radius
      #
      for i_seq, site_cart in enumerate(sites_cart):
        if(self.selection[i_seq]):
          sel = maptbx.grid_indices_around_sites(
            unit_cell  = self.xray_structure.unit_cell(),
            fft_n_real = real_map_unpadded.focus(),
            fft_m_real = real_map_unpadded.all(),
            sites_cart = flex.vec3_double([site_cart]),
            site_radii = flex.double([atom_radii[i_seq]]))
          self.gifes[i_seq] = sel
          m1 = map_1.select(sel)
          ed1 = map_1.eight_point_interpolation(scatterers[i_seq].site)
          sigma_occ.append(ed1)
          a = None
          if(self.atoms_with_labels is not None):
            a = self.atoms_with_labels[i_seq]
          self._result[i_seq] = group_args(atom = a, m1 = m1, ed1 = ed1,
            xyz=site_cart)
      self.xray_structure.set_occupancies(sigma_occ)

      ### For testing other residue averaging options
      residues = self.extract_residues(model_i = model_i,
                                       number_previous_scatters = number_previous_scatters)
      self.xray_structure.residue_selections = residues

    # Residue detail
    if self.residue_detail:
      assert self.pdb_hierarchy is not None
      residues = self.extract_residues(model_i = model_i,
                                       number_previous_scatters = number_previous_scatters)
      self.gifes = [None,]*len(residues)
      self._result = [None,]*len(residues)
      for i_seq, residue in enumerate(residues):
        residue_sites_cart = sites_cart.select(residue.selection)
        if 0: print(i_seq, list(residue.selection)) # DEBUG
        sel = maptbx.grid_indices_around_sites(
          unit_cell  = self.xray_structure.unit_cell(),
          fft_n_real = real_map_unpadded.focus(),
          fft_m_real = real_map_unpadded.all(),
          sites_cart = residue_sites_cart,
          site_radii = flex.double(residue.selection.size(), atom_radius))
        self.gifes[i_seq] = sel
        m1 = map_1.select(sel)
        ed1 = flex.double()
        for i_seq_r in residue.selection:
          ed1.append(map_1.eight_point_interpolation(scatterers[i_seq_r].site))
        self._result[i_seq] = \
          group_args(residue = residue, m1 = m1, ed1 = flex.mean(ed1),
            xyz=residue_sites_cart.mean(), n_atoms=residue_sites_cart.size())

        residue_scatterers = scatterers.select(residue.selection)
        residue_ed1 = flex.double()
        for n,scatter in enumerate(residue_scatterers):
          if ignore_hd:
            if scatter.element_symbol() not in ['H', 'D']:
              residue_ed1.append(ed1[n])
          else:
            residue_ed1.append(ed1[n])
        for x in range(ed1.size()):
          sigma_occ.append(flex.mean(residue_ed1))

      self.xray_structure.set_occupancies(sigma_occ)
      self.xray_structure.residue_selections = residues

    del map_1

  def extract_residues(self, model_i, number_previous_scatters, combine = True):
    result = []
    model = self.pdb_hierarchy.models()[model_i]
    rm = []
    for chain in model.chains():
      for rg in chain.residue_groups():
        rg_i_seqs = []
        r_name = None
        for ag in rg.atom_groups():
          if(r_name is None): r_name = ag.resname
          for atom in ag.atoms():
            if(self.selection[atom.i_seq - number_previous_scatters]):
              rg_i_seqs.append(atom.i_seq - number_previous_scatters)
        if(len(rg_i_seqs) != 0):
          rm.append(group_args(
            selection = flex.size_t(rg_i_seqs),
            name      = r_name,
            model_id  = model_i,
            resid     = rg.resid(),
            chain_id  = chain.id))
    result.append(rm)

    if(combine):
      r0 = result[0]
      for r in result[1:]:
        for i, ri in enumerate(r):
          r0[i].selection.extend(ri.selection)
          assert r0[i].name == ri.name
    else:
      r0 = result[0]
      for r in result[1:]:
        r0.extend(r)

    return r0


def get_map_sigma(ens_pdb_hierarchy,
                           ens_pdb_xrs,
                           log,
                           model_i,
                           number_previous_scatters,
                           residue_detail = True,
                           ignore_hd = True,
                           map_coeffs_1 = None,
                           fft_map_1 = None):
  assert [map_coeffs_1, fft_map_1].count(None) == 1
  if fft_map_1 == None:
    fft_map_1 = map_coeffs_1.fft_map()
    fft_map_1.apply_sigma_scaling()
  map_1 = fft_map_1.real_map_unpadded()
  atom_radius = 1.5
  hydrogen_atom_radius = 1.0
  map_cc_obj = map_cc_funct(
    map_1                = map_1,
    xray_structure       = ens_pdb_xrs,
    model_i              = model_i,
    number_previous_scatters = number_previous_scatters,
    fft_map              = fft_map_1,
    atom_radius          = atom_radius,
    hydrogen_atom_radius = hydrogen_atom_radius,
    selection            = None,
    residue_detail       = residue_detail,
    ignore_hd            = ignore_hd,
    pdb_hierarchy        = ens_pdb_hierarchy)
  return map_cc_obj.xray_structure

def write_ensemble_pdb(filename,
                       xrs_list,
                       ens_pdb_hierarchy):
    out = open(filename, 'w')
    crystal_symmetry = xrs_list[0].crystal_symmetry()
    print("REMARK   3  TIME-AVERAGED ENSEMBLE REFINEMENT", file=out)
    print("REMARK   3  OCCUPANCY = MAP SIGMA LEVEL", file=out)
    print(pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry), file=out)
    print(pdb.format_scale_records(unit_cell = crystal_symmetry.unit_cell()), file=out)
    atoms_reset_serial = True

    for i_model, xrs in enumerate(xrs_list):
      scatterers = xrs.scatterers()
      sites_cart = xrs.sites_cart()
      u_isos = xrs.extract_u_iso_or_u_equiv()
      occupancies = scatterers.extract_occupancies()
      u_carts = scatterers.extract_u_cart_plus_u_iso(xrs.unit_cell())
      scat_types = scatterers.extract_scattering_types()
      i_model_ens_pdb_hierarchy = ens_pdb_hierarchy.models()[i_model]
      pdb_atoms = i_model_ens_pdb_hierarchy.atoms()
      for j_seq, atom in enumerate(pdb_atoms):
        if j_seq < len(sites_cart):
          atom.xyz = sites_cart[j_seq]
          atom.occ = occupancies[j_seq]
          atom.b = adptbx.u_as_b(u_isos[j_seq])
          e = scat_types[j_seq]
          if (len(e) > 1 and "+-0123456789".find(e[1]) >= 0):
            atom.element = "%2s" % e[:1]
            atom.charge = "%-2s" % e[1:]
          elif (len(e) > 2):
            atom.element = "%2s" % e[:2]
            atom.charge = "%-2s" % e[2:]
          else:
            atom.element = "%2s" % e
            atom.charge = "  "
          if (scatterers[j_seq].flags.use_u_aniso()):
            atom.uij = u_carts[j_seq]
          elif(False):
            atom.uij = self.u_cart
          else:
            atom.uij = (-1,-1,-1,-1,-1,-1)
      if (atoms_reset_serial):
        atoms_reset_serial_first_value = 1
      else:
        atoms_reset_serial_first_value = None
    out.write(ens_pdb_hierarchy.as_pdb_string(
      append_end=False,
      atoms_reset_serial_first_value=atoms_reset_serial_first_value))

def reflection_file_server(crystal_symmetry, reflection_files, log):
  from iotbx import reflection_file_utils
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=log)

class ensemble_probability(object):
  def run(self, args, command_name, out=sys.stdout):
    command_line = (iotbx_option_parser(
      usage="%s [options]" % command_name,
      description='Example: %s data.mtz data.mtz ref_model.pdb'%command_name)
      .option(None, "--show_defaults",
        action="store_true",
        help="Show list of parameters.")
      ).process(args=args)

    cif_file = None
    processed_args = utils.process_command_line_args(
                       args          = args,
                       log           = sys.stdout,
                       master_params = master_phil)
    params = processed_args.params
    if(params is None): params = master_phil
    self.params = params.extract().ensemble_probability
    pdb_file_names = processed_args.pdb_file_names
    if len(pdb_file_names) != 1 :
      raise Sorry("Only one PDB structure may be used")
    pdb_file = file_reader.any_file(pdb_file_names[0])
    self.log = multi_out()
    self.log.register(label="stdout", file_object=sys.stdout)
    self.log.register(
      label="log_buffer",
      file_object=StringIO(),
      atexit_send_to=None)
    sys.stderr = self.log
    log_file = open(pdb_file_names[0].split('/')[-1].replace('.pdb','') + '_pensemble.log', "w")

    self.log.replace_stringio(
        old_label="log_buffer",
        new_label="log",
        new_file_object=log_file)
    utils.print_header(command_name, out = self.log)
    params.show(out = self.log)
    #
    f_obs = None
    r_free_flags = None
    reflection_files = processed_args.reflection_files

    if self.params.fobs_vs_fcalc_post_nll:
      if len(reflection_files) == 0:
        raise Sorry("Fobs from input MTZ required for fobs_vs_fcalc_post_nll")

    if len(reflection_files) > 0:
      crystal_symmetry = processed_args.crystal_symmetry
      print('Reflection file : ', processed_args.reflection_file_names[0], file=self.log)
      utils.print_header("Model and data statistics", out = self.log)
      rfs = reflection_file_server(
        crystal_symmetry = crystal_symmetry,
        reflection_files = processed_args.reflection_files,
        log              = self.log)

      parameters = extract_xtal_data.data_and_flags_master_params().extract()
      determine_data_and_flags_result = extract_xtal_data.run(
        reflection_file_server = rfs,
        parameters             = parameters,
        keep_going             = True)
      f_obs = determine_data_and_flags_result.f_obs
      number_of_reflections = f_obs.indices().size()
      r_free_flags = determine_data_and_flags_result.r_free_flags
      test_flag_value = determine_data_and_flags_result.test_flag_value
      if(r_free_flags is None):
        r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))

    # process PDB
    pdb_file.assert_file_type("pdb")
    #
    pdb_in = pdb.input(file_name=pdb_file.file_name)
    ens_pdb_hierarchy = pdb_in.construct_hierarchy()
    ens_pdb_hierarchy.atoms().reset_i_seq()
    ens_pdb_xrs_s = pdb_in.xray_structures_simple()
    number_structures = len(ens_pdb_xrs_s)
    print('Number of structure in ensemble : ', number_structures, file=self.log)

    # Calculate sigmas from input map only
    if self.params.assign_sigma_from_map and self.params.ensemble_sigma_map_input is not None:
      # process MTZ
      input_file = file_reader.any_file(self.params.ensemble_sigma_map_input)
      if input_file.file_type == "hkl" :
        if input_file.file_object.file_type() != "ccp4_mtz" :
           raise Sorry("Only MTZ format accepted for map input")
        else:
          mtz_file = input_file
      else:
        raise Sorry("Only MTZ format accepted for map input")
      miller_arrays = mtz_file.file_server.miller_arrays
      map_coeffs_1 = miller_arrays[0]
      #
      xrs_list = []
      for n, ens_pdb_xrs in enumerate(ens_pdb_xrs_s):
        # get sigma levels from ensemble fc for each structure
        xrs = get_map_sigma(ens_pdb_hierarchy = ens_pdb_hierarchy,
                          ens_pdb_xrs       = ens_pdb_xrs,
                          map_coeffs_1      = map_coeffs_1,
                          residue_detail    = self.params.residue_detail,
                          ignore_hd         = self.params.ignore_hd,
                          log               = self.log)
        xrs_list.append(xrs)
      # write ensemble pdb file, occupancies as sigma level
      filename = pdb_file_names[0].split('/')[-1].replace('.pdb','') + '_vs_' + self.params.ensemble_sigma_map_input.replace('.mtz','') + '_pensemble.pdb'
      write_ensemble_pdb(filename = filename,
                         xrs_list = xrs_list,
                         ens_pdb_hierarchy = ens_pdb_hierarchy
                         )

    # Do full analysis vs Fobs
    else:
      model_map_coeffs = []
      fmodel = None
      # Get <fcalc>
      for model, ens_pdb_xrs in enumerate(ens_pdb_xrs_s):
        ens_pdb_xrs.set_occupancies(1.0)
        if model == 0:
          # If mtz not supplied get fobs from xray structure...
          # Use input Fobs for scoring against nll
          if self.params.fobs_vs_fcalc_post_nll:
            dummy_fobs = f_obs
          else:
            if f_obs == None:
              if self.params.fcalc_high_resolution == None:
                raise Sorry("Please supply high resolution limit or input mtz file.")
              dummy_dmin = self.params.fcalc_high_resolution
              dummy_dmax = self.params.fcalc_low_resolution
            else:
              print('Supplied mtz used to determine high and low resolution cuttoffs', file=self.log)
              dummy_dmax, dummy_dmin = f_obs.d_max_min()
            #
            dummy_fobs = abs(ens_pdb_xrs.structure_factors(d_min = dummy_dmin).f_calc())
            dummy_fobs.set_observation_type_xray_amplitude()
            # If mtz supplied, free flags are over written to prevent array size error
            r_free_flags = dummy_fobs.array(data=flex.bool(dummy_fobs.data().size(),False))
          #
          fmodel = utils.fmodel_simple(
                     scattering_table         = "wk1995",
                     xray_structures          = [ens_pdb_xrs],
                     f_obs                    = dummy_fobs,
                     target_name              = 'ls',
                     bulk_solvent_and_scaling = False,
                     r_free_flags             = r_free_flags
                     )
          f_calc_ave = fmodel.f_calc().array(data = fmodel.f_calc().data()*0).deep_copy()
          # XXX Important to ensure scale is identical for each model and <model>
          fmodel.set_scale_switch = 1.0
          f_calc_ave_total = fmodel.f_calc().data().deep_copy()
        else:
          fmodel.update_xray_structure(xray_structure  = ens_pdb_xrs,
                                       update_f_calc   = True,
                                       update_f_mask   = False)
          f_calc_ave_total += fmodel.f_calc().data().deep_copy()
        print('Model :', model+1, file=self.log)
        print("\nStructure vs real Fobs (no bulk solvent or scaling)", file=self.log)
        print('Rwork          : %5.4f '%fmodel.r_work(), file=self.log)
        print('Rfree          : %5.4f '%fmodel.r_free(), file=self.log)
        print('K1             : %5.4f '%fmodel.scale_k1(), file=self.log)
        fcalc_edm        = fmodel.electron_density_map()
        fcalc_map_coeffs = fcalc_edm.map_coefficients(map_type = 'Fc')
        fcalc_mtz_dataset = fcalc_map_coeffs.as_mtz_dataset(column_root_label ='Fc')
        if self.params.output_model_and_model_ave_mtz:
          fcalc_mtz_dataset.mtz_object().write(file_name = str(model+1)+"_Fc.mtz")
        model_map_coeffs.append(fcalc_map_coeffs.deep_copy())

      fmodel.update(f_calc = f_calc_ave.array(f_calc_ave_total / number_structures))
      print("\nEnsemble vs real Fobs (no bulk solvent or scaling)", file=self.log)
      print('Rwork          : %5.4f '%fmodel.r_work(), file=self.log)
      print('Rfree          : %5.4f '%fmodel.r_free(), file=self.log)
      print('K1             : %5.4f '%fmodel.scale_k1(), file=self.log)

      # Get <Fcalc> map
      fcalc_ave_edm        = fmodel.electron_density_map()
      fcalc_ave_map_coeffs = fcalc_ave_edm.map_coefficients(map_type = 'Fc').deep_copy()
      fcalc_ave_mtz_dataset = fcalc_ave_map_coeffs.as_mtz_dataset(column_root_label ='Fc')
      if self.params.output_model_and_model_ave_mtz:
        fcalc_ave_mtz_dataset.mtz_object().write(file_name = "aveFc.mtz")
      fcalc_ave_map_coeffs = fcalc_ave_map_coeffs.fft_map()
      fcalc_ave_map_coeffs.apply_volume_scaling()
      fcalc_ave_map_data   = fcalc_ave_map_coeffs.real_map_unpadded()
      fcalc_ave_map_stats  = maptbx.statistics(fcalc_ave_map_data)

      print("<Fcalc> Map Stats :", file=self.log)
      fcalc_ave_map_stats.show_summary(f = self.log)
      offset = fcalc_ave_map_stats.min()
      model_neg_ll = []

      number_previous_scatters = 0

      # Run through structure list again and get probability
      xrs_list = []
      for model, ens_pdb_xrs in enumerate(ens_pdb_xrs_s):
        if self.params.verbose:
          print('\n\nModel                   : ', model+1, file=self.log)
        # Get model atom sigmas vs Fcalc
        fcalc_map = model_map_coeffs[model].fft_map()
        fcalc_map.apply_volume_scaling()
        fcalc_map_data  = fcalc_map.real_map_unpadded()
        fcalc_map_stats  = maptbx.statistics(fcalc_map_data)
        if self.params.verbose:
          print("Fcalc map stats         :", file=self.log)
        fcalc_map_stats.show_summary(f = self.log)

        xrs = get_map_sigma(ens_pdb_hierarchy = ens_pdb_hierarchy,
                            ens_pdb_xrs       = ens_pdb_xrs,
                            fft_map_1         = fcalc_map,
                            model_i           = model,
                            residue_detail    = self.params.residue_detail,
                            ignore_hd         = self.params.ignore_hd,
                            number_previous_scatters = number_previous_scatters,
                            log               = self.log)
        fcalc_sigmas = xrs.scatterers().extract_occupancies()
        del fcalc_map
        # Get model atom sigmas vs <Fcalc>
        xrs = get_map_sigma(ens_pdb_hierarchy = ens_pdb_hierarchy,
                            ens_pdb_xrs       = ens_pdb_xrs,
                            fft_map_1         = fcalc_ave_map_coeffs,
                            model_i           = model,
                            residue_detail    = self.params.residue_detail,
                            ignore_hd         = self.params.ignore_hd,
                            number_previous_scatters = number_previous_scatters,
                            log               = self.log)

        ### For testing other residue averaging options
        #print xrs.residue_selections

        fcalc_ave_sigmas = xrs.scatterers().extract_occupancies()
        # Probability of model given <model>
        prob = fcalc_ave_sigmas / fcalc_sigmas
        # XXX debug option
        if False:
          for n,p in enumerate(prob):
            print(' {0:5d} {1:5.3f}'.format(n,p), file=self.log)
        # Set probabilty between 0 and 1
        # XXX Make Histogram / more stats
        prob_lss_zero = flex.bool(prob <= 0)
        prob_grt_one = flex.bool(prob > 1)
        prob.set_selected(prob_lss_zero, 0.001)
        prob.set_selected(prob_grt_one, 1.0)
        xrs.set_occupancies(prob)
        xrs_list.append(xrs)
        sum_neg_ll = sum(-flex.log(prob))
        model_neg_ll.append((sum_neg_ll, model))
        if self.params.verbose:
          print('Model probability stats :', file=self.log)
          print(prob.min_max_mean().show(), file=self.log)
          print('  Count < 0.0 : ', prob_lss_zero.count(True), file=self.log)
          print('  Count > 1.0 : ', prob_grt_one.count(True), file=self.log)

        # For averaging by residue
        number_previous_scatters += ens_pdb_xrs.sites_cart().size()

      # write ensemble pdb file, occupancies as sigma level
      write_ensemble_pdb(filename = pdb_file_names[0].split('/')[-1].replace('.pdb','') + '_pensemble.pdb',
                       xrs_list = xrs_list,
                       ens_pdb_hierarchy = ens_pdb_hierarchy
                       )

      # XXX Test ordering models by nll
      # XXX Test removing nth percentile atoms
      if self.params.sort_ensemble_by_nll_score or self.params.fobs_vs_fcalc_post_nll:
        for percentile in [1.0,0.975,0.95,0.9,0.8,0.6,0.2]:
          model_neg_ll = sorted(model_neg_ll)
          f_calc_ave_total_reordered = None
          print_list = []
          for i_neg_ll in model_neg_ll:
            xrs = xrs_list[i_neg_ll[1]]
            nll_occ = xrs.scatterers().extract_occupancies()

            # Set q=0 nth percentile atoms
            sorted_nll_occ = sorted(nll_occ, reverse=True)
            number_atoms = len(sorted_nll_occ)
            percentile_prob_cutoff = sorted_nll_occ[int(number_atoms * percentile)-1]
            cutoff_selections = flex.bool(nll_occ < percentile_prob_cutoff)
            cutoff_nll_occ = flex.double(nll_occ.size(), 1.0).set_selected(cutoff_selections, 0.0)
            #XXX Debug
            if False:
              print('\nDebug')
              for x in range(len(cutoff_selections)):
                print(cutoff_selections[x], nll_occ[x], cutoff_nll_occ[x])
              print(percentile)
              print(percentile_prob_cutoff)
              print(cutoff_selections.count(True))
              print(cutoff_selections.size())
              print(cutoff_nll_occ.count(0.0))
              print('Count q = 1           : ', cutoff_nll_occ.count(1.0))
              print('Count scatterers size : ', cutoff_nll_occ.size())

            xrs.set_occupancies(cutoff_nll_occ)
            fmodel.update_xray_structure(xray_structure  = xrs,
                                         update_f_calc   = True,
                                         update_f_mask   = True)

            if f_calc_ave_total_reordered == None:
              f_calc_ave_total_reordered = fmodel.f_calc().data().deep_copy()
              f_mask_ave_total_reordered = fmodel.f_masks()[0].data().deep_copy()
              cntr = 1
            else:
              f_calc_ave_total_reordered += fmodel.f_calc().data().deep_copy()
              f_mask_ave_total_reordered += fmodel.f_masks()[0].data().deep_copy()
              cntr+=1
            fmodel.update(f_calc = f_calc_ave.array(f_calc_ave_total_reordered / cntr).deep_copy(),
                          f_mask = f_calc_ave.array(f_mask_ave_total_reordered / cntr).deep_copy()
                          )

            # Update solvent and scale
            # XXX Will need to apply_back_trace on latest version
            fmodel.set_scale_switch = 0
            fmodel.update_all_scales()

            # Reset occ for outout
            xrs.set_occupancies(nll_occ)
            # k1 updated vs Fobs
            if self.params.fobs_vs_fcalc_post_nll:
              print_list.append([cntr, i_neg_ll[0], i_neg_ll[1], fmodel.r_work(), fmodel.r_free()])

          # Order models by nll and print summary
          print('\nModels ranked by nll <Fcalc> R-factors recalculated', file=self.log)
          print('Percentile cutoff : {0:5.3f}'.format(percentile), file=self.log)
          xrs_list_sorted_nll = []
          print('      |      NLL     <Rw>     <Rf>    Ens Model', file=self.log)
          for info in print_list:
            print(' {0:4d} | {1:8.1f} {2:8.4f} {3:8.4f} {4:12d}'.format(
              info[0],
              info[1],
              info[3],
              info[4],
              info[2]+1,
              ), file=self.log)
            xrs_list_sorted_nll.append(xrs_list[info[2]])

        # Output nll ordered ensemble

        write_ensemble_pdb(filename = 'nll_ordered_' + pdb_file_names[0].split('/')[-1].replace('.pdb','') + '_pensemble.pdb',
                       xrs_list = xrs_list_sorted_nll,
                       ens_pdb_hierarchy = ens_pdb_hierarchy
                       )


if __name__ == "__main__":
  ep = ensemble_probability()
  ep.run(args         = sys.argv[1:],
         command_name = 'ensemble_probability')
