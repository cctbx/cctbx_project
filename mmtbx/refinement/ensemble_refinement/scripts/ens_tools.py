from __future__ import absolute_import, division, print_function
# Ensemble tools for pymol
# Tom Burnley

import math, sys
from pymol import cmd, stored
from six.moves import zip
from six.moves import range

class LogWriter:
      def __init__(self, stdout, filename):
          self.stdout = stdout
          self.logfile = open(filename, 'a')

      def write(self, text):
          self.stdout.write(text)
          self.logfile.write(text)

      def close(self):
          self.stdout.close()
          self.logfile.close()

def print_array_stats(array = None,
                      log = None):
    if array is not None:
        if len(array) == 0: array = [0.0]
        n         = len(array)
        mean      = sum(array) / len(array)
        maximum   = max(array)
        minimum   = min(array)
        print("n         : %4d" % n, file=log)
        print("mean      : %4.3f" % mean, file=log)
        print("min       : %4.3f" % minimum, file=log)
        print("max       : %4.3f" % maximum, file=log)

# Return distance between two coords
def distance(x,y):
    return math.sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]))

def ens_measure(pk1 = None,
                pk2 = None,
                pk3 = None,
                pk4 = None,
                name = None,
                log = None,
                verbose = True):
    '''
DESCRIPTION

    Statistics from ensemble structure measurements.  If:
      2 selections give = distance
      3 selections give = angle
      4 selections give = dihedral angle

USAGE

    ens_measure pk1, pk2, pk3, pk4, name, log, verbose

ARGUMENTS

    log = name of log file
    verbose = prints individual measurements

 EXAMPLE

    ens_measure atom1, atom2, name = 'measure', log 'ens.log'
  '''
    print('\nEnsemble measurement', file=log)

    if [pk1, pk2, pk3, pk4].count(None) > 2:
        print('\nERROR: Please supply at least 2 seletions')
        return
    number_models = cmd.count_states(pk1)
    measurements = []

    # distance
    if [pk1, pk2, pk3, pk4].count(None) == 2:
        print('Distance', file=log)
        if name == None: name = 'ens_distance'
        # display as object
        cmd.distance(name = name,
                     selection1 = pk1,
                     selection2 = pk2)

        # get individual values
        for n in range(number_models):
            measurements.append( cmd.get_distance(pk1, pk2, n+1) )
        assert len(measurements) == number_models

    # angle
    if [pk1, pk2, pk3, pk4].count(None) == 1:
        print('Angle', file=log)
        # display as object
        if name == None: name = 'ens_angle'
        cmd.angle(name = name,
                  selection1 = pk1,
                  selection2 = pk2,
                  selection3 = pk3)

        # get individual values
        for n in range(number_models):
            measurements.append( cmd.get_angle(atom1 = pk1,
                                               atom2 = pk2,
                                               atom3 = pk3,
                                               state = n+1) )
        assert len(measurements) == number_models

    # Dihedral angle
    if [pk1, pk2, pk3, pk4].count(None) == 0:
        print('Dihedral angle', file=log)
        # display as object
        if name == None: name = 'ens_dihedral'
        cmd.dihedral(name = name,
                  selection1 = pk1,
                  selection2 = pk2,
                  selection3 = pk3,
                  selection4 = pk4)

        # get individual values
        for n in range(number_models):
            measurements.append( cmd.get_dihedral(atom1 = pk1,
                                                  atom2 = pk2,
                                                  atom3 = pk3,
                                                  atom4 = pk4,
                                                  state = n+1) )
        assert len(measurements) == number_models

    # print stats
    if verbose:
        print(' State  Value', file=log)
        for n, measurement in enumerate(measurements):
            print('  %4d  %3.3f '%(n+1, measurement), file=log)

    print('\nMeasurement statistics', file=log)
    print_array_stats(array                 = measurements,
                      log                   = log)

def ens_rmsd(ens_selection,
             ref_selection,
             log_name = None):
    '''
DESCRIPTION

    Prints RMSD per structure in ensemble w.r.t. a reference structure

USAGE

    ens_rmsd ensemble_selection, reference_selection, name, log,

ARGUMENTS

    log = name of log file
    verbose = calculates structure by structure RMSD for ensemble w.r.t. a single reference structure

 EXAMPLE

    ens_rmsd ensemble_selection, reference_selection, name = 'rmsd', log 'ens.log'
  '''

    if log_name == None:
      log = LogWriter(sys.stdout, 'log.txt')
    else:
      log = LogWriter(sys.stdout, log_name+'.txt')

    # initialize arrays
    ens_selection = ens_selection + ' and not resn hoh'
    ref_selection = ref_selection + ' and not resn hoh'
    rmsd_states = []
    mean_coords = None
    number_models = cmd.count_states(ens_selection)

    # get models, mean coords
    print('\nRMSD by state', file=log)
    print('\n State | RMSD', file=log)
    for i in range(number_models):
      ens_coords = cmd.get_model(ens_selection,state=i+1).get_coord_list()
      ref_coords = cmd.get_model(ref_selection,state=1).get_coord_list()
      atom_sqr_dev = []
      for atom in range(len(ref_coords)):
        x = ref_coords[atom]
        y = ens_coords[atom]
        atom_sqr_dev.append(distance(x,y)**2)
      rmsd = math.sqrt(sum(atom_sqr_dev) / len(atom_sqr_dev))
      rmsd_states.append(rmsd)
      print(' %5d | %5.3f '%(i+1, rmsd), file=log)


    print_array_stats(array                 = rmsd_states,
                      log                   = log)
    print('\nRMSD all states : %5.3f '%(cmd.rms(ens_selection, ref_selection)))

def ens_rmsf(selection,
             rmsf_spectrum = False,
             mean_structure = False,
             mean_per_resi = True,
             log_name = None):
    '''
DESCRIPTION

    Generates and colors structure by RMSF statistics, calculated from mulistate ensembles

USAGE

    ens_rmsf selection, sigma_rmsf_spectrum, mean_structure, histogram_number_bins, log_name

ARGUMENTS

    sigma_rmsf_spectrum = overwrite q col with sigma RMSF (per atom), color by q
    mean_structure = generate new single state structure with mean coords
    histogram_number_bins = number of bins for output histograms
    log_name = name of log file

EXAMPLE

    ens_rmsf protein, True, True, 25, 'ens_rmsf.log'
  '''
    if log_name == None:
      log = LogWriter(sys.stdout, 'log.txt')
    else:
      log = LogWriter(sys.stdout, log_name+'.txt')
    print("\nEnsemble RMSF", file=log)
    print("N.B. waters excluded", file=log)
    print("N.B. B column information from first model", file=log)
    print("Rmsf (Angstrom) [w.r.t mean structure]", file=log)
    print("B_atom (Angstrom^2) [atomic Bfactor used in simulation]", file=log)
    print("B_rmsf (Angstrom^2) [rmsf converted to Bfactor]", file=log)

    # initialize arrays
    selection = selection + ' and not resn hoh'
    models = []
    mean_coords = None
    number_models = cmd.count_states(selection)
    r_number_models = 1.0 / number_models

    # get models, mean coords
    for i in range(number_models):
      models.append(cmd.get_model(selection,state=i+1))
      coords_for_mean = [[x[0]*r_number_models,x[1]*r_number_models,x[2]*r_number_models] for x in models[i].get_coord_list()]
      if mean_coords == None:
        mean_coords = coords_for_mean
      else:
        n = []
        for mean, for_mean in zip(mean_coords, coords_for_mean):
          mean = [sum(_x) for _x in zip(mean,for_mean)]
          n.append(mean)
        mean_coords = n

    # calculate RMSF w.r.t. mean structure
    rmsf_coord = [0.0]*len(mean_coords)
    for i in range(number_models):
      coord_array = models[i].get_coord_list()
      for i_seq, xyz in enumerate(coord_array):
        rmsf_coord[i_seq] += distance(xyz, mean_coords[i_seq])**2
    rmsf_coord = [(x / number_models)**0.5 for x in rmsf_coord]

    # Generate new model object with average xyz coord
    if mean_structure:
      mean_structure_name = 'mean_xyz_'+selection
      cmd.create(mean_structure_name, selection, 1)
      stored.xyz = [[v[0],v[1],v[2]] for v in mean_coords]
      cmd.alter_state(1, mean_structure_name,'(x,y,z) = stored.xyz.pop(0)')

    # convert to B factor (A^2)
    rmsf_as_b_coord = [x**2 * ((8.0 * math.pi**2) / 3.0) for x in rmsf_coord]

    # get atomic b factor info
    atom_b = [at.b for at in models[0].atom]
    b_atom_plus_b_rsmf = [sum(_x) for _x in zip(atom_b,rmsf_as_b_coord)]

    # Colour by rmsf
    if rmsf_spectrum:
#      cmd.color('grey', 'all')
#      cmd.alter('all','q=0')
      atom = models[0].atom

      for n, atom in enumerate(models[0].atom):
        atom_sel = 'id ' + str(atom.id)
        atom_action = 'q = ' + str(rmsf_sigma[n])
        cmd.alter(atom_sel,atom_action)

      print("\n\nQ infomation updated with RMSF sigma\n")
      cmd.spectrum('q',selection = selection)

      if mean_structure:
        cmd.spectrum('q', selection = mean_structure_name)

    else:
      for n, atom in enumerate(models[0].atom):
        atom_sel = 'id ' + str(atom.id)
        atom_action = 'q = ' + str(rmsf_coord[n])
        cmd.alter(atom_sel,atom_action)

      print("\n\nQ infomation updated with RMSF (Angstrom)\n")
      cmd.spectrum('q',selection = selection)

    # array stats
    print('\nB_atom (A^2): ', file=log)
    print_array_stats(array                 = atom_b,
                      log                   = log)
    print('\nRmsf (A): ', file=log)
    print_array_stats(array                 = rmsf_coord,
                      log                   = log)
    print('\nB_rmsf (A^2): ', file=log)
    print_array_stats(array                 = rmsf_as_b_coord,
                      log                   = log)
    print('\nB_atom + B_rmsf (A^2):', file=log)
    print_array_stats(array                 = b_atom_plus_b_rsmf,
                      log                   = log)

    # individual atom stats
    print('\n\n Resi | Name | Chain   | Rmsf | B_rmsf | B_atom |    B_rmsf+B_atom\n', file=log)
    for i_seq, b_factor in enumerate(atom_b):
      print(' %7s %7s %7s | %8.3f | %8.3f %8.3f | %8.3f'%(
                    models[0].atom[i_seq].resi,
                    models[0].atom[i_seq].name,
                    models[0].atom[i_seq].chain,
                    rmsf_coord[i_seq],
                    rmsf_as_b_coord[i_seq],
                    b_factor,
                    b_atom_plus_b_rsmf[i_seq]), file=log)

    if mean_per_resi:
      print('\nMean per residue', file=log)
      # update atom to include info
      for n, atom in enumerate(models[0].atom):
        atom.rmsf               = rmsf_coord[n]
        atom.b_rmsf             = rmsf_as_b_coord[n]
        atom.b_atom_plus_b_rsmf = atom.b + atom.b_rmsf

      print(' Resi Resn | Atoms | Atom_rmsf | B_atom B_rmsf B_atom+B_rmsf', file=log)
      def print_mean_residue():
          print(' %4s %5s|   %3d |  %8.3f | %8.3f %8.3f %8.3f '%(
                current_resi,
                current_resn,
                len(res_b),
                sum(res_rmsf) / len(res_rmsf),
                sum(res_b) / len(res_b),
                sum(res_b_rmsf) / len(res_b_rmsf),
                sum(res_b_atom_plus_b_rsmf) / len(res_b_atom_plus_b_rsmf) ), file=log)

      current_resi = models[0].atom[0].resi
      current_resn = models[0].atom[0].resn
      res_rmsf = []
      res_b = []
      res_b_rmsf = []
      res_b_atom_plus_b_rsmf = []

      for atom in models[0].atom:
        if atom.resi == current_resi:
          assert atom.resn == current_resn
          res_rmsf.append(atom.rmsf)
          res_b.append(atom.b)
          res_b_rmsf.append(atom.b_rmsf)
          res_b_atom_plus_b_rsmf.append(atom.b_atom_plus_b_rsmf)
        else:
          print_mean_residue()
          current_resi = atom.resi
          current_resn = atom.resn
          res_rmsf = [atom.rmsf]
          res_b = [atom.b]
          res_b_rmsf = [atom.b_rmsf]
          res_b_atom_plus_b_rsmf = [atom.b_atom_plus_b_rsmf]
      print_mean_residue()

def ens_prob():
  print('\n\nEnsemble probability options')
  # get models, mean coords
  models = []
  selection = 'all'
  for i in range(cmd.count_states(selection)):
    models.append(cmd.get_model(selection,state=i+1))

  for n, model in enumerate(models):
    residues = model.get_residues()
    for residue in residues:
      # Get individual atom info
      q_list = []
      q_list_mc  = []
      q_list_sc = []
      for i_seq in range (residue[0], residue[1]):
        # Ignore hydrogens
        if model.atom[i_seq].symbol != 'H':
          q_list.append(float(model.atom[i_seq].q))
          if model.atom[i_seq].name in ['N','CA','C','O']:
            q_list_mc.append(float(model.atom[i_seq].q))
          else:
            q_list_sc.append(float(model.atom[i_seq].q))

      # Set probability per residue
      # Mean p
      if len(q_list) > 0:    p_new = sum(q_list) / len(q_list)
      if len(q_list_mc) > 0: p_new_mc = sum(q_list_mc) / len(q_list_mc)
      if len(q_list_sc) > 0: p_new_sc = sum(q_list_sc) / len(q_list_sc)

#      # Joint
      p_new = q_list[0]
      for p in q_list[1:]:
        p_new *= p

#      # nll
#      p_new = math.log(q_list[0])
#      for p in q_list[1:]:
#        p_new += math.log(max(p,0.001))
#      p_new *= -1

      if i_seq == residue[1]-1:
        for i_seq in range (residue[0], residue[1]):
          if True:
            atom_sel = 'id ' + str(model.atom[i_seq].id) + ' and state ' + str(n+1)
            atom_action = 'b = ' + str(p_new)
            cmd.alter(atom_sel, atom_action)
          else:
            atom_sel = 'id ' + str(model.atom[i_seq].id) + ' and state ' + str(n+1)
            if model.atom[i_seq].name in ['N','CA','C','O']:
              atom_action = 'b = ' + str(p_new_mc)
            else:
              atom_action = 'b = ' + str(p_new_sc)
            cmd.alter(atom_sel, atom_action)

def print_names(selection):
  print('\n\nSelection:\n\n\n\n\n')
  selection_string = 'select sel_name, id '
  for x in cmd.identify(selection,0):
    print(x)
    selection_string += string.strip(str(x) + '+')
  print(selection_string[:-1])

cmd.extend('ens_measure',ens_measure)
cmd.extend('ens_rmsf',ens_rmsf)
cmd.extend('ens_rmsd',ens_rmsd)
cmd.extend('ens_prob',ens_prob)
cmd.extend('print_names',print_names)
