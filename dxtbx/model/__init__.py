from __future__ import absolute_import, division
import boost.python
from cctbx import sgtbx # import dependency
from cctbx.crystal_orientation import crystal_orientation # import dependency
from dxtbx_model_ext import *
from dxtbx.model.beam import *
from dxtbx.model.goniometer import *
from dxtbx.model.detector import *
from dxtbx.model.scan import *
from dxtbx.model.crystal import *
from dxtbx.model.profile import *
from libtbx.containers import OrderedSet, OrderedDict

class DetectorAux(boost.python.injector, Detector):
  def iter_panels(self):
    ''' Iterate through just the panels depth-first. '''
    for obj in self.iter_preorder():
      if obj.is_panel():
        yield obj

  def iter_preorder(self):
    ''' Iterate through the groups and panels depth-first. '''
    stack = [self.hierarchy()]
    while (len(stack) > 0):
      node = stack.pop()
      yield node
      if node.is_group():
        for child in reversed(node):
          stack.append(child)

  def iter_levelorder(self):
    ''' Iterate through the groups and panels depth-first. '''
    from collections import deque
    queue = deque([self.hierarchy()])
    while (len(queue) > 0):
      node = queue.popleft()
      yield node
      if node.is_group():
        for child in node:
          queue.append(child)


class CrystalAux(boost.python.injector, Crystal):

  def show(self, show_scan_varying=False, out=None):
    CrystalAux._show(self, show_scan_varying, out)

  @staticmethod
  def _show(self, show_scan_varying=False, out=None):
    from scitbx import matrix
    if out is None:
      import sys
      out = sys.stdout
    uc = self.get_unit_cell().parameters()
    uc_sd = self.get_cell_parameter_sd()
    sg = str(self.get_space_group().info())
    umat = matrix.sqr(self.get_U()).mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()
    bmat = matrix.sqr(self.get_B()).mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()
    amat = (matrix.sqr(self.get_U()) * matrix.sqr(self.get_B())).mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()

    msg =  ["Crystal:"]
    from libtbx.utils import format_float_with_standard_uncertainty
    if len(uc_sd) != 0:
      cell_str = [format_float_with_standard_uncertainty(
          v, e, minimum=1.e-5) for (v, e) in zip(uc, uc_sd)]
      msg.append("    Unit cell: (" + ", ".join(cell_str) + ")")
    else:
      msg.append("    Unit cell: " + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc)
    msg.append("    Space group: " + sg)
    msg.append("    U matrix:  " + umat[0])
    msg.append("               " + umat[1])
    msg.append("               " + umat[2])
    msg.append("    B matrix:  " + bmat[0])
    msg.append("               " + bmat[1])
    msg.append("               " + bmat[2])
    msg.append("    A = UB:    " + amat[0])
    msg.append("               " + amat[1])
    msg.append("               " + amat[2])
    if self.num_scan_points > 0:
      msg.append("    A sampled at " + str(self.num_scan_points) \
                 + " scan points")
      if show_scan_varying:
        for i in range(self.num_scan_points):
          A = matrix.sqr(self.get_A_at_scan_point(i))
          B = matrix.sqr(self.get_B_at_scan_point(i))
          U = matrix.sqr(self.get_U_at_scan_point(i))
          uc = self.get_unit_cell_at_scan_point(i).parameters()
          umat = U.mathematica_form(format="% 5.4f",
                                    one_row_per_line=True).splitlines()
          bmat = B.mathematica_form(format="% 5.4f",
                                    one_row_per_line=True).splitlines()
          amat = A.mathematica_form(format="% 5.4f",
                                    one_row_per_line=True).splitlines()
          msg.append("  Scan point #%i:" %(i+1))
          msg.append("    Unit cell: " + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc)
          msg.append("    U matrix:  " + umat[0])
          msg.append("               " + umat[1])
          msg.append("               " + umat[2])
          msg.append("    B matrix:  " + bmat[0])
          msg.append("               " + bmat[1])
          msg.append("               " + bmat[2])
          msg.append("    A = UB:    " + amat[0])
          msg.append("               " + amat[1])
          msg.append("               " + amat[2])
    print >> out, "\n".join(msg)

  def __str__(self):
    from cStringIO import StringIO
    s = StringIO()
    msg = self.show(out=s)
    s.seek(0)
    return s.read()

  @staticmethod
  def _to_dict(crystal):
    ''' Convert the crystal model to a dictionary

    Params:
        crystal The crystal model

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    from scitbx import matrix

    # Get the real space vectors
    A = matrix.sqr(crystal.get_A()).inverse()
    real_space_a = (A[0], A[1], A[2])
    real_space_b = (A[3], A[4], A[5])
    real_space_c = (A[6], A[7], A[8])

    # Get the space group Hall symbol
    hall = crystal.get_space_group().info().type().hall_symbol()

    # Isoforms used for stills
    try:
      identified_isoform = crystal.identified_isoform
    except AttributeError:
      identified_isoform = None

    # Collect the information as a python dictionary
    xl_dict = OrderedDict([
      ('__id__', 'crystal'),
      ('real_space_a', real_space_a),
      ('real_space_b', real_space_b),
      ('real_space_c', real_space_c),
      ('space_group_hall_symbol', hall)])

    if identified_isoform is not None:
      xl_dict['identified_isoform'] = identified_isoform

    # Add in scan points if present
    if crystal.num_scan_points > 0:
      A_at_scan_points = tuple([crystal.get_A_at_scan_point(i) \
                                for i in range(crystal.num_scan_points)])
      xl_dict['A_at_scan_points'] = A_at_scan_points

    # Add in covariance of B if present
    cov_B = tuple(crystal.get_B_covariance())
    if len(cov_B) != 0:
      xl_dict['B_covariance'] = cov_B

    return xl_dict

  def to_dict(crystal):
    return CrystalAux._to_dict(crystal)

  @staticmethod
  def from_dict(d):
    ''' Convert the dictionary to a crystal model

    Params:
        d The dictionary of parameters

    Returns:
        The crystal model

    '''
    from dxtbx.model import Crystal

    # If None, return None
    if d is None:
      return None

    # Check the version and id
    if str(d['__id__']) != "crystal":
      raise ValueError("\"__id__\" does not equal \"crystal\"")

    # Extract from the dictionary
    real_space_a = d['real_space_a']
    real_space_b = d['real_space_b']
    real_space_c = d['real_space_c']
    # str required to force unicode to ascii conversion
    space_group  = str("Hall:" + d['space_group_hall_symbol'])
    xl = Crystal(real_space_a, real_space_b, real_space_c,
                       space_group_symbol=space_group)

    # Isoforms used for stills
    try:
      xl.identified_isoform = d['identified_isoform']
    except KeyError:
      pass

    # Extract scan point setting matrices, if present
    try:
      A_at_scan_points = d['A_at_scan_points']
      xl.set_A_at_scan_points(A_at_scan_points)
    except KeyError:
      pass

    # Extract covariance of B, if present
    try:
      cov_B = d['B_covariance']
      xl.set_B_covariance(cov_B)
    except KeyError:
      pass

    return xl

class MosaicCrystalKabsch2010Aux(CrystalAux, MosaicCrystalKabsch2010):
  def show(self, show_scan_varying=False, out=None):
    from scitbx import matrix
    CrystalAux._show(self, show_scan_varying, out)

    if out is None:
      import sys
      out = sys.stdout

    msg = []
    msg.append("    Mosaicity:  %.6f"%self.get_mosaicity())

    print >> out, "\n".join(msg)

  def __str__(self):
    from cStringIO import StringIO
    s = StringIO()
    msg = self.show(out=s)
    s.seek(0)
    return s.read()

  def to_dict(crystal):
    ''' Convert the crystal model to a dictionary

    Params:
        crystal The crystal model

    Returns:
        A dictionary of the parameters

    '''
    xl_dict = CrystalAux._to_dict(crystal)


    # Get the mosaicity
    mosaicity = crystal.get_mosaicity()
    xl_dict['mosaicity'] = mosaicity

    return xl_dict

  @staticmethod
  def from_dict(d):
    ''' Convert the dictionary to a crystal model

    Params:
        d The dictionary of parameters

    Returns:
        The crystal model

    '''
    xl = MosaicCrystalKabsch2010(CrystalAux.from_dict(d))

    # This parameter doesn't survive the Crystal copy constructor so has to be re-set.
    # Isoforms used for stills
    try:
      xl.identified_isoform = d['identified_isoform']
    except KeyError:
      pass

    # Extract mosaicity, if present
    try:
      mosaicity = d['mosaicity']
      xl.set_mosaicity(mosaicity)
    except KeyError:
      pass

    return xl

class MosaicCrystalSauter2014Aux(CrystalAux, MosaicCrystalSauter2014):
  def show(self, show_scan_varying=False, out=None):
    from scitbx import matrix
    CrystalAux._show(self, show_scan_varying, out)

    if out is None:
      import sys
      out = sys.stdout

    msg = []
    msg.append("    Half mosaic angle (degrees):  %.6f"%self.get_half_mosaicity_deg())
    msg.append("    Domain size (Angstroms):  %.6f"%self.get_domain_size_ang())

    print >> out, "\n".join(msg)

  def get_A_as_sqr(self): # required for lunus
    from scitbx.matrix import sqr
    return sqr(self.get_A())

  def get_A_inverse_as_sqr(self):
    return self.get_A_as_sqr().inverse()

  def __str__(self):
    from cStringIO import StringIO
    s = StringIO()
    msg = self.show(out=s)
    s.seek(0)
    return s.read()

  def to_dict(crystal):
    ''' Convert the crystal model to a dictionary

    Params:
        crystal The crystal model

    Returns:
        A dictionary of the parameters

    '''
    xl_dict = CrystalAux._to_dict(crystal)

    # Get the mosaic parameters
    half_mosaicity = crystal.get_half_mosaicity_deg()
    xl_dict['ML_half_mosaicity_deg'] = half_mosaicity

    domain_size = crystal.get_domain_size_ang()
    xl_dict['ML_domain_size_ang'] = domain_size

    return xl_dict

  @staticmethod
  def from_dict(d):
    ''' Convert the dictionary to a crystal model

    Params:
        d The dictionary of parameters

    Returns:
        The crystal model

    '''
    xl = MosaicCrystalSauter2014(CrystalAux.from_dict(d))

    # Parameters for maximum likelihood values
    try:
      xl.set_half_mosaicity_deg(d['ML_half_mosaicity_deg'])
    except KeyError:
      pass
    try:
      xl.set_domain_size_ang(d['ML_domain_size_ang'])
    except KeyError:
      pass

    # This parameter doesn't survive the Crystal copy constructor so have to be re-set
    # Isoforms used for stills
    try:
      xl.identified_isoform = d['identified_isoform']
    except KeyError:
      pass

    return xl

class ExperimentListAux(boost.python.injector, ExperimentList):

  def beams(self):
    ''' Get a list of the unique beams (includes None). '''
    return list(OrderedSet(e.beam for e in self))

  def detectors(self):
    ''' Get a list of the unique detectors (includes None). '''
    return list(OrderedSet(e.detector for e in self))

  def goniometers(self):
    ''' Get a list of the unique goniometers (includes None). '''
    return list(OrderedSet(e.goniometer for e in self))

  def scans(self):
    ''' Get a list of the unique scans (includes None). '''
    return list(OrderedSet(e.scan for e in self))

  def crystals(self):
    ''' Get a list of the unique crystals (includes None). '''
    return list(OrderedSet(e.crystal for e in self))

  def profiles(self):
    ''' Get a list of the unique profile models (includes None). '''
    return list(OrderedSet(e.profile for e in self))

  def scaling_models(self):
    ''' Get a list of the unique scaling models (includes None). '''
    return list(OrderedSet(e.scaling_model for e in self))


  def imagesets(self):
    ''' Get a list of the unique imagesets (includes None).

    This returns unique complete sets rather than partial.
    '''
    return list(OrderedSet([e.imageset for e in self if e.imageset is not None]))
    # temp = OrderedDict([(e.imageset.reader(), i)
    #   for i, e in enumerate(self) if e.imageset is not None])
    # return OrderedDict([(self[i].imageset.complete_set(), None)
    #   for i in temp.itervalues()]).keys()

  def all_stills(self):
    ''' Check if all the experiments are stills '''
    def is_still(e):
      if e.goniometer is None or e.scan is None or e.scan.get_oscillation()[1] == 0:
        return True
      return False
    assert(len(self) > 0)
    result = True
    for e in self:
      if not is_still(e):
        result = False
        break

    return result

  def all_sweeps(self):
    ''' Check if all the experiments are from sweeps '''
    def is_still(e):
      if e.goniometer is None or e.scan is None or e.scan.get_oscillation()[1] == 0:
        return True
      return False
    assert(len(self) > 0)
    result = True
    for e in self:
      if is_still(e):
        result = False
        break
    return result

  def to_dict(self):
    ''' Serialize the experiment list to dictionary. '''
    from dxtbx.imageset import ImageSet, ImageSweep, ImageGrid

    # Check the experiment list is consistent
    assert(self.is_consistent())

    # Get the list of unique models
    blist = self.beams()
    dlist = self.detectors()
    glist = self.goniometers()
    slist = self.scans()
    clist = self.crystals()
    ilist = self.imagesets()
    plist = self.profiles()
    scalelist = self.scaling_models()

    # Create the output dictionary
    result = OrderedDict()
    result['__id__'] = 'ExperimentList'
    result['experiment'] = []

    # Function to find in list by reference
    def find_index(l, m):
      for i, mm in enumerate(l):
        if mm is m:
          return i
      return -1

    # Add the experiments to the dictionary
    for e in self:
      obj = OrderedDict()
      obj['__id__'] = 'Experiment'
      if e.beam is not None:
        obj['beam'] = find_index(blist, e.beam)
      if e.detector is not None:
        obj['detector'] = find_index(dlist, e.detector)
      if e.goniometer is not None:
        obj['goniometer'] = find_index(glist, e.goniometer)
      if e.scan is not None:
        obj['scan'] = find_index(slist, e.scan)
      if e.crystal is not None:
        obj['crystal'] = find_index(clist, e.crystal)
      if e.profile is not None:
        obj['profile'] = find_index(plist, e.profile)
      if e.scaling_model is not None:
        obj['scaling_model'] = find_index(scalelist, e.scaling_model)
      if e.imageset is not None:
        same = [e.imageset is imset for imset in ilist]
        assert(same.count(True) == 1)
        obj['imageset'] = same.index(True)
        if e.scan is None and not isinstance(e.imageset, ImageSweep):
          if len(e.imageset) != len(e.imageset.complete_set()):
            obj['imageset'] = obj['imageset']
      result['experiment'].append(obj)

    def get_template(imset):
      if imset.reader().is_single_file_reader():
        return imset.reader().master_path()
      else:
        return imset.get_template()

    # Serialize all the imagesets
    result['imageset'] = []
    for imset in ilist:
      if isinstance(imset, ImageSweep):
        # FIXME_HACK
        template = get_template(imset)
        r = OrderedDict([
          ('__id__', 'ImageSweep'),
          ('template', template)])
      # elif isinstance(imset, MemImageSet):
      #   r = OrderedDict([
      #     ('__id__', 'MemImageSet')])
        if imset.reader().is_single_file_reader():
          r['single_file_indices'] = list(imset.indices())
      elif isinstance(imset, ImageSet):
        r = OrderedDict([
          ('__id__', 'ImageSet'),
          ('images', imset.paths())])
        if imset.reader().is_single_file_reader():
          r['single_file_indices'] = list(imset.indices())
      elif isinstance(imset, ImageGrid):
        r = OrderedDict([
          ('__id__', 'ImageGrid'),
          ('images', imset.paths()),
          ('grid_size', imset.get_grid_size())])
        if imset.reader().is_single_file_reader():
          r['single_file_indices'] = list(imset.indices())
      else:
        raise TypeError('expected ImageSet or ImageSweep, got %s' % type(imset))
      r['mask'] = imset.external_lookup.mask.filename
      r['gain'] = imset.external_lookup.gain.filename
      r['pedestal'] = imset.external_lookup.pedestal.filename
      r['dx'] = imset.external_lookup.dx.filename
      r['dy'] = imset.external_lookup.dy.filename
      r['params'] = imset.params()
      result['imageset'].append(r)

    # Extract all the model dictionaries
    result['beam'] = [b.to_dict() for b in blist if b is not None]
    result['detector'] = [d.to_dict() for d in dlist if d is not None]
    result['goniometer'] = [g.to_dict() for g in glist if g is not None]
    result['scan'] = [s.to_dict() for s in slist if s is not None]
    result['crystal'] = [c.to_dict() for c in clist if c is not None]
    result['profile'] = [p.to_dict() for p in plist if p is not None]
    result['scaling_model'] = [s.to_dict() for s in scalelist if s is not None]

    # Return the dictionary
    return result

  def to_datablocks(self):
    ''' Return the experiment list as a datablock list.
    This assumes that the experiment contains 1 datablock.'''
    from dxtbx.datablock import DataBlockFactory

    # Convert the experiment list to dict
    obj = self.to_dict()

    # Convert the dictionary to a datablock dictionary
    obj['__id__'] = 'DataBlock'
    for e in obj['experiment']:
      iid = e['imageset']
      imageset = obj['imageset'][iid]
      if 'beam' in e:
        imageset['beam'] = e['beam']
      if 'detector' in e:
        imageset['detector'] = e['detector']
      if 'goniometer' in e:
        imageset['goniometer'] = e['goniometer']
      if 'scan' in e:
        imageset['scan'] = e['scan']

    # Remove the experiments
    del obj['experiment']

    # Create the datablock
    return DataBlockFactory.from_dict([obj])

  def nullify_all_single_file_reader_format_instances(self):
    '''
    Parallel reading of HDF5 from the same handle is not allowed. Python
    multiprocessing is a bit messed up and used fork on linux so need to
    close and reopen file.

    '''
    for experiment in self:
      if experiment.imageset.reader().is_single_file_reader():
        experiment.imageset.reader().nullify_format_instance()
