from __future__ import division
from xfel.ui.db import db_proxy
from scitbx.array_family import flex

class Event(db_proxy):
  def __init__(self, app, event_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_event" % app.params.experiment_tag, id = event_id, **kwargs)
    self.event_id = self.id

class Experiment(db_proxy):
  def __init__(self, app, experiment_id = None, experiment = None, cell = None, **kwargs):
    assert [experiment_id, experiment].count(None) == 1
    if experiment is not None:
      self.imageset = Imageset(app)
      self.beam = Beam(app, beam = experiment.beam)
      self.detector = Detector(app, detector = experiment.detector)
      self.crystal = Crystal(app, crystal = experiment.crystal, cell = cell)

      kwargs['imageset_id'] = self.imageset.id
      kwargs['beam_id'] = self.beam.id
      kwargs['detector_id'] = self.detector.id
      kwargs['crystal_id'] = self.crystal.id
      kwargs['crystal_cell_id'] = self.crystal.cell_id

    db_proxy.__init__(self, app, "%s_experiment" % app.params.experiment_tag, id = experiment_id, **kwargs)
    self.experiment_id = self.id

    if experiment is None:
      self.imageset = Imageset(app, imageset_id=self.imageset_id)
      self.beam = Beam(app, beam_id=self.beam_id)
      self.detector = Detector(app, self.detector_id)
      self.crystal = Crystal(app, self.crystal_id)

class Imageset(db_proxy):
  def __init__(self, app, imageset_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_imageset" % app.params.experiment_tag, id=imageset_id, **kwargs)
    self.imageset_id = self.id

class Beam(db_proxy):
  def __init__(self, app, beam_id = None, beam = None, **kwargs):
    assert [beam_id, beam].count(None) == 1
    if beam is not None:
      u_s0 = beam.get_unit_s0()
      kwargs['direction_1'] = u_s0[0]
      kwargs['direction_2'] = u_s0[1]
      kwargs['direction_3'] = u_s0[2]
      kwargs['wavelength'] = beam.get_wavelength()

    db_proxy.__init__(self, app, "%s_beam" % app.params.experiment_tag, id=beam_id, **kwargs)
    self.beam_id = self.id

class Detector(db_proxy):
  def __init__(self, app, detector_id = None, detector = None, **kwargs):
    assert [detector_id, detector].count(None) == 1
    if detector is not None:
      kwargs['distance'] = flex.mean(flex.double([p.get_distance() for p in detector]))

    db_proxy.__init__(self, app, "%s_detector" % app.params.experiment_tag, id=detector_id, **kwargs)
    self.detector_id = self.id

class Crystal(db_proxy):
  def __init__(self, app, crystal_id = None, crystal = None, cell = None, **kwargs):
    from scitbx import matrix
    assert [crystal_id, crystal].count(None) == 1
    if crystal is not None:
      u = matrix.sqr(crystal.get_U())  # orientation matrix
      for i in xrange(len(u)):
        kwargs['ori_%d' % (i + 1)] = u[i]
      kwargs['mosaic_block_rotation'] = crystal._ML_half_mosaicity_deg
      kwargs['mosaic_block_size'] = crystal._ML_domain_size_ang

      if cell is not None:
        self.cell = cell
      else:
        try:
          isoform_name = crystal.identified_isoform
        except AttributeError:
          self.cell = Cell(app, crystal=crystal, isoform_id = None)
        else:
          tag = app.params.experiment_tag
          query = """SELECT cell.id from `%s_cell` cell
                     JOIN `%s_isoform` isoform ON cell.isoform_id = isoform.id
                     JOIN `%s_trial` trial ON isoform.trial_id = trial.id
                     WHERE isoform.name = '%s' AND trial.trial = %d""" % (
            tag, tag, tag, isoform_name, app.params.input.trial)
          cursor = app.execute_query(query)
          results = cursor.fetchall()
          assert len(results) == 1
          self.cell = Cell(app, cell_id = results[0][0])
      kwargs['cell_id'] = self.cell.id

    db_proxy.__init__(self, app, "%s_crystal" % app.params.experiment_tag, id=crystal_id, **kwargs)
    self.crystal_id = self.id

    if crystal is None:
      assert cell is None
      self.cell = Cell(app, cell_id = self.cell_id)

class Isoform(db_proxy):
  def __init__(self, app, isoform_id=None, **kwargs):
    db_proxy.__init__(self, app, "%s_isoform" % app.params.experiment_tag, id=isoform_id, **kwargs)
    self.isoform_id = self.id

  def __getattr__(self, key):
    if key == "cell":
      cells = self.app.get_all_x(Cell, 'cell', where = "WHERE isoform_id = %d"%self.id)
      assert len(cells) == 1
      return cells[0]
    else:
      return super(Isoform, self).__getattr__(key)

class Cell(db_proxy):
  def __init__(self, app, cell_id = None, crystal = None, **kwargs):
    assert [cell_id, crystal].count(None) in [1,2]
    if crystal is not None:
      for key, p in zip(['a', 'b', 'c', 'alpha', 'beta', 'gamma'], crystal.get_unit_cell().parameters()):
        kwargs['cell_%s'%key] = p
      kwargs['lookup_symbol'] = crystal.get_space_group().type().lookup_symbol()
    db_proxy.__init__(self, app, "%s_cell" % app.params.experiment_tag, id=cell_id, **kwargs)
    self.cell_id = self.id

    assert [self.isoform_id, self.trial_id].count(None) <= 1
    if self.isoform_id is not None:
      self.isoform = Isoform(app, isoform_id = self.isoform_id)
      self.bins = app.get_cell_bins(self.id)
    elif self.trial_id is not None:
      self.bins = app.get_cell_bins(self.id)
    else:
      self.isoform = None
      self.bins = []

  def __getattr__(self, key):
    if key == "unit_cell":
      from cctbx.uctbx import unit_cell
      return unit_cell([self.cell_a, self.cell_b, self.cell_c,
                        self.cell_alpha, self.cell_beta, self.cell_gamma])
    else:
      return super(Cell, self).__getattr__(key)

class Bin(db_proxy):
  def __init__(self, app, bin_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_bin" % app.params.experiment_tag, id=bin_id, **kwargs)
    self.bin_id = self.id

class Cell_Bin(db_proxy):
  def __init__(self, app, cell_bin_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_cell_bin" % app.params.experiment_tag, id=cell_bin_id, **kwargs)
    self.cell_bin_id = self.id
