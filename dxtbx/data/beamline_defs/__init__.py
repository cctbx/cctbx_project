from __future__ import division

def get_beamline_definition(detector_sn_or_key, **kwargs):
  import unicodedata
  import string
  import types
  if not isinstance(detector_sn_or_key, types.UnicodeType):
    detector_sn_or_key = unicode(detector_sn_or_key, 'utf-8', 'ignore')

  valid_chars = frozenset("_.%s%s" % (string.ascii_letters, string.digits))
  filename = unicodedata.normalize('NFKD', detector_sn_or_key).encode('ASCII', 'ignore')
  filename = ''.join(c if c in valid_chars else '_' for c in filename)
  while '__' in filename:
    filename = filename.replace('__', '_')
    # http://stackoverflow.com/questions/1546226/a-simple-way-to-remove-multiple-spaces-in-a-string-in-python/15913564#15913564

  try:
    import importlib
    beamline = importlib.import_module('dxtbx.data.beamline_defs.%s' % filename)
    return beamline.get_definition(**kwargs)
  except ImportError:
    return Dummy("Dummy CIF generator. No information for detector %s (%s.py)" % (detector_sn_or_key, filename))

class template(object):
  def CIF_block(self):
    '''Interface function to generate a CIF block for this detector.'''
    raise RuntimeError('This needs to be overridden.')

  def mmCIF_block(self):
    '''Interface function to generate an mmCIF block for this detector.'''
    raise RuntimeError('This needs to be overridden.')

  def _lookup(self, mmCIFsemantics):
    keys = {
      # Entries can either be
      #
      # a tuple
      #    ( CIF-string, mmCIF-string )
      #
      # or a string such as
      #    "_something?something_else"
      # where
      #   _something_something_else is CIF, and
      #   _something.something_else is mmCIF

      'df.detector':
        ('_diffrn_detector', '_diffrn.detector'),
      'df.rad.mono':       '_diffrn_radiation?monochromator',
      'df.rad.source':     '_diffrn_radiation?source',
      'df.rad.type':       '_diffrn_radiation?type',
      'df.m.dev':          '_diffrn?measurement_device',
      'df.m.dev_type':     '_diffrn?measurement_device_type',
      'df.m.method':       '_diffrn?measurement_method',
      'df.m.spec_supp':    '_diffrn?measurement_specimen_support',
    }
    column = 1 if mmCIFsemantics else 0
    def __lookup(key):
      entry = keys[key]
      if isinstance(entry, basestring):
        return entry.replace('?', '.' if mmCIFsemantics else '_')
      else:
        return entry[column]
    return __lookup

  def write_block_to_file(self, block, filename):
    import iotbx.cif.model
    c = iotbx.cif.model.cif()
    c['cif'] = block
    with open(filename, 'w') as fh:
      c.show(out=fh)

  def _date_to_epoch(self, year, month, day):
    import datetime
    return (datetime.datetime(year,month,day,0,0)
            - datetime.datetime(1970,1,1)).total_seconds()

  def __str__(self):
    if hasattr(self, '_str_override'):
      return self._str_override
    else:
      return "CIF block generator"


class Dummy(template):
  def __init__(self, str_override):
    self._str_override = str_override

  def _dummy_CIF(self):
    import iotbx.cif.model
    return iotbx.cif.model.block()
  CIF_block = _dummy_CIF
  mmCIF_block = _dummy_CIF

