from __future__ import division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 11/05/2018
Description : wxPython 3-4 compatibility tools
'''

from wx import __version__ as wxv
from contextlib import contextmanager
import importlib

wx4 = wxv[0] == '4'

modnames = [
  ('PyControl', 'Control'),
  ('PyDataObjectSimple', 'DataObjectSimple'),
  ('PyDropTarget', 'DropTarget'),
  ('PyEvtHandler', 'EvtHandler'),
  ('PyImageHandler', 'ImageHandler'),
  ('PyLocale', 'Locale'),
  ('PyLog', 'Log'),
  ('PyPanel', 'Panel'),
  ('PyPickerBase', 'PickerBase'),
  ('PyPreviewControlBar', 'PreviewControlBar'),
  ('PyPreviewFrame', 'PreviewFrame'),
  ('PyPrintPreview', 'PrintPreview'),
  ('PyScrolledWindow', 'ScrolledWindow'),
  ('PySimpleApp', 'App'),
  ('PyTextDataObject', 'TextDataObject'),
  ('PyTimer', 'Timer'),
  ('PyTipProvider', 'adv.TipProvider'),
  ('PyValidator', 'Validator'),
  ('PyWindow'', Window')
]

def find_module(module):
  for m in modnames:
    if module.__name__ in m:
      return m

def get_wx_mod(base, module):
  mname = find_module(module)[1] if wx4 else find_module(module)[0]
  bname = base.__name__
  if '.' in mname:
    spl = [i for i in mname.split('.') if i != bname]
    modname = '.'.join(spl[:-1])
    mod = importlib.import_module('{}.{}'.format(bname, modname))
    return getattr(mod, spl[-1])
  else:
    return getattr(base, mname)

@contextmanager
def wx_mod(base, module):
  ''' Identify and import the appropriate wxPython module '''
  yield get_wx_mod(base, module)