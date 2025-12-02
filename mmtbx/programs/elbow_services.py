from __future__ import absolute_import, division, print_function

import os

from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry

# from mmtbx.ligands import elbow_utils
from mmtbx.monomer_library import geostd_utils
from mmtbx.monomer_library import mmtbx_utils
from mmtbx import chemical_components

sources='geostd monlib ccd user'

sources_long = {'geostd' : 'GeoStd',
                'monlib' : 'Monomer Library',
                'ccd' : 'Chem. Comp. Lib.',
                'user' : 'User defined',
                }

class Program_where_is_that_cif_file(ProgramTemplate):

  description = '''
mmtbx.where_is_that_cif_file:

Usage examples:
  mmtbx.where_is_that_cif_file <code>
  '''

  datatypes = ['phil']

  master_phil_str = """
  where {
    code = None
      .type = str
    sources = *%s all
      .type = choice(multi=True)
    display = *full compact ultra oneline
      .type = choice(multi=False)
      .short_caption = Amount of text to "compress" from file
    format = *cif mol2 frcmod amber
      .type = choice(multi=False)
      .short_caption = Check if there are Amber files
    pH = *None low
      .type = choice(multi=False)
      .short_caption = Check for neutron possible protonations
    user_dir = None
      .type = path
      .short_caption = directory to look for restraints
    name = None
      .type = str
  }
""" % sources

  # ---------------------------------------------------------------------------
  def get_user_dir(self):
    d=os.environ.get('TEST', None)
    if self.params.where.user_dir:
      d=self.params.where.user_dir
    return d

  def validate(self):
    print('Validating inputs', file=self.logger)
    assert self.params.where.sources
    code=self.params.where.code
    print(code)
    for i in ['*', '?']:
      if code.find(i)>-1:
        raise Sorry(f'Character "{i}" not recognized. Use "." for wildcard.')
    if len(code)==1 and code=='.':
      raise Sorry('Use a start character to search codes with that start. eg "A."')
    if self.params.where.pH=='None':
      self.params.where.pH=None
    if self.params.where.pH and 'cif' not in self.params.where.format:
      raise Sorry(self.params.where.format=='cif', 'pH=%s only for cif format' % self.params.where.pH)
    if 'user' in self.params.where.sources:
      if not self.get_user_dir():
        raise Sorry('Directory for "user" source not found')

  def _find_filename_from_source(self, code, source, pH=None, file_format='cif'):
    filename=None
    if source=='geostd':
      filename = geostd_utils.get_geostd_file(code, pH=pH, file_format=file_format)
    elif source=='monlib':
      filename = mmtbx_utils.get_monomer_cif_file(code)
    elif source=='user':
      filename=None
      d=self.get_user_dir()
      if d is None: return filename
      d=os.path.abspath(d)
      for preamble in ['', 'data_']:
        filename = os.path.join(d,code[0].lower(), '%s%s.%s' % (preamble, code.upper(), file_format))
        if os.path.exists(filename):
          break
      else:
        filename=None
    elif source=='ccd' and file_format=='cif':
      filename = chemical_components.get_cif_filename(code)
      if not os.path.exists(filename):
        filename = None
    return filename

  def update_results(self, filename, source):
    if filename is None: return
    f=open(filename, 'r')
    lines=f.read()
    del f
    self.results.setdefault(source, {})
    self.results[source]['filename']=filename
    self.results[source]['lines']=lines
    print('%s\n  %s %s\n%s' % ('-'*80, sources_long[source], filename, '-'*80),
          file=self.logger)
    return lines

  # ---------------------------------------------------------------------------
  def run(self, log=None):
    from iotbx import cif
    self.results={}
    code=self.params.where.code
    displayed=False
    if 'all' in self.params.where.sources:
      self.params.where.sources=sources.split()
      # self.params.where.display='oneline'
    if code.find('.')>-1:
      codes=chemical_components.get_filenames_from_start(code[:-1])
      self.results['codes']=codes
      outl='  '
      for i, code in enumerate(codes):
        outl+=f'  {code:5}'
        if i%10==9:
          print(outl, file=self.logger)
          outl='  '
      print(outl, file=self.logger)
      return
    for i, source in enumerate(self.params.where.sources):
      if self.params.where.format=='amber':
        filename=self._find_filename_from_source(code,
                                                 source,
                                                 pH=self.params.where.pH,
                                                 file_format='mol2',
                                                 )
        if filename:
          print('Amber files', file=self.logger)
          print('  mol2   : %s'%filename, file=self.logger)
          filename=self._find_filename_from_source(code,
                                                   source,
                                                   pH=self.params.where.pH,
                                                   file_format='frcmod',
                                                   )
          print('  frcmod : %s'%filename, file=self.logger)
          continue
      else:
        filename=self._find_filename_from_source(code,
                                                 source,
                                                 pH=self.params.where.pH,
                                                 file_format=self.params.where.format,
                                                 )
      # display
      for ts, to in [['monlib', '  Checking "Monomer Library"'],
                     ['ccd',    '  Checking "Chem. Comp. Lib."'],
                    ]:
        if (filename is None and
            self.params.where.format == 'cif' and
            self.params.where.pH is None):
          print('  File not found in "%s"' % sources_long[source], file=self.logger)
          if ts not in self.params.where.sources and len(self.params.where.sources)==i+1:
            print(to, file=self.logger)
            source=ts
            filename=self._find_filename_from_source(code, source)
          if filename: break
      if filename is None:
        continue
      displayed=True
      lines=self.update_results(filename, source)
      if self.params.where.display=='full':
        if self.params.where.name:
          for line in lines.splitlines():
            if line.find('loop_')>-1: print(line, file=self.logger)
            if line.find('_chem_comp.')>-1: continue
            if line.strip().find('_')==0: print(line, file=self.logger)
            if line.find(' %s ' % self.params.where.name)>-1:
              print(line, file=self.logger)
        else:
          print(lines, file=self.logger)
      elif self.params.where.display in ['compact', 'ultra']:
        outl=''
        loop=''
        for line in lines.splitlines():
          if not line.strip(): continue
          if line.find('#')==0: continue
          if self.params.where.display=='ultra':
            if line.strip()[0]=='_':
              # if outl and outl.split()[0]!='|':
              #   print(outl, file=self.logger)
              #   outl=''
              if len(line.split())>1:
                loop += ' %s\n' % line
              elif loop:
                attr = line.split('.')[-1]
                if len(loop+attr)>79:
                  print(loop, file=self.logger)
                  loop=line
                else:
                  loop += ' %s' % attr
              else:
                loop=line
              continue
            else:
              if loop:
                print(loop, file=self.logger)
                loop=''
            if line=='loop_': continue
            if len(outl+line)>75:
              print(outl, file=self.logger)
              outl=f'  {line}'
            else:
              outl+=f' | {line}'
          else:
            print(line, file=self.logger)
      elif self.params.where.display in ['oneline']:
        pass
      if self.params.where.display not in ['oneline']:
        print('%s\n  %s %s\n%s' % ('-'*80, sources_long[source], filename, '-'*80), file=self.logger)

    s, r = chemical_components.get_obsolete_status_from_chemical_components(code)
    if s == "OBS":
      print("\n  Code '%s' has been replaced by '%s'\n" % (code.upper(), r),
            file=self.logger)
      return

    if self.params.where.format!='cif': return
    if self.params.where.pH is not None: return

    if not displayed:
      tmp={}
      geostd_keys=None
      monlib_keys=None
      for filename in [
        "geostd_list.cif",
        "mon_lib_list.cif",
        ]:
        if filename=="geostd_list.cif":
          cif_list = geostd_utils.get_cif_list(filename)
          geostd_keys=list(cif_list.keys())
        else:
          cif_list = geostd_utils.get_cif_list(filename)
          monlib_keys=list(cif_list.keys())
        tmp.update(cif_list)
        hits=[]
        for key in tmp.keys():
          if key.upper().find(code.upper())>-1:
            hits.append(key)
      if len(hits)==0:
        print(f'\nNo hits found for "{code}".', file=self.logger)
        outl='  '
        for i, key in enumerate(tmp.keys()):
          outl+=f' {key:25}'
          if i%3==2: outl+='\n  '
        print(outl, file=self.logger)
      elif len(hits)==1:
        if hits[0] in geostd_keys:
          source='GeoStd'
        else:
          source='Monomer Library'
        print('\nFound %s in %s list.cif\n' % (hits[0], source), file=self.logger)
        print(tmp[hits[0]], file=self.logger)
        print('\nFound %s in %s list.cif\n' % (hits[0], source), file=self.logger)
      else:
        print('\n  Multiple choices: Make a more precise search.', file=self.logger)
        for key in hits:
          print(f'    {key}', file=self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results
