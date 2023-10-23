from __future__ import absolute_import, division, print_function
import os

def classify_histidine(hierarchy, resname='HIS'):
  from mmtbx.validation.rotalyze import rotalyze
  result = rotalyze(
      pdb_hierarchy=hierarchy,
      quiet=False)
  names = []
  for rot in result.results:
    if rot.resname!=resname: continue
    names.append(rot.rotamer_name)
  hs=0
  ha=None
  for atom in hierarchy.atoms():
    if atom.parent().resname!=resname: continue
    if atom.name.strip() in ['HD1', 'HE2']:
      hs+=1
      ha=atom.name.strip()
  assert len(names)==1
  if hs==2: ha = 'HD1, HE2'
  return names[0], ha

def run_hbond(args):
  from iotbx.cli_parser import run_program
  from mmtbx.programs.hbond import Program
  hbonds = run_program(program_class=Program,
                       args=tuple(args),
                       )
  return hbonds

def run_serial_or_parallel(func, argstuples, nproc=1, log=None):
  import time
  from libtbx import easy_mp
  if nproc==-1: assert 0, 'testing using %s' % nproc
  rc = []
  name = func.__name__.replace('_',' ')
  if nproc==1:
    for i, args in enumerate(argstuples):
      t0=time.time()
      # print('  Running "%s" job %d' % (name, i+1), file=log)
      res = func(*args)
      rc.append(res)
      # print('    Time for job %d: %0.1fs' % (i+1, time.time()-t0), file=log)
  elif nproc>1:
    print('  Running %d jobs on %d procs' % (len(argstuples), nproc), file=log)
    i=0
    t0=time.time()
    for args, res, err_str in easy_mp.multi_core_run( func,
                                                      argstuples,
                                                      nproc,
                                                      keep_input_order=True):
      assert not err_str, '\n\nDebug in serial :\n%s' % err_str
      print('  Running "%s" job %d : %0.1fs' % (name, i+1, time.time()-t0), file=log)
      rc.append(res)
      i+=1
  return rc

def get_hbonds_via_filenames(filenames, nq_or_h, nproc=1, restraint_filenames=None):
  assert nproc>0
  argstuples = []
  for i, filename in enumerate(filenames):
    assert os.path.exists(filename), '"%s"' % filename
    argstuples.append([[filename,
                       '--quiet',
                       'output_pymol_file=True',
                       'output_restraint_file=False',
                       'output_skew_kurtosis_plot=False',
                       'min_data_size=0',
                       'prefix=%s' % filename.replace('.pdb',''),
                       ]])
    if restraint_filenames:
      argstuples[-1][-1]+=restraint_filenames

  rc = run_serial_or_parallel(run_hbond, argstuples, nproc=nproc)

  i=0
  hbondss=[]
  pymols = ''
  for i, filename in enumerate(filenames):
    hbondss.append(rc[i])
    pf = filename.replace('.pdb', '.pml')
    assert os.path.exists(pf), 'file not found %s' % pf
    f=open(pf, 'a')
    f.write('\n')
    f.write('show sticks, resn %s\n' % nq_or_h)
    del f
    pymols += '  phenix.pymol %s &\n' % pf
  return hbondss, pymols

def get_rotamers_via_filenames(filenames, selection, resname='HIS'):
  from iotbx import pdb
  rotamers=[]
  for i, filename in enumerate(filenames):
    hierarchy = pdb.input(filename).construct_hierarchy()
    asc1 = hierarchy.atom_selection_cache()
    sel = asc1.selection(selection)
    hierarchy = hierarchy.select(sel)
    rc = classify_histidine(hierarchy, resname=resname)
    rotamers.append(rc[0])
    i+=1
  return rotamers
