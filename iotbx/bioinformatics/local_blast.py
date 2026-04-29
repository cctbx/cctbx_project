"""
This is a tool to run BLAST locally against selected databased such
as PDBaa ...etc.  The executables and databases are all distributed
with Phenix so no extra installation is required. PDBaa has been
implemted.  More classes, e.g. ScopE, PDBstructure ... will be
added later.
LWH 4/12/19

Useage:
>>>from iotbx.bioinformatics import local_blast
>>>myxml=local_blast.pdbaa(seq=myseq).run()

where
myseq is the query protein sequence string. You can use 'X' to fill gaps.
myxml is the stdout_lines object of the blast XML output.
"""

from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import os,sys,libtbx.load_env

#setup lib dir
phenixpath=os.getenv('PHENIX')
ligand_lib_dir = libtbx.env.find_in_repositories(
  relative_path=os.path.join("chem_data","ligand_lib"),
  test=os.path.isdir)
cwd=os.getcwd()

#make sure blast exists.  Never had a problem but probably won't hurt.
def checkblast(binary="blastp"):
  phenix_blast_exe=''
  systype=sys.platform #linux2,darwin,win32
  sysname='Linux'
  phenix_blast_exe='%s_%s'%(binary,systype)
  if systype=='win32':
    phenix_blast_exe='%s_%s.exe'%(binary,systype)
    sysname='Windows'
  elif systype=='darwin':
    phenix_blast_exe='%s_%s'%(binary,systype)
    sysname='OSX'
  elif systype.startswith('linux') and sys.version_info.major == 3:
    systype = 'linux2'
    phenix_blast_exe='%s_%s'%(binary,systype)
  else:
    pass
  phenix_blast=os.path.join(ligand_lib_dir, phenix_blast_exe)
  blastexe=None
  blastpath=None
  if os.path.exists(phenix_blast):
    blastpath=phenix_blast
    blastexe=phenix_blast_exe
    #print('%s version is running...\n'%sysname)
  else:
    print('BLAST executable does not exist. please check your Phenix installation.')
    sys.exit(0)
  return blastpath



class pdbaa(object):
  def __init__(self, workdir=None, prefix=None, seq=None, output=None):
    self.workdir=workdir
    self.prefix=prefix
    self.seq=seq
    self.output=output

  def run(self, debug=False, binary="blastp"):
    blastpath=checkblast(binary)
    curdir=os.getcwd()
    if self.workdir is None:
      self.workdir=curdir
    elif os.path.exists(self.workdir) is False:
      print("Input working directory does not exist. Try to work in current directory instead.")
      self.workdir=curdir
    fasta_file='myprotein.fasta'
    fasta_path=os.path.join(self.workdir,fasta_file)
    if self.prefix is None:
       self.prefix="myprotein"
    fastaline=">%s\n%s\n"%(self.prefix,self.seq)
    f=open(fasta_path,'w').writelines(fastaline)
    dbname="pdbaa.00"
    outfmt="-outfmt 5" #xml_out
    blastdb=os.path.join(ligand_lib_dir,dbname)
    if binary=="blastall":
      blastrun_seq=" -p blastp -i %s -a 8 -F F -W 3 -G 11 -E 2 \
          -V F -e 1E-3 -m 7 -d %s"%(fasta_path, blastdb)
    else:
      blastrun_seq= "-query %s -matrix BLOSUM62 -num_threads 8 -word_size 3 -gapopen 11\
        -gapextend 2 -evalue 1E-3 %s -db %s"%(fasta_path, outfmt, blastdb)

    cmds="%s %s"%(blastpath,blastrun_seq)
    #print(cmds)
    try:
      result = easy_run.fully_buffered(
        command=cmds,
        stdin_lines='')
    except KeyboardInterrupt :
      raise KeyboardInterrupt
    else :
      if debug:
        output='myprotein.xml'
        open(output, "w").write("\n".join(result.stdout_lines))
      return result.stdout_lines
