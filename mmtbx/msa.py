
from libtbx import easy_run
import libtbx.load_env
import sys, os

class align_pdb_residues (object) :
  def __init__ (self,
                pdb_sequences,
                pdb_names,
                reference_sequence=None,
                reference_sequence_name=None,
                reference_index=None,
                group_sequences=True) :
    adopt_init_args(self, locals())
    n_models = len(pdb_sequences)
    assert ((n_models > 1) and (n_models == len(pdb_names)))
    if reference_index is None :
      reference_sequence_name = pdb_names[0]
    else :
      assert (reference_sequence is None)
      reference_sequence_name = pdb_names[reference_index]
    fasta = "\n".join([ ">%s\n%s" % (n,s) for n,s in zip(pdb_sequences,
                        pdb_names) ])
    if (reference_sequence is not None) :
      if (reference_sequence_name is None) :
        reference_sequence_name = "seq_file"
      ref_seq_fasta = ">%s\n%s" % (reference_sequence_name, reference_sequence)
      fasta = ref_seq_fasta + fasta
    self.muscle_aln = get_muscle_alignment(fasta, group_sequences)
    self._lookup_table = {}
    for i, name in self.muscle_aln.names :
      pass

  def convert_residue_number (self, pdb_name, resseq) :
    pass

def run_muscle (fasta_sequences, group_sequences=True) :
  if not libtbx.env.has_module(name="muscle") :
    raise RuntimeError("MUSCLE not available or not configured.")
  exe_path = libtbx.env.under_build("muscle/exe/muscle")
  if (os.name == "nt") :
    exe_path += ".exe"
  if (not os.path.isfile(exe_path)) :
    raise RuntimeError("muscle executable is not available.")
  assert isinstance(fasta_sequences, str)
  cmd = "%s -quiet -clw" % exe_path
  if group_sequences :
    cmd += " -group"
  else :
    cmd += " -stable"
  muscle_out = easy_run.fully_buffered(cmd,
    stdin_lines=fasta_sequences)
  return muscle_out.stdout_lines

def get_muscle_alignment (fasta_sequences, group_sequences=True) :
  muscle_out = run_muscle(fasta_sequences, group_sequences)
  from iotbx.bioinformatics import clustal_alignment_parse
  alignment, null = clustal_alignment_parse("\n".join(muscle_out))
  return alignment

def exercise () :
  fasta_sequences = """\
>1MRU_A
GSHMTTPSHLSDRYELGEILGFGGMSEVHLARDLRLHRDVAVKVLRADLARDPSFYLRFRREAQNAAALNHPAIVAVYDT
GEAETPAGPLPYIVMEYVDGVTLRDIVHTEGPMTPKRAIEVIADACQALNFSHQNGIIHRDVKPANIMISATNAVKVMDF
GIARAIADSGNSVTQTAAVIGTAQYLSPEQARGDSVDARSDVYSLGCVLYEVLTGEPPFTGDSPVSVAYQHVREDPIPPS
ARHEGLSADLDAVVLKALAKNPENRYQTAAEMRADLVRVHNGEPPEAPKVLTDAERTSLLSSAAGNLSGPR
>1L3R_E
GNAAAAKKGSEQESVKEFLAKAKEDFLKKWETPSQNTAQLDQFDRIKTLGTGSFGRVMLVKHKESGNHYAMKILDKQKVV
KLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVAGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSL
DLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFA
DQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPKFKG
PGDTSNFDDYEEEEIRVSINEKCGKEFTEF
>2A19_B
GAHTVDKRFGMDFKEIELIGSGGFGQVFKAKHRIDGKTYVIKRVKYNNEKAEREVKALAKLDHVNIVHYNGCWDGFDYDP
ETSSKNSSRSKTKCLFIQMEFCDKGTLEQWIEKRRGEKLDKVLALELFEQITKGVDYIHSKKLINRDLKPSNIFLVDTKQ
VKIGDFGLVTSLKNDGKRTRSKGTLRYMSPEQISSQDYGKEVDLYALGLILAELLHVCDTAFETSKFFTDLRDGIISDIF
DKKEKTLLQKLLSKKPEDRPNTSEILRTLTVWKKSPEKNERHTA"""
  if (not libtbx.env.has_module(name="muscle")) :
    print "Skipping MUSCLE tests: muscle module not available."
    return
  clustal_aln = get_muscle_alignment(fasta_sequences)
  assert (clustal_aln.names == ['1MRU_A', '2A19_B', '1L3R_E'])
  assert (len(clustal_aln.alignments[0]) == 371)
  print "OK"
  return clustal_aln

if __name__ == "__main__" :
  exercise()
