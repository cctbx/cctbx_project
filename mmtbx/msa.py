
from libtbx import easy_run
import libtbx.load_env
from libtbx import adopt_init_args
import sys, os

class align_pdb_residues (object) :
  def __init__ (self,
                pdb_sequences,
                pdb_names,
                pdb_offsets,
                reference_sequence=None,
                reference_sequence_name="sequence",
                reference_sequence_offset=0,
                reference_index=None) :
#                group_sequences=True) :
    adopt_init_args(self, locals())
    n_models = len(pdb_sequences)
    assert ((n_models >= 1) and (n_models==len(pdb_names)==len(pdb_offsets)))
    use_pdb_sequence = False
    if (reference_sequence is not None) :
      assert (reference_index is None)
    else :
      use_pdb_sequence = True
      if (reference_index is None) :
        reference_index = 0
      reference_sequence = pdb_sequences[reference_index]
      reference_sequence_name = pdb_names[reference_index]
      reference_sequence_offset = pdb_offsets[reference_index]
    fasta = "\n".join([ ">%s\n%s" % (n,s) for n,s in zip(pdb_names,
                        pdb_sequences) ])
    if (not use_pdb_sequence) :
      ref_seq_fasta = ">%s\n%s" % (reference_sequence_name, reference_sequence)
      fasta = ref_seq_fasta + "\n" + fasta
    self.muscle_aln = get_muscle_alignment(fasta)#, group_sequences)
    assert (self.muscle_aln is not None)
    self._lookup_table = {}
    i_ref = self.muscle_aln.names.index(reference_sequence_name)
    self._reference_alignment = self.muscle_aln.alignments[i_ref]
    assert (i_ref is not None)
    for i, name in enumerate(self.muscle_aln.names) :
      if (i == i_ref) :
        self._lookup_table[name] = None
      else :
        i_pdb = self.pdb_names.index(name)
        aln = self.muscle_aln.alignments[i]
        j = k = 0
        new_resseqs = []
        for res1, res2 in zip(aln, self._reference_alignment) :
          if (res2 != '-') :
            k += 1
          if (res1 == '-') :
            continue
          else :
            if (res2 == '-') :
              new_resseqs.append(None)
            else :
              new_resseqs.append(k - reference_sequence_offset)
            j += 1
        assert (j == len(new_resseqs))
        self._lookup_table[name] = new_resseqs

  def convert_residue_number (self, pdb_name, resseq) :
    assert isinstance(resseq, int)
    if (self._lookup_table[pdb_name] is None) :
      return resseq
    i_pdb = self.pdb_names.index(pdb_name)
    offset = self.pdb_offsets[i_pdb]
    i_res = resseq + offset - 1
    return self._lookup_table[pdb_name][i_res]

def run_muscle (fasta_sequences, group_sequences=True) :
  assert group_sequences # XXX this isn't actually optional!
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
  #else :
  #  cmd += " -stable"
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
  pdb_sequences = [
    "XXATKGQNEPKKHFMVILSTCYGWSDD",
    "VTKGNDPKRHYIVLLTGAGSCFGWSDDARY",
    "GSHXXATKGQNDPKRHYMILSCYGWDDRY",
  ]
  pdb_names = ["pdb1","pdb2","pdb3"]
  pdb_offsets = [0, 0, 3]
  reference_sequence = "MIATKGQNEPKKHFMVILSTCYGWSDDAKFGG"
  m = align_pdb_residues(
    pdb_sequences=pdb_sequences,
    pdb_names=pdb_names,
    pdb_offsets=pdb_offsets,
    reference_sequence=reference_sequence)
  assert (m.convert_residue_number("pdb3", -2) is None)
  assert (m.convert_residue_number("pdb1", 10) == 10)
  assert (m.convert_residue_number("pdb2", 5) == 8)
  assert (m.convert_residue_number("pdb2", 16) == 19)
  assert (m.convert_residue_number("pdb2", 18) is None)
  assert (m.convert_residue_number("pdb1", 27) == 27)
  assert (m.convert_residue_number("pdb2", 30) == 30)
  #print m.muscle_aln.format(80,10)
  print "OK"

if __name__ == "__main__" :
  exercise()
