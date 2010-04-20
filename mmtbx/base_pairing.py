
dna_rna_params_str = """
base_pair
  .multiple = True
  .optional = True
  .style = noauto
{
  base1 = None
    .type = str
  base2 = None
    .type = str
  pair_type = *wc
    .type = choice
    .caption = Watson-Crick
# TODO
#  planar = False
#    .type = bool
}
"""
