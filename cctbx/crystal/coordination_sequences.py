from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex

simple = crystal.coordination_sequences_simple
simple_sym = crystal.coordination_sequences_simple_sym
shell_asu_tables = crystal.coordination_sequences_shell_asu_tables

class node(object):

  def __init__(self, asu_mappings, i_seq, rt_mx):
    self.i_seq = i_seq
    self.rt_mx = rt_mx
    self.rt_mx_unique = str(rt_mx.multiply(asu_mappings.special_op(i_seq)))

def find_node(test_node, node_list):
  for list_node in node_list:
    if (    list_node.i_seq == test_node.i_seq
        and list_node.rt_mx_unique == test_node.rt_mx_unique):
      return True
  return False

def simple_and_slow(pair_asu_table, max_shell=10):
  asu_mappings = pair_asu_table.asu_mappings()
  term_table = []
  for i_seq_pivot,pair_asu_dict_pivot in enumerate(pair_asu_table.table()):
    rt_mx_pivot = asu_mappings.get_rt_mx(i_seq=i_seq_pivot, i_sym=0)
    if (pair_asu_dict_pivot.size() == 0):
      term_table.append([])
      continue
    nodes_middle = []
    nodes_next = [node(
      asu_mappings=asu_mappings,
      i_seq=i_seq_pivot,
      rt_mx=sgtbx.rt_mx())]
    terms = [1]
    for i_shell_minus_1 in xrange(max_shell):
      nodes_prev = nodes_middle
      nodes_middle = nodes_next
      nodes_next = []
      for node_m in nodes_middle:
        rt_mx_i = asu_mappings.get_rt_mx(i_seq=node_m.i_seq, i_sym=0)
        rt_mx_ni = node_m.rt_mx.multiply(rt_mx_i.inverse())
        for j_seq,j_sym_groups in pair_asu_table.table()[node_m.i_seq].items():
          for j_sym_group in j_sym_groups:
            for j_sym in j_sym_group:
              rt_mx_j = asu_mappings.get_rt_mx(i_seq=j_seq, i_sym=j_sym)
              new_node = node(
                asu_mappings=asu_mappings,
                i_seq=j_seq,
                rt_mx=rt_mx_ni.multiply(rt_mx_j))
              if (    not find_node(test_node=new_node, node_list=nodes_prev)
                  and not find_node(test_node=new_node, node_list=nodes_middle)
                  and not find_node(test_node=new_node, node_list=nodes_next)):
                nodes_next.append(new_node)
      terms.append(len(nodes_next))
    term_table.append(terms)
  return term_table

def get_kriber_coseq_file(file_name):
  result = {}
  for line in open(file_name):
    flds = line.split()
    tag = flds[0]
    terms = [int(f) for f in flds[1:]]
    result.setdefault(tag, []).append(terms)
  return result

def show_terms(structure, term_table, coseq_dict=None):
  assert len(term_table) == structure.scatterers().size()
  for scatterer,terms in zip(structure.scatterers(), term_table):
    print scatterer.label, list(terms),
    if (coseq_dict is not None):
      terms_to_match = list(terms[1:])
      have_match = False
      tags = coseq_dict.keys()
      tags.sort()
      for tag in tags:
        for coseq_terms in coseq_dict[tag]:
          n = min(len(coseq_terms), len(terms_to_match))
          if (coseq_terms[:n] == terms_to_match[:n]):
            print tag,
            have_match = True
      if (not have_match):
        print "Unknown",
    print
  sums_terms = flex.double()
  multiplicities = flex.double()
  for scatterer,terms in zip(structure.scatterers(), term_table):
    sums_terms.append(flex.sum(flex.size_t(list(terms))))
    multiplicities.append(scatterer.multiplicity())
  print "TD%d: %.2f" % (
    len(terms)-1, flex.mean_weighted(sums_terms, multiplicities))
