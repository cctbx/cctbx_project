from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import sys
from libtbx import adopt_init_args
from libtbx.utils import Sorry

class manager(object):
  def __init__(self, individual_sites       = False,
                     torsion_angles         = False,
                     rigid_body             = False,
                     individual_adp         = False,
                     group_adp              = False,
                     tls                    = False,
                     occupancies            = False,
                     group_anomalous        = False,
                     sites_individual       = None,
                     sites_torsion_angles   = None,
                     sites_rigid_body       = None,
                     adp_individual_iso     = None,
                     adp_individual_aniso   = None,
                     adp_group              = None,
                     group_h                = None,
                     adp_tls                = None,
                     s_occupancies          = None):
                     # XXX group_anomalous should be here
    adopt_init_args(self, locals())
    self.sites_individual       = self._deep_copy(self.sites_individual)
    self.sites_torsion_angles   = self._deep_copy(self.sites_torsion_angles)
    self.sites_rigid_body       = self._deep_copy(self.sites_rigid_body)
    self.adp_individual_iso     = self._deep_copy(self.adp_individual_iso)
    self.adp_individual_aniso   = self._deep_copy(self.adp_individual_aniso)
    self.adp_group              = self._deep_copy(self.adp_group)
    self.group_h                = self._deep_copy(self.group_h)
    self.adp_tls                = self._deep_copy(self.adp_tls)
    self.s_occupancies          = self._deep_copy(self.s_occupancies)
    self.check_all()

  def _deep_copy(self, x):
    result = []
    if(x is None): result = x
    elif(self.is_bool(x) or self.is_size_t(x)):
      result = x.deep_copy()
    elif(len(x)==0): result = x
    elif(self.is_size_t(x[0])):
      for item in x:
        if(self.is_size_t(item)):
          result.append(item.deep_copy())
        else: raise RuntimeError("Bad selection array type.")
    elif(self.is_size_t(x[0][0])):
      for item1 in x:
        tmp = []
        for item2 in item1:
          if(self.is_size_t(item2)):
            tmp.append(item2.deep_copy())
          else: raise RuntimeError("Bad selection array type.")
        result.append(tmp)
    else: raise RuntimeError("Bad selection array type.")
    return result

  def deep_copy(self):
    return manager(
      individual_sites       = self.individual_sites,
      torsion_angles         = self.torsion_angles,
      rigid_body             = self.rigid_body,
      individual_adp         = self.individual_adp,
      group_adp              = self.group_adp,
      tls                    = self.tls,
      occupancies            = self.occupancies,
      group_anomalous        = self.group_anomalous,
      sites_individual       = self._deep_copy(self.sites_individual),
      sites_torsion_angles   = self._deep_copy(self.sites_torsion_angles),
      sites_rigid_body       = self._deep_copy(self.sites_rigid_body),
      adp_individual_iso     = self._deep_copy(self.adp_individual_iso),
      adp_individual_aniso   = self._deep_copy(self.adp_individual_aniso),
      adp_group              = self._deep_copy(self.adp_group),
      group_h                = self._deep_copy(self.group_h),
      adp_tls                = self._deep_copy(self.adp_tls),
      s_occupancies          = self._deep_copy(self.s_occupancies))
      # XXX group_anomalous should be here

  def is_size_t(self, x):
    return ("%s"%x.__class__).count("array_family_flex_ext.size_t") > 0

  def is_bool(self, x):
    return ("%s"%x.__class__).count("array_family_flex_ext.bool") > 0

  def _count_selected(self, selections):
    assert selections is not None
    selections = self._deep_copy(selections)
    result = True
    try: lx = len(selections)
    except Exception: lx = selections.size()
    if(lx == 0): return result
    if(self.is_bool(selections)):
      if(selections.count(True) == 0): result = False
    elif(self.is_size_t(selections)):
      if(selections.size() == 0): result = False
    elif(self.is_size_t(selections[0]) or self.is_bool(selections[0])):
      as_1d = []
      for sel in selections:
        if(self.is_size_t(sel)):
          as_1d.extend(list(sel))
        elif(self.is_bool(sel)):
          as_1d.extend(list(sel.iselection()))
        else:
          result = False
          break
      if(len(as_1d) == 0 or len(as_1d) != len(set(as_1d))):
        result = False
    else:
      as_1d = []
      for i in selections:
        for j in i:
          for k in j:
            as_1d.append(k)
      if(len(as_1d) == 0 or len(as_1d) != len(set(as_1d))):
        result = False
    return result

  def check_all(self):
    prefix = "\nBad (empty or mixed) selection in %s"
    if(self.individual_sites):
      if(not self._count_selected(self.sites_individual)):
        raise Sorry(prefix%"sites_individual.")
    if(self.torsion_angles):
      if(not self._count_selected(self.sites_torsion_angles)):
        raise Sorry(prefix%"sites_torsion_angles.")
    if(self.rigid_body):
      if(not self._count_selected(self.sites_rigid_body)):
        raise Sorry(prefix%"sites_rigid_body.")
    if(self.individual_adp):
      if(not self.tls and [self.adp_individual_iso,
         self.adp_individual_aniso].count(None)==0):
        if(not (self.adp_individual_aniso &
           self.adp_individual_iso).all_eq(False)):
          raise Sorry("Same atoms selected for iso and aniso ADP refinement.")
      elif(self.adp_individual_iso is not None):
        if(not self._count_selected(self.adp_individual_iso)):
          if(not self._count_selected(self.adp_individual_aniso)):
            raise Sorry(prefix%"adp_individual_iso.")
      elif(self.adp_individual_aniso is not None):
        if(not self._count_selected(self.adp_individual_aniso)):
          raise Sorry(prefix%"adp_individual_aniso.")
      else: raise Sorry("No selection for individual_adp.")
    if(self.group_adp):
      if(not self._count_selected(self.adp_group)):
        raise Sorry(prefix%"adp_group.")
    if(self.tls):
      if(not self._count_selected(self.adp_tls)):
        raise Sorry(prefix%"adp_tls.")
    if(self.occupancies and self.s_occupancies is not None):
      if(not self._count_selected(self.s_occupancies)):
        raise Sorry(prefix%"occupancies.")
    if(self.group_anomalous):
      pass # XXX selections not used in common framework

  def szs(self, x):
    if(x is not None): return str(len(x))
    else: return str(0)

  def ca(self, x):
    if(x is None):           return str(0)
    elif(self.is_bool(x)):   return str(x.count(True))
    elif(self.is_size_t(x)): return str(x.size())
    elif(len(x)==0):         return str(0)
    elif(self.is_size_t(x[0])):
      return str(flex.sum(flex.size_t([i.size() for i in x])))
    else: raise RuntimeError("Bad selection array type.")

  def count_occupancies(self, x):
    result = flex.size_t()
    if(x is not None):
      for i in x:
        for j in i:
          for k in j:
            result.append(k)
    return str(result.size())

  def show(self, log = None):
    if(log is None): log = sys.stdout
    print("Refinement flags and selection counts:", file=log)
    print("  individual_sites       = %5s (%s atoms)"%(
      str(self.individual_sites), self.ca(self.sites_individual)), file=log)
    print("  torsion_angles         = %5s (%s atoms)"%(
      str(self.torsion_angles), self.ca(self.sites_torsion_angles)), file=log)
    print("  rigid_body             = %5s (%s atoms in %s groups)"%(
      str(self.rigid_body), self.ca(self.sites_rigid_body),
      self.szs(self.sites_rigid_body)), file=log)
    print("  individual_adp         = %5s (iso = %s aniso = %s)"%(
      str(self.individual_adp), self.ca(self.adp_individual_iso),
      self.ca(self.adp_individual_aniso)), file=log)
    print("  group_adp              = %5s (%s atoms in %s groups)"%(
      str(self.group_adp), self.ca(self.adp_group), self.szs(self.adp_group)), file=log)
    print("  tls                    = %5s (%s atoms in %s groups)" % (
      str(self.tls), self.ca(self.adp_tls), self.szs(self.adp_tls)), file=log)
    co_res = self.count_occupancies(self.s_occupancies)
    print("  occupancies            = %5s (%s atoms)"%(
      str(self.occupancies), co_res), file=log)
    print("  group_anomalous        = %5s"%self.group_anomalous, file=log) # XXX selections not available
    log.flush()

  def _select(self, x, selection):
    try: lx = len(x)
    except Exception: lx = x.size()
    if(lx == 0): return x
    if(self.is_bool(x)):
      x = x.select(selection)
    elif(self.is_size_t(x[0])):
      x_new = []
      for i_seq, item in enumerate(x):
        val = flex.bool(selection.size(), item).select(selection).iselection()
        if(val.size() > 0): x_new.append(val)
      x = x_new
    elif(self.is_size_t(x[0][0])):
      x_new = []
      for item1 in x:
        tmp = []
        for item2 in item1:
          if(self.is_size_t(item2)):
            v=flex.bool(selection.size(),item2).select(selection).iselection()
            if(v.size() > 0): tmp.append(v)
          else: raise RuntimeError("Bad selection array type.")
        if(len(tmp) > 0): x_new.append(tmp)
      x = x_new
    else: raise RuntimeError("Bad selection array type.")
    return x

  def select_detached(self, selection):
    new_sites_individual     = None
    new_sites_torsion_angles = None
    new_sites_rigid_body     = None
    new_adp_individual_iso   = None
    new_adp_individual_aniso = None
    new_adp_group            = None
    new_group_h              = None
    new_adp_tls              = None
    new_s_occupancies        = None
    if(self.sites_individual is not None):
      new_sites_individual = self._select(self.sites_individual, selection)
    if(self.sites_torsion_angles is not None):
      new_sites_torsion_angles = self._select(
        self.sites_torsion_angles, selection)
    if(self.adp_individual_iso is not None):
      new_adp_individual_iso= self._select(self.adp_individual_iso, selection)
    if(self.adp_individual_aniso is not None):
      new_adp_individual_aniso= self._select(self.adp_individual_aniso,
        selection)
    if(self.sites_rigid_body is not None):
      new_sites_rigid_body = self._select(self.sites_rigid_body, selection)
    if(self.adp_group is not None):
      new_adp_group = self._select(self.adp_group, selection)
    if(self.group_h is not None):
      new_group_h = self._select(self.group_h, selection)
    if(self.adp_tls is not None):
      new_adp_tls = self._select(self.adp_tls, selection)
    if(self.s_occupancies is not None):
      new_s_occupancies = self._select(self.s_occupancies, selection)
    return manager(
      individual_sites       = self.individual_sites,
      torsion_angles         = self.torsion_angles  ,
      rigid_body             = self.rigid_body      ,
      individual_adp         = self.individual_adp  ,
      group_adp              = self.group_adp       ,
      tls                    = self.tls             ,
      occupancies            = self.occupancies     ,
      group_anomalous        = self.group_anomalous ,
      sites_individual       = new_sites_individual    ,
      sites_torsion_angles   = new_sites_torsion_angles,
      sites_rigid_body       = new_sites_rigid_body    ,
      adp_individual_iso     = new_adp_individual_iso  ,
      adp_individual_aniso   = new_adp_individual_aniso,
      adp_group              = new_adp_group           ,
      group_h                = new_group_h             ,
      adp_tls                = new_adp_tls             ,
      s_occupancies          = new_s_occupancies       )

  def select(self, selection):
    assert self.is_bool(selection)
    if(self.sites_individual is not None):
      self.sites_individual = self._select(self.sites_individual, selection)
    if(self.sites_torsion_angles is not None):
      self.sites_torsion_angles = self._select(
        self.sites_torsion_angles, selection)
    if(self.adp_individual_iso is not None):
      self.adp_individual_iso= self._select(self.adp_individual_iso, selection)
    if(self.adp_individual_aniso is not None):
      self.adp_individual_aniso= self._select(self.adp_individual_aniso,
        selection)
    if(self.sites_rigid_body is not None):
      self.sites_rigid_body = self._select(self.sites_rigid_body, selection)
    if(self.adp_group is not None):
      self.adp_group = self._select(self.adp_group, selection)
    if(self.group_h is not None):
      self.group_h = self._select(self.group_h, selection)
    if(self.adp_tls is not None):
      self.adp_tls = self._select(self.adp_tls, selection)
    if(self.s_occupancies is not None):
      self.s_occupancies = self._select(self.s_occupancies, selection)
    return self

  def inflate(self, sites_individual       = None,
                    sites_torsion_angles   = None,
                    sites_rigid_body       = None,
                    adp_individual_iso     = None,
                    adp_individual_aniso   = None,
                    adp_group              = None,
                    group_h                = None,
                    adp_tls                = None,
                    s_occupancies          = None,
                    occupancies_group      = None,
                    size_all               = None):
                    # XXX group_anomalous selection should be added
    if(sites_individual is not None and self.sites_individual is not None):
      assert self.is_bool(sites_individual)
      self.sites_individual.extend(sites_individual)
    if(    sites_torsion_angles is not None
       and self.sites_torsion_angles is not None):
      assert self.is_bool(sites_torsion_angles)
      self.sites_torsion_angles.extend(sites_torsion_angles)
    if(adp_individual_iso is not None):
      assert self.is_bool(adp_individual_iso)
      if(self.adp_individual_iso is None):
        assert size_all is not None
        self.adp_individual_iso = flex.bool(size_all, False)
      # inflate existing iso flags if present
      elif (size_all and self.adp_individual_iso.size() < size_all):
        n_new = size_all - self.adp_individual_iso.size()
        self.adp_individual_iso.extend(flex.bool(n_new, False))
      self.adp_individual_iso.extend(adp_individual_iso)
    if(adp_individual_aniso is not None):
      assert self.is_bool(adp_individual_aniso)
      if(self.adp_individual_aniso is None):
        assert size_all is not None
        self.adp_individual_aniso = flex.bool(size_all, False)
      # inflate existing aniso flags if present
      elif (size_all and self.adp_individual_aniso.size() < size_all):
        n_new = size_all - self.adp_individual_aniso.size()
        self.adp_individual_aniso.extend(flex.bool(n_new, False))
      self.adp_individual_aniso.extend(adp_individual_aniso)
    if(sites_rigid_body is not None):
      assert hasattr(sites_rigid_body, 'count')
      self.sites_rigid_body.extend(sites_rigid_body)
    if(adp_group is not None):
      assert hasattr(adp_group, 'count')
      self.adp_group.extend(adp_group)
    if(group_h is not None):
      assert hasattr(group_h, 'count')
      self.group_h.extend(group_h)
    if(adp_tls is not None):
      assert hasattr(adp_tls, 'count')
      self.adp_tls.extend(adp_tls)
    if(s_occupancies is not None):
      assert hasattr(s_occupancies, 'count')
      if(self.s_occupancies is not None):
        self.s_occupancies.extend(s_occupancies)
    self.check_all()
    return self

  def _add_to_single_size_t(self, x, next_to_i_seq, squeeze_in, mode=1):
    tmp_x = []
    new_independent_element = None
    added = 0
    for sel in x:
      if(sel < next_to_i_seq):
        tmp_x.append(sel)
      elif(sel == next_to_i_seq):
        tmp_x.append(sel)
        if(squeeze_in):
          if(mode == 1):
            if(x.size() > 1):
              tmp_x.append(next_to_i_seq+1)
              added = 1
            elif(x.size() == 1):
              new_independent_element = [next_to_i_seq+1]
              added = 1
            else: raise RuntimeError
          elif(mode == 2):
            new_independent_element = [next_to_i_seq+1]
            added = 1
          else: raise RuntimeError
      else:
        tmp_x.append(sel+1)
    if(new_independent_element is not None):
      new_independent_element = flex.size_t(new_independent_element)
    return flex.size_t(tmp_x), new_independent_element, added

  def _add(self, x, next_to_i_seq, squeeze_in):
    if(x is None): return x
    elif(self.is_bool(x)):
      x_new = []
      result = self._add_to_single_size_t(x.iselection(),
        next_to_i_seq, squeeze_in)
      if(result[1] is None): x_new.extend([result[0]])
      else: x_new.extend([result[0], result[1]])
      if(result[2] == 0 and squeeze_in):
        x_new.extend([flex.size_t([next_to_i_seq+1])])
      x_new_new = flex.size_t()
      for i_x_new in x_new: x_new_new.extend(i_x_new)
      return flex.bool(x.size()+1, x_new_new)
    #elif(self.is_size_t(x)): # XXX not testes so disabled
    #  return self._add_to_single_size_t(x, next_to_i_seq, squeeze_in)
    elif(len(x)==0): raise RuntimeError
    elif(self.is_size_t(x[0])):
      x_new = []
      added = 0
      for x_ in x:
        result = self._add_to_single_size_t(x_, next_to_i_seq, squeeze_in)
        added += result[2]
        if(result[1] is None): x_new.extend([result[0]])
        else: x_new.extend([result[0], result[1]])
      if(added == 0 and squeeze_in):
        x_new.extend([flex.size_t([next_to_i_seq+1])])
      return x_new
    elif(self.is_size_t(x[0][0])):
      xx_new = []
      added = 0
      result_1 = None
      for xx in x:
        x_new = []
        for x_ in xx:
          result = self._add_to_single_size_t(x_, next_to_i_seq, squeeze_in, 2)
          added += result[2]
          if(result[1] is None):
            x_new.extend([result[0]])
          else:
            x_new.extend([result[0]])
          if(result[1] is not None):
            xx_new.append([result[1]])
        if(len(x_new) > 0): xx_new.append(x_new)
      if(added == 0 and squeeze_in):
        xx_new.append([flex.size_t([next_to_i_seq+1])])
      return xx_new
    else: raise RuntimeError("Bad selection array type.")

  def add(self, next_to_i_seqs,
                sites_individual     = False,
                sites_torsion_angles = False,
                sites_rigid_body     = False,
                adp_individual_iso   = False,
                adp_individual_aniso = False,
                adp_group            = False,
                group_h              = False,
                adp_tls              = False,
                s_occupancies        = False):
                # XXX group_anomalous selection should be added
    next_to_i_seqs = flex.size_t(next_to_i_seqs)
    perm = flex.sort_permutation(next_to_i_seqs, reverse = True)
    next_to_i_seqs = next_to_i_seqs.select(perm)
    for next_to_i_seq in next_to_i_seqs:
      if(self.sites_individual is not None):
        self.sites_individual = self._add(
          x             = self.sites_individual,
          next_to_i_seq = next_to_i_seq,
          squeeze_in    = sites_individual)
      if(self.sites_torsion_angles is not None):
        self.sites_torsion_angles = self._add(
          x             = self.sites_torsion_angles,
          next_to_i_seq = next_to_i_seq,
          squeeze_in    = sites_torsion_angles)
      if(self.sites_rigid_body is not None):
        self.sites_rigid_body = self._add(
          x             = self.sites_rigid_body,
          next_to_i_seq = next_to_i_seq,
          squeeze_in    = sites_rigid_body)
      if(self.adp_individual_iso is not None):
        self.adp_individual_iso = self._add(
          x             = self.adp_individual_iso,
          next_to_i_seq = next_to_i_seq,
          squeeze_in    = adp_individual_iso)
      if(self.adp_individual_aniso is not None):
        self.adp_individual_aniso = self._add(
          x             = self.adp_individual_aniso,
          next_to_i_seq = next_to_i_seq,
          squeeze_in    = adp_individual_aniso)
      if(self.adp_group is not None):
        self.adp_group = self._add(
          x             = self.adp_group,
          next_to_i_seq = next_to_i_seq,
          squeeze_in    = adp_group)
      if(self.group_h is not None):
        self.group_h = self._add(
          x             = self.group_h,
          next_to_i_seq = next_to_i_seq,
          squeeze_in    = group_h)
      if(self.adp_tls is not None):
        self.adp_tls = self._add(
          x             = self.adp_tls,
          next_to_i_seq = next_to_i_seq,
          squeeze_in    = adp_tls)
      if(self.s_occupancies is not None):
        self.s_occupancies = self._add(
          x             = self.s_occupancies,
          next_to_i_seq = next_to_i_seq,
          squeeze_in    = s_occupancies)
    return self
