
from __future__ import absolute_import, division, print_function
from gltbx.gl import *
import scitbx.math
from scitbx.array_family import flex
import scitbx.matrix
import time
import sys
from six.moves import range

class cartoon (object) :
  def __init__ (self, pdb_hierarchy, sec_str, selection_cache=None) :
    self._segments = []
    if selection_cache is None :
      selection_cache = pdb_hierarchy.atom_selection_cache()
    last_atom     = None
    ca_sele = selection_cache.selection("name ' CA '")
    c_sele = selection_cache.selection("name ' C  '")
    o_sele = selection_cache.selection("name ' O  '")
    last_i_seq = None
    last_labels = None
    last_site = None
    last_ss = None
    last_resseq = - sys.maxsize
    current_segment = None
    t1 = time.time()
    for model in pdb_hierarchy.models() :
      for chain in pdb_hierarchy.chains() :
        main_conf = chain.conformers()[0]
        for residue in main_conf.residues() :
          resseq = residue.resseq_as_int()
          for atom in residue.atoms() :
            i_seq = atom.i_seq
            site_cart = atom.xyz
            if ca_sele[i_seq] :
              ss_type = sec_str[i_seq]
              if (resseq > (last_resseq + 1)) :
                current_segment = self.new_segment(site_cart=site_cart,
                  i_seq=i_seq,
                  ss_type=ss_type)
              elif (ss_type != last_ss) :
                if (current_segment is not None) :
                  current_segment.set_next_site(site_cart)
                current_segment = self.new_segment(site_cart=site_cart,
                  i_seq=i_seq,
                  ss_type=ss_type,
                  prev_site=last_site)
              else :
                current_segment.add_residue(site_cart, i_seq)
              last_ss = ss_type
              last_resseq = resseq
              last_site = site_cart
              break
    print("%d segments" % len(self._segments))
    t2 = time.time()
    print("Extract backbone: %.3fms" % ((t2-t1) * 1000))

  def new_segment (self, site_cart, i_seq, ss_type, prev_site=None) :
    new_seg = None
    if (ss_type == 0) :
      new_seg = loop(site_cart, i_seq, prev_site)
    elif (ss_type == 1) :
      new_seg = helix(site_cart, i_seq, prev_site)
    elif (ss_type == 2) :
      new_seg = strand(site_cart, i_seq, prev_site)
    assert (new_seg is not None)
    self._segments.append(new_seg)
    return new_seg

  def construct_geometry (self, smoothness=5) :
    t1 = time.time()
    for segment in self._segments :
      segment.construct_geometry(smoothness)
    t2 = time.time()
    print("Construct geometry: %.3fms" % ((t2-t1) * 1000))

  def draw_ribbon (self, atom_colors, atoms_visible) :
    for segment in self._segments :
      segment.draw_ribbon(atom_colors, atoms_visible)

  def get_points_and_lines (self) :
    points = flex.vec3_double()
    line_i_seqs = []
    k = 0
    for segment in self._segments :
      vertices = segment.get_vertices()
      points.extend(vertices)
      for i in range(k, k + vertices.size() - 1) :
        line_i_seqs.append((i, i+1))
      k += vertices.size()
    return (points, line_i_seqs)

class segment (object) :
  def __init__ (self, site_cart, i_seq, prev_site=None) :
    self._prev_site = prev_site
    self._next_site = None
    self._n_sites = 1
    self._anchors = flex.vec3_double()
    self._anchors.append(site_cart)
    self._indices = flex.size_t()
    self._indices.append(i_seq)
    self._vertices = None
    self._vertex_i_seqs = None

  def add_residue (self, site_cart, i_seq) :
    self._anchors.append(site_cart)
    self._indices.append(i_seq)

  def set_next_site (self, next_site) :
    self._next_site = next_site

  def get_vertices (self) :
    if (self._vertices is None) :
      self.construct_geometry()
    return self._vertices

  def construct_geometry (self, smoothness=5) :
    n_sites = len(self._anchors)
    anchors = self._anchors
    indices = self._indices
    if (self._prev_site is not None) :
      anchors.insert(0, self._prev_site)
      indices.insert(0, indices[0])
    else :
      v01 = vec3(anchors[0]) - vec3(anchors[1])
      v00 = vec3(anchors[0]) + v01
      anchors.insert(0, v00.elems)
      indices.insert(0, indices[0])
    if (self._next_site is not None) :
      anchors.append(self._next_site)
      indices.append(indices[-1])
    else :
      vlast = vec3(anchors[-2]) - vec3(anchors[-1])
      vXX = vec3(anchors[-1]) - vlast
      anchors.append(vXX.elems)
      indices.append(indices[-1])
    vertices = flex.vec3_double()
    vertex_i_seqs = flex.size_t()
    for i in range(1, n_sites) :
      v0 = vec3(anchors[i-1])
      v1 = vec3(anchors[i])
      v2 = vec3(anchors[i+1])
      v3 = vec3(anchors[i+2])
      v10 = (v1 - v0).normalize()
      v23 = (v2 - v3).normalize()
      v05 = (v0 + v1) / 2.0
      _v15 = (v1 + v2) / 2.0
      v25 = (v2 + v3) / 2.0
      v15 = _v15 + ((v10 + v23) / 2.0)
      if (i == 1) :
        vertices.append(v05.elems)
        vertex_i_seqs.append(indices[0])
        new_vertices_0 = crspline_interp(v0.elems, v05.elems, v1.elems, v15.elems, smoothness)
        vertices.extend(new_vertices_0)
        for k in range(new_vertices_0.size()) :
          vertex_i_seqs.append(indices[i])
      vertices.append(anchors[i])
      vertex_i_seqs.append(indices[i])
      new_vertices_1 = crspline_interp(v05.elems, v1.elems, v15.elems, v2.elems,        smoothness)
      vertices.extend(new_vertices_1)
      for k in range(new_vertices_1.size()) :
        vertex_i_seqs.append(indices[i])
      vertices.append(v15.elems)
      vertex_i_seqs.append(indices[i])
      new_vertices_2 = crspline_interp(v1.elems, v15.elems, v2.elems, v25.elems,
        smoothness)
      vertices.extend(new_vertices_2)
      for k in range(new_vertices_2.size()) :
        vertex_i_seqs.append(indices[i+1])
      if (i == (n_sites - 1)) :
        new_vertices_3 = crspline_interp(v15.elems, v2.elems, v25.elems, v3.elems, smoothness)
        vertices.extend(new_vertices_3)
        for k in range(new_vertices_3.size()) :
          vertex_i_seqs.append(indices[n_sites])
    assert (vertices.size() == vertex_i_seqs.size())
    self._vertices = vertices
    self._vertex_i_seqs = vertex_i_seqs

  def draw_ribbon (self, atom_colors, atoms_visible) :
    vertices = self._vertices
    indices = self._vertex_i_seqs
    in_line_loop = False
    glLineWidth(2.0)
    for k, point in enumerate(vertices) :
      i_seq = indices[k]
      if atoms_visible[i_seq] :
        if (not in_line_loop) :
          glBegin(GL_LINE_STRIP)
          in_line_loop = True
        glColor3f(*atom_colors[i_seq])
        glVertex3f(*point)
      elif in_line_loop :
        glEnd()
        in_line_loop = False
    if in_line_loop :
      glEnd()

class loop (segment) :
  pass

class helix (segment) :
  pass

class strand (segment) :
  pass

def vec3 (xyz) :
  return scitbx.matrix.rec(xyz, (1,3))

def crspline_interp (*args) :
  return scitbx.math.interpolate_catmull_rom_spline(*args)
