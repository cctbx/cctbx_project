from __future__ import absolute_import, division, print_function

import smtbx.refinement
from libtbx.test_utils import approx_equal
from smtbx.regression.test_data import fnames

expected_str_reparametrisation = """\
digraph dependencies {
204 -> 0;
205 -> 1;
0 [label="independent_occupancy_parameter (N3) #0"];
1 [label="independent_occupancy_parameter (C7A) #1"];
2 [label="independent_site_parameter (F1) #2"];
5 [label="independent_u_star_parameter (F1) #5"];
11 [label="independent_site_parameter (F2) #11"];
14 [label="independent_u_star_parameter (F2) #14"];
20 [label="independent_site_parameter (N8) #20"];
23 [label="independent_u_star_parameter (N8) #23"];
29 [label="independent_site_parameter (N3) #29"];
32 [label="independent_u_iso_parameter (N3) #32"];
33 [label="independent_site_parameter (C9) #33"];
36 [label="independent_u_star_parameter (C9) #36"];
42 [label="independent_site_parameter (C4) #42"];
45 [label="independent_u_star_parameter (C4) #45"];
51 [label="independent_site_parameter (N5) #51"];
54 [label="independent_u_star_parameter (N5) #54"];
60 [label="independent_site_parameter (C2) #60"];
63 [label="independent_u_star_parameter (C2) #63"];
69 [label="independent_site_parameter (C10) #69"];
72 [label="independent_u_star_parameter (C10) #72"];
78 [label="independent_site_parameter (C1) #78"];
81 [label="independent_u_star_parameter (C1) #81"];
87 [label="independent_site_parameter (C11) #87"];
90 [label="independent_u_star_parameter (C11) #90"];
96 [label="independent_site_parameter (C13) #96"];
99 [label="independent_u_star_parameter (C13) #99"];
105 [label="independent_site_parameter (C6) #105"];
108 [label="independent_u_star_parameter (C6) #108"];
114 [label="independent_site_parameter (N12) #114"];
117 [label="independent_u_star_parameter (N12) #117"];
123 [label="independent_site_parameter (C7A) #123"];
126 [label="independent_u_star_parameter (C7A) #126"];
132 [label="independent_site_parameter (C14) #132"];
135 [label="independent_u_star_parameter (C14) #135"];
141 [label="independent_site_parameter (C7B) #141"];
144 [label="independent_u_star_parameter (C7B) #144"];
150 [label="independent_site_parameter (C3) #150"];
153 [label="independent_u_iso_parameter (C3) #153"];
154 [label="independent_occupancy_parameter [cst] (F1) #154"];
155 [label="independent_fp_parameter [cst] (F1) #155"];
156 [label="independent_fdp_parameter [cst] (F1) #156"];
157 [label="independent_occupancy_parameter [cst] (F2) #157"];
158 [label="independent_fp_parameter [cst] (F2) #158"];
159 [label="independent_fdp_parameter [cst] (F2) #159"];
160 [label="independent_occupancy_parameter [cst] (N8) #160"];
161 [label="independent_fp_parameter [cst] (N8) #161"];
162 [label="independent_fdp_parameter [cst] (N8) #162"];
163 [label="independent_fp_parameter [cst] (N3) #163"];
164 [label="independent_fdp_parameter [cst] (N3) #164"];
165 [label="independent_occupancy_parameter [cst] (C9) #165"];
166 [label="independent_fp_parameter [cst] (C9) #166"];
167 [label="independent_fdp_parameter [cst] (C9) #167"];
168 [label="independent_occupancy_parameter [cst] (C4) #168"];
169 [label="independent_fp_parameter [cst] (C4) #169"];
170 [label="independent_fdp_parameter [cst] (C4) #170"];
171 [label="independent_occupancy_parameter [cst] (N5) #171"];
172 [label="independent_fp_parameter [cst] (N5) #172"];
173 [label="independent_fdp_parameter [cst] (N5) #173"];
174 [label="independent_occupancy_parameter [cst] (C2) #174"];
175 [label="independent_fp_parameter [cst] (C2) #175"];
176 [label="independent_fdp_parameter [cst] (C2) #176"];
177 [label="independent_occupancy_parameter [cst] (C10) #177"];
178 [label="independent_fp_parameter [cst] (C10) #178"];
179 [label="independent_fdp_parameter [cst] (C10) #179"];
180 [label="independent_occupancy_parameter [cst] (C1) #180"];
181 [label="independent_fp_parameter [cst] (C1) #181"];
182 [label="independent_fdp_parameter [cst] (C1) #182"];
183 [label="independent_occupancy_parameter [cst] (C11) #183"];
184 [label="independent_fp_parameter [cst] (C11) #184"];
185 [label="independent_fdp_parameter [cst] (C11) #185"];
186 [label="independent_occupancy_parameter [cst] (C13) #186"];
187 [label="independent_fp_parameter [cst] (C13) #187"];
188 [label="independent_fdp_parameter [cst] (C13) #188"];
189 [label="independent_occupancy_parameter [cst] (C6) #189"];
190 [label="independent_fp_parameter [cst] (C6) #190"];
191 [label="independent_fdp_parameter [cst] (C6) #191"];
192 [label="independent_occupancy_parameter [cst] (N12) #192"];
193 [label="independent_fp_parameter [cst] (N12) #193"];
194 [label="independent_fdp_parameter [cst] (N12) #194"];
195 [label="independent_fp_parameter [cst] (C7A) #195"];
196 [label="independent_fdp_parameter [cst] (C7A) #196"];
197 [label="independent_occupancy_parameter [cst] (C14) #197"];
198 [label="independent_fp_parameter [cst] (C14) #198"];
199 [label="independent_fdp_parameter [cst] (C14) #199"];
200 [label="independent_fp_parameter [cst] (C7B) #200"];
201 [label="independent_fdp_parameter [cst] (C7B) #201"];
202 [label="independent_fp_parameter [cst] (C3) #202"];
203 [label="independent_fdp_parameter [cst] (C3) #203"];
204 [label="affine_asu_occupancy_parameter (C3) #204"];
205 [label="affine_asu_occupancy_parameter (C7B) #205"]
}
"""

def exercise_simple_disorder():
  ins = fnames.thpp_ins
  model = smtbx.refinement.model.from_shelx(ins)
  ls = model.least_squares()
  assert str(ls.reparametrisation).strip() == \
          expected_str_reparametrisation.strip()
  ls.build_up()
  covann = ls.covariance_matrix_and_annotations()
  assert approx_equal(covann.variance_of('C7B.occ'),
                      covann.variance_of('C7A.occ'))
  assert approx_equal(covann.variance_of('C3.occ'),
                      covann.variance_of('N3.occ'))

if __name__ == '__main__':
  exercise_simple_disorder()
  print('OK')
