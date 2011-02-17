from cctbx.array_family import flex # import dependency
from cctbx.eltbx import wavelengths, tiny_pse
from cctbx import adptbx

def generator(xray_structure,
              data_are_intensities=True,
              title=None,
              wavelength=None,
              temperature=None,
              least_square_cyles=None,
              overall_scale_factor=None,
              weighting_scheme_params=None,
              sort_scatterers = True,
              ):
  space_group = xray_structure.space_group()
  assert not space_group.is_centric() or space_group.is_origin_centric()
  if title is None:
    title = '????'
  if wavelength is None:
    wavelength = wavelengths.characteristic('Mo').as_angstrom()
  sgi = xray_structure.space_group_info()
  uc = xray_structure.unit_cell()

  yield 'TITL %s in %s\n' % (title, sgi.type().lookup_symbol())
  yield 'CELL %.5f %s\n' % (
    wavelength,
    ' '.join(('%.4f ',)*3 + ('%.3f',)*3) % uc.parameters())
  yield 'ZERR %i 0. 0. 0. 0. 0. 0.\n' % sgi.group().order_z()

  latt = 1 + 'PIRFABC'.find(sgi.group().conventional_centring_type_symbol())
  if space_group.is_origin_centric(): latt = -latt
  yield 'LATT %i\n' % latt
  for i in xrange(space_group.n_smx()):
    rt_mx = space_group(0, 0, i)
    yield 'SYMM %s\n' % rt_mx
  yield '\n'

  uc_content = xray_structure.unit_cell_content()
  for e in uc_content:
    uc_content[e] = "%.1f" % uc_content[e]
  sfac = []
  unit = []
  prior = ('C', 'H')
  for e in prior:
    if e in uc_content:
      sfac.append(e)
      unit.append(uc_content[e])
  dsu = [ (tiny_pse.table(e).atomic_number(), e) for e in uc_content ]
  dsu.sort()
  sorted = [ item[-1] for item in dsu ]
  for e in sorted:
    sfac.append(e)
    unit.append(uc_content[e])
  yield 'SFAC %s\n' % ' '.join(sfac)
  yield 'UNIT %s\n' % ' '.join(unit)
  sf_idx = dict([ (e, i + 1) for i, e in enumerate(sfac) ])
  yield '\n'

  if temperature:
    yield 'TEMP %.0f\n' % temperature

  if least_square_cyles:
    yield 'L.S. %i\n' % least_square_cyles

  yield '\n'

  if weighting_scheme_params is not None:
    a, b = weighting_scheme_params
    if b is None:
      yield 'WGHT %.6f\n' % a
    else:
      yield 'WGHT %.6f %.6f\n' % (a, b)

  if overall_scale_factor is not None:
    yield 'FVAR %.5f' % overall_scale_factor

  fmt_tmpl = ('%-4s', '%-2i') + ('%.6f',)*3 + ('%.5f',)
  fmt_iso = ' '.join(fmt_tmpl + ('%.5f', '\n'))
  fmt_aniso = ' '.join(
    fmt_tmpl + ('%.5f',)*2 + ('=\n ',) + ('%.5f',)*4 + ('\n',))
  if sort_scatterers:
    dsu = [ (tiny_pse.table(sc.scattering_type).atomic_number(), sc)
            for sc in xray_structure.scatterers() ]
    dsu.sort(reverse=True)
    scatterers = flex.xray_scatterer([ item[-1] for item in dsu ])
  else:
    scatterers = xray_structure.scatterers()
  for sc in scatterers:
    assert sc.flags.use_u_iso() ^ sc.flags.use_u_aniso()
    params = (sc.label, sf_idx[sc.scattering_type]) + sc.site
    occ = sc.occupancy
    if not sc.flags.grad_occupancy(): occ += 10
    params += (occ, )
    if sc.flags.use_u_iso():
      yield fmt_iso %  (params + (sc.u_iso,))
    else:
      u11, u22, u33, u12, u13, u23 = adptbx.u_star_as_u_cif(uc, sc.u_star)
      yield fmt_aniso % (params + (u11, u22, u33, u23, u13, u12))

  if data_are_intensities: hklf = 4
  else: hklf = 3
  yield 'HKLF %i\n' % hklf
