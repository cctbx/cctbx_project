from __future__ import absolute_import, division, print_function
import math
from cctbx import adptbx
import libtbx
import os
import json
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

import time

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

def load_table(element):
  element = element.strip().upper()
  assert element in ["C","H","N","O","S"]
  path=libtbx.env.find_in_repositories("cctbx/maptbx/bcr/tables")
  file_name = "%s/%s.json"%(path, element)
  assert os.path.isfile(file_name)
  with open(file_name, 'r') as file:
    return json.load(file)

def compute(hierarchy, unit_cell, n_real, resolution=None, resolutions=None,
            debug=False):
  assert [resolution, resolutions].count(None)==1
  atoms = hierarchy.atoms()
  if resolutions is None:
    resolutions = [resolution,] * atoms.size()
  RadFact = 2.0
  RadAdd  = 0.5
  RadMu   = RadFact + RadAdd
  sites_cart = atoms.extract_xyz()
  adp_as_u = atoms.extract_b()*adptbx.b_as_u(1.)
  occupancy= atoms.extract_occ()
  arrays = {}
  for e in list(set(atoms.extract_element())):
    e = e.strip().upper()
    d = load_table(element=e)
    arrays[e] = d
  ScaleB = 1.0 / (8.0 * math.pi**2)
  kscale = math.pi**1.5
  bcr_scatterers = []
  for s, u, o, e, r in zip(sites_cart, adp_as_u, occupancy,
                           atoms.extract_element(), resolutions):
    e = e.strip().upper()
    entry = arrays[e]
    keys = [float(x) for x in entry.keys()]
    key = str(min(keys, key=lambda x: abs(x - r)))
    vals = entry[key]
    R = flex.double(vals['R'])
    B = flex.double(vals['B'])
    C = flex.double(vals['C'])
    sel = R < (r*RadMu)
    mu    = R.select(sel)
    nu    = B.select(sel) * ScaleB
    kappa = C.select(sel)
    musq  = mu * mu
    kappi = kappa/kscale
    bcr_scatterer = ext.bcr_scatterer(
      site_cart = s,
      u_iso     = u,
      occ       = o,
      radius    = r*RadFact, # atomic radius = atomic_resolution * RadFact
      resolution=r,
      mu        = mu,
      kappa     = kappa,
      nu        = nu,
      musq      = musq,
      kappi     = kappi)
    bcr_scatterers.append(bcr_scatterer)
  #
  if debug:
    OmegaMap_py = CalcOmegaMap(Ncrs=n_real, Scrs=[0,0,0], Nxyz=n_real,
      unit_cell=unit_cell, bcr_scatterers=bcr_scatterers)
  #
  o = ext.vrm(
    Ncrs           = n_real,
    Scrs           = [0,0,0],
    Nxyz           = n_real,
    unit_cell      = unit_cell,
    bcr_scatterers = bcr_scatterers)
  #
  #t0 = time.time()
  OmegaMap_cpp = o.compute_map(arg_value=False)
  #print("cpp: ", time.time()-t0)

  OmegaMap_cpp_2 = o.map  # alternative way to get the map
  # Re-order
  nx,ny,nz = n_real
  mpy  = flex.double(flex.grid(n_real))
  mcpp = flex.double(flex.grid(n_real))
  for iz in range(0, nz):
    for iy in range(0, ny):
      for ix in range(0, nx):
        if debug:
          mpy[ix,iy,iz] = OmegaMap_py[iz][iy][ix]
        mcpp[ix,iy,iz] = OmegaMap_cpp[iz,iy,ix]
  #
  if debug:
    cc = flex.linear_correlation(x=mcpp.as_1d(), y=mpy.as_1d()).coefficient()
    assert approx_equal(cc, 1.0)
    cc = flex.linear_correlation(
       x=OmegaMap_cpp.as_1d(), y=OmegaMap_cpp_2.as_1d()).coefficient()
    assert approx_equal(cc, 1.0)
    for func in [flex.min, flex.max, flex.mean]:
      assert approx_equal(func(mcpp), func(mpy))
  #
  if debug: return mcpp, OmegaMap_py, bcr_scatterers
  else:     return mcpp

#-------------------------------------------
def CalcFuncMap(OmegaMap, ControlMap, Ncrs) :

    Mx,My,Mz = Ncrs

    FuncMap = 0.0

    for iz in range(Mz) :
       for iy in range(My) :
          for ix in range(Mx) :
              FuncMap = FuncMap + (OmegaMap[iz][iy][ix] - ControlMap[iz][iy][ix]) **2

    FuncMap = 0.5 * FuncMap

    return FuncMap

#-------------------------------------------
def CalcGradMap(OmegaMap ,ControlMap, Ncrs) :

    Mx,My,Mz = Ncrs

    GradMap = [[[OmegaMap[iz][iy][ix]-ControlMap[iz][iy][ix] \
                 for ix in range(Mx)] for iy in range(My)] for iz in range(Mz)]

    return GradMap

#=====================================
def CalcOmegaMap(Ncrs, Scrs, Nxyz, unit_cell, bcr_scatterers) :

    acell, bcell, ccell, alpha, beta,  gamma = unit_cell.parameters()
    OrthMatrix  = unit_cell.orthogonalization_matrix()
    DeortMatrix = unit_cell.fractionalization_matrix()

    Mx, My, Mz = Ncrs
    Sx, Sy, Sz = Scrs
    Nx, Ny, Nz = Nxyz
    Fx, Fy, Fz = Sx + Mx, Sy + My, Sz + Mz
    StepX, StepY, StepZ = 1.0/float(Nx), 1.0/float(Ny), 1.0/float(Nz)

    orthxx, orthxy, orthxz = OrthMatrix[0] , OrthMatrix[3] , OrthMatrix[6]
    orthyy, orthyz         =                 OrthMatrix[4] , OrthMatrix[7]
    orthzz                 =                                 OrthMatrix[8]

    dortxx, dortxy, dortxz = DeortMatrix[0], DeortMatrix[3], DeortMatrix[6]
    dortyy, dortyz         =                 DeortMatrix[4], DeortMatrix[7]
    dortzz                 =                                 DeortMatrix[8]

    StepXX, StepXY, StepXZ = orthxx * StepX, orthxy * StepY, orthxz * StepZ
    StepYY, StepYZ         =                 orthyy * StepY, orthyz * StepZ
    StepZZ                 =                                 orthzz * StepZ

    StepXX2, StepXXS = StepXX * 2, StepXX * StepXX
    StepYY2, StepYYS = StepYY * 2, StepYY * StepYY
    StepZZ2, StepZZS = StepZZ * 2, StepZZ * StepZZ

    StepXXS2, StepYYS2, StepZZS2 = StepXXS * 2, StepYYS * 2, StepZZS * 2

    cosalp, cosbet, cosgam = math.cos(alpha), math.cos(beta), math.cos(gamma)
    vol0 = math.sqrt(1.0 - cosalp*cosalp - cosbet*cosbet - cosgam*cosgam + 2.0*cosalp*cosbet*cosgam)

    RprojX = math.sin(alpha) / (vol0 * acell)
    RprojY = math.sin(beta)  / (vol0 * bcell)
    RprojZ = math.sin(gamma) / (vol0 * ccell)

    OmegaMap = [[[ 0.0 for ix in range(Mx) ] for iy in range(My)] for iz in range(Mz)]

    for iatom in range(len(bcr_scatterers)) :
        bcr_scatterer = bcr_scatterers[iatom]

        RadAtom = bcr_scatterer.radius
        RadAtom2  = RadAtom   * RadAtom
        RadAtomX  = RadAtom   * RprojX
        RadAtomY  = RadAtom   * RprojY
        RadAtomZ  = RadAtom   * RprojZ

        r = bcr_scatterer.site_cart
        xat, yat, zat = r[0], r[1], r[2]
        cat, bat = bcr_scatterer.occ, bcr_scatterer.u_iso

#       get parameters of the box around the atom

        xfrac = dortxx * xat + dortxy * yat + dortxz * zat
        yfrac =                dortyy * yat + dortyz * zat
        zfrac =                               dortzz * zat

        Kx1,Kx2,Ky1,Ky2,Kz1,Kz2,dxf,dyf,dzf  = AtomBox(xfrac, yfrac, zfrac,
                   RadAtomX,RadAtomY,RadAtomZ,StepX,StepY,StepZ,Sx,Sy,Sz,Fx,Fy,Fz)
        dx0 = orthxx * dxf + orthxy * dyf + orthxz * dzf
        dy0 =                orthyy * dyf + orthyz * dzf
        dz0 =                               orthzz * dzf

#       identify grid points to which the atom contributes

        AtomGrids = []
        AtomGrid0 = [-1, -1, -1]
        ix0       = -1

        dz  = dz0
        r2z = dz * dz
        cz  = dz * StepZZ2 + StepZZS
        for iz in range(Kz1,Kz2) :
            dy   = dy0
            dt0  = dx0
            r2zy = dy * dy + r2z
            cy   = dy * StepYY2 + StepYYS

            for iy in range(Ky1,Ky2) :
                dx    = dt0
                r2zyx = dx * dx + r2zy
                cx    = dx * StepXX2 + StepXXS

                for ix in range(Kx1,Kx2) :
                    if r2zyx == 0.0 :
                       AtomGrid0 = [ix,iy,iz]
                       ix0       = ix
                    elif r2zyx <= RadAtom2 :
                       rzyx  = math.sqrt(r2zyx)
                       rzyx2 = rzyx * 2
                       AtomGrids.append([r2zyx,rzyx,rzyx2,ix,iy,iz])

                    r2zyx = r2zyx + cx
                    cx    = cx    + StepXXS2
                    dx    = dx    + StepXX

                r2zy = r2zy + cy
                cy   = cy   + StepYYS2
                dt0  = dt0  + StepXY
                dy   = dy   + StepYY

            r2z = r2z + cz
            cz  = cz  + StepZZS2
            dx0 = dx0 + StepXZ
            dy0 = dy0 + StepYZ
            dz  = dz  + StepZZ

        dz  = dz0
        r2z = dz * dz
        cz  = dz * StepZZ2 + StepZZS

#       calculte contributions

        Ngrids     = len(AtomGrids)
        GridValues = [0.0 for i in range(Ngrids)]
        GridValue0 = 0.0

        for iterm in range(0, bcr_scatterer.mu.size()) :
            mu    = bcr_scatterer.mu[iterm]
            nu    = bcr_scatterer.nu[iterm]
            musq  = bcr_scatterer.musq[iterm]
            kappi = bcr_scatterer.kappi[iterm]

            nuatom  = nu + bat
            nuatom2 = nuatom + nuatom
            fact1   = kappi / nuatom2**1.5
            munuat  = mu / nuatom

#           contribution to the node coinciding with the atomic center

            if ix0 >= 0 :
               argg  = musq / nuatom2
               fact2 = math.exp(-argg)
               GridValue0 = GridValue0 + fact1 * fact2

#           contribution to nodes different from the atomic center

            if mu == 0.0 :
               for ig in range(Ngrids) :
                   r2zyx,rzyx,rzyx2,ix,iy,iz = AtomGrids[ig]
                   argg  = r2zyx / nuatom2
                   fact2 = math.exp(-argg)
                   GridValues[ig] = GridValues[ig] + fact1 * fact2

            else :
               for ig in range(Ngrids) :
                   r2zyx,rzyx,rzyx2,ix,iy,iz = AtomGrids[ig]
                   tterm   = rzyx2 * munuat
                   argg  = (rzyx-mu)**2 / nuatom2
                   fact2 = math.exp(-argg) * (1.0 - math.exp(-tterm)) / tterm
                   GridValues[ig] = GridValues[ig] + fact1 * fact2
#       unpack values array

        for ig in range(Ngrids) :
            r2zyx,rzyx,rzyx2,ix,iy,iz = AtomGrids[ig]
            OmegaMap[iz][iy][ix] += GridValues[ig] * cat

        ix0,iy0,iz0 = AtomGrid0
        if ix0 >= 0 :
           OmegaMap[iz0][iy0][ix0] += GridValue0 * cat

    return OmegaMap

#=====================================
def CalcGradAtom(GradMap, Ncrs, Scrs, Nxyz, unit_cell, bcr_scatterers):

    Natoms = len(bcr_scatterers)

    acell, bcell, ccell, alpha, beta,  gamma = unit_cell.parameters()
    OrthMatrix  = unit_cell.orthogonalization_matrix()
    DeortMatrix = unit_cell.fractionalization_matrix()

    Mx, My, Mz = Ncrs
    Sx, Sy, Sz = Scrs
    Nx, Ny, Nz = Nxyz
    Fx, Fy, Fz = Sx + Mx, Sy + My, Sz + Mz
    StepX, StepY, StepZ = 1.0/float(Nx), 1.0/float(Ny), 1.0/float(Nz)

    orthxx, orthxy, orthxz = OrthMatrix[0] , OrthMatrix[3] , OrthMatrix[6]
    orthyy, orthyz         =                 OrthMatrix[4] , OrthMatrix[7]
    orthzz                 =                                 OrthMatrix[8]

    dortxx, dortxy, dortxz = DeortMatrix[0], DeortMatrix[3], DeortMatrix[6]
    dortyy, dortyz         =                 DeortMatrix[4], DeortMatrix[7]
    dortzz                 =                                 DeortMatrix[8]

    StepXX, StepXY, StepXZ = orthxx * StepX, orthxy * StepY, orthxz * StepZ
    StepYY, StepYZ         =                 orthyy * StepY, orthyz * StepZ
    StepZZ                 =                                 orthzz * StepZ

    StepXX2, StepXXS = StepXX * 2, StepXX * StepXX
    StepYY2, StepYYS = StepYY * 2, StepYY * StepYY
    StepZZ2, StepZZS = StepZZ * 2, StepZZ * StepZZ

    StepXXS2, StepYYS2, StepZZS2 = StepXXS * 2, StepYYS * 2, StepZZS * 2

    cosalp, cosbet, cosgam = math.cos(alpha), math.cos(beta), math.cos(gamma)
    vol0 = math.sqrt(1.0 - cosalp*cosalp - cosbet*cosbet - cosgam*cosgam + 2.0*cosalp*cosbet*cosgam)

    RprojX = math.sin(alpha) / (vol0 * acell)
    RprojY = math.sin(beta)  / (vol0 * bcell)
    RprojZ = math.sin(gamma) / (vol0 * ccell)

    GradAtom = [[0.0, 0.0, 0.0, 0.0, 0.0] for iat in range(Natoms)]

    for iatom in range(len(bcr_scatterers)) :
        bcr_scatterer = bcr_scatterers[iatom]

        RadAtom = bcr_scatterer.radius
        RadAtom2  = RadAtom   * RadAtom
        RadAtomX  = RadAtom   * RprojX
        RadAtomY  = RadAtom   * RprojY
        RadAtomZ  = RadAtom   * RprojZ

        r = bcr_scatterer.site_cart
        xat, yat, zat = r[0], r[1], r[2]
        cat, bat = bcr_scatterer.occ, bcr_scatterer.u_iso

#       get parameters of the box around the atom

        xfrac = dortxx * xat + dortxy * yat + dortxz * zat
        yfrac =                dortyy * yat + dortyz * zat
        zfrac =                               dortzz * zat

        Kx1,Kx2,Ky1,Ky2,Kz1,Kz2,dxf,dyf,dzf  = AtomBox(xfrac, yfrac, zfrac,
                   RadAtomX,RadAtomY,RadAtomZ,StepX,StepY,StepZ,Sx,Sy,Sz,Fx,Fy,Fz)

        dx0 = orthxx * dxf + orthxy * dyf + orthxz * dzf
        dy0 =                orthyy * dyf + orthyz * dzf
        dz0 =                               orthzz * dzf

        AtomGrids = []
        AtomGrid0 = [-1, -1, -1]
        ix0       = -1
        Grad0     = 0.0

        dz  = dz0
        r2z = dz * dz
        cz  = dz * StepZZ2 + StepZZS
        for iz in range(Kz1,Kz2) :
            dy   = dy0
            dt0  = dx0
            r2zy = dy * dy + r2z
            cy   = dy * StepYY2 + StepYYS

            for iy in range(Ky1,Ky2) :
                dx    = dt0
                r2zyx = dx * dx + r2zy
                cx    = dx * StepXX2 + StepXXS

                for ix in range(Kx1,Kx2) :

                    if r2zyx == 0.0 :
                       AtomGrid0 = [ix,iy,iz]
                       ix0       = ix
                       Grad0     = GradMap[iz][iy][ix]
                    elif r2zyx <= RadAtom2 :
                       rzyx  = math.sqrt(r2zyx)
                       rzyx2 = rzyx * 2
                       Gradxyz = GradMap[iz][iy][ix]
                       AtomGrids.append([r2zyx,rzyx,rzyx2,dx,dy,dz,Gradxyz])

                    r2zyx = r2zyx + cx
                    cx    = cx    + StepXXS2
                    dx    = dx    + StepXX

                r2zy = r2zy + cy
                cy   = cy   + StepYYS2
                dt0  = dt0  + StepXY
                dy   = dy   + StepYY

            r2z = r2z + cz
            cz  = cz  + StepZZS2
            dx0 = dx0 + StepXZ
            dy0 = dy0 + StepYZ
            dz  = dz  + StepZZ

#       calculate contributions

        Ngrids     = len(AtomGrids)
        GridValues = [0.0 for i in range(Ngrids)]
        GridValue0 = 0.0

        GradN  = [0.0 for ig in range(Ngrids)]
        GradC  = [0.0 for ig in range(Ngrids)]
        GradR  = [0.0 for ig in range(Ngrids)]
        GradR0 = [0.0 for ig in range(Ngrids)]
        GradC0, GradN0 = 0.0, 0.0

#        Nterms1, Nterms2 = TermsAtom[iatom]

        for iterm in range(0, bcr_scatterer.mu.size()) :
            mu    = bcr_scatterer.mu[iterm]
            nu    = bcr_scatterer.nu[iterm]
            musq  = bcr_scatterer.musq[iterm]
            kappi = bcr_scatterer.kappi[iterm]

            nuatom  = nu + bat
            nuatom2 = nuatom + nuatom
            munuat  = mu / nuatom
            fact1   = kappi / nuatom2**1.5

#           contribution from the node coinciding with the atomic center

            if AtomGrid0[0] >= 0 :
               argg  = musq / nuatom2
               fact2 = fact1 * math.exp(-argg)
               fact3 = fact2 / nuatom
               GradC0 = GradC0 + fact2
               GradN0 = GradN0 + fact3 * (argg - 1.5)

#           contribution from nodes different from the atomic center

            if mu == 0.0 :
               for ig in range(Ngrids) :
                   r2zyx = AtomGrids[ig][0]
                   argg  = r2zyx / nuatom2
                   fact2 = fact1 * math.exp(-argg)
                   fact3 = fact2 / nuatom
                   GradC[ig]  = GradC[ig]  + fact2
                   GradR0[ig] = GradR0[ig] + fact3
                   GradN[ig]  = GradN[ig]  + fact3 * (argg - 1.5)

            else :
               for ig in range(Ngrids) :
                   r2zyx,rzyx,rzyx2 = AtomGrids[ig][0],AtomGrids[ig][1],AtomGrids[ig][2]
                   arg0  = (rzyx-mu) / nuatom
                   argg  = (rzyx-mu) * arg0 * 0.5
                   fact2 = fact1 * math.exp(-argg)
                   tterm = rzyx2 * munuat
                   factt = math.exp(-tterm)
                   factc = (1.0 - factt) / tterm
                   GradC[ig] = GradC[ig] + fact2 *  factc
                   GradR[ig] = GradR[ig] + fact2 * (factc*(arg0 * rzyx + 1.0) - factt)
                   GradN[ig] = GradN[ig] + fact2 * (factc*(argg - 0.5)        - factt) / nuatom

#       collect scaled sums from different grid nodes

        GradXt, GradYt, GradZt, GradCt, GradNt = 0.0, 0.0, 0.0, 0.0, 0.0
        for ig in range(Ngrids) :
            r2zyx,rzyx,rzyx2,dx,dy,dz,Gradxyz = AtomGrids[ig]
            GradCt = Gradxyz *  GradC[ig] + GradCt
            GradNt = Gradxyz *  GradN[ig] + GradNt
            GradRg = Gradxyz * (GradR[ig]/r2zyx + GradR0[ig])
            GradXt = GradRg * dx + GradXt
            GradYt = GradRg * dy + GradYt
            GradZt = GradRg * dz + GradZt

        ix0,iy0,iz0 = AtomGrid0
        if ix0 >= 0 :
           GradCt = Grad0 * GradC0 + GradCt
           GradNt = Grad0 * GradN0 + GradNt

        GradAtom[iatom] = [GradXt * cat, GradYt * cat, GradZt * cat, GradCt, GradNt * cat]

    return GradAtom

#=====================================
def AtomBox(xfrac,yfrac,zfrac,RadAtomX,RadAtomY,RadAtomZ,StepX,StepY,StepZ,Sx,Sy,Sz,Fx,Fy,Fz) :

    x1 , y1 , z1  = (xfrac-RadAtomX)/StepX, (yfrac-RadAtomY)/StepY, (zfrac-RadAtomZ)/StepZ
    Kx1, Ky1, Kz1 = int(x1), int(y1), int(z1)
    if x1 >= 0.0 : Kx1 = Kx1 + 1
    if y1 >= 0.0 : Ky1 = Ky1 + 1
    if z1 >= 0.0 : Kz1 = Kz1 + 1

    x2 , y2 , z2  = (xfrac+RadAtomX)/StepX, (yfrac+RadAtomY)/StepY, (zfrac+RadAtomZ)/StepZ
    Kx2, Ky2, Kz2 = int(x2)+1, int(y2)+1, int(z2)+1

    if Kx1 < Sx : Kx1 = Sx
    if Ky1 < Sy : Ky1 = Sy
    if Kz1 < Sz : Kz1 = Sz

    if Kx2 > Fx : Kx2 = Fx
    if Ky2 > Fy : Ky2 = Fy
    if Kz2 > Fz : Kz2 = Fz

    dxf, dyf, dzf = Kx1*StepX - xfrac, Ky1*StepY - yfrac, Kz1*StepZ - zfrac

    Kx1, Ky1, Kz1 = Kx1 - Sx, Ky1 - Sy, Kz1 - Sz
    Kx2, Ky2, Kz2 = Kx2 - Sx, Ky2 - Sy, Kz2 - Sz

    return Kx1,Kx2,Ky1,Ky2,Kz1,Kz2,dxf,dyf,dzf
