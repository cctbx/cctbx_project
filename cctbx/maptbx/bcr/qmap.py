import math
import time
from scitbx.array_family import flex
import iotbx.map_manager
from libtbx.test_utils import approx_equal

######################################################
#
#    Calculating Q-Map from an atomic model with local resolution
#    and calculating respective gradients            by A.G. Urzhumtsev
#
######################################################

def InputData(nflog) :

#   come from the FileData

    FileAtoms = 'Atom2_Orth.pdb'
    FileMap   = 'Map_Orth_d25_B40.mrc'

#    FileMap   = 'none'
#    FileAtoms = 'Atom2_Orth_init.pdb'

    FileOut   = 'NewMap.mrc'

    RadFact = 2.0
    RadAdd  = 0.5
    Ngridx, Ngridy, Ngridz = 32,   40,   48

#   default values of the parameters

    acell,  bcell,  ccell  = 0.0, 0.0, 0.0
    alpha,  beta ,  gamma  = 0.0, 0.0, 0.0
    Mcol,   Mrow,   Msec   = 1,      2,      3
    Startx, Starty, Startz = 0,      0,      0
    Mmapx,  Mmapy,  Mmapz  = Ngridx, Ngridy, Ngridz
    Dmin, Dmax, Dmean, Rmsd   = 0.0, 0.0, 0.0, 0.0
    NumSpGr, BytesSG, FlagSkw = 0, 0, 0
    stamp                     = b'DA\x00\x00'

    cont_all = [0.0 for icon in range(66)]

    Ncrs = [Mmapx , Mmapy , Mmapz]
    Scrs = [Startx, Starty, Startz]
    Nxyz = [Ngridx, Ngridy, Ngridz]
    Mcrs = [Mcol,   Mrow,   Msec]
    Uabc = [acell,  bcell,  ccell]
    Uang = [alpha,  beta,   gamma]
    symm = []

    title = b' '*80

    cont_all[0] , cont_all[1] , cont_all[2]  = Ncrs
    cont_all[3]                              = 2
    cont_all[4] , cont_all[5] , cont_all[6]  = Scrs
    cont_all[7] , cont_all[8] , cont_all[9]  = Nxyz
    cont_all[10], cont_all[11], cont_all[12] = Uabc
    cont_all[13], cont_all[14], cont_all[15] = Uang
    cont_all[16], cont_all[17], cont_all[18] = Mcrs
    cont_all[19], cont_all[20], cont_all[21] = Dmin, Dmax, Dmean
    cont_all[22], cont_all[23], cont_all[24] = NumSpGr, BytesSG, FlagSkw
    cont_all[52]                             = b'MAP '
    cont_all[53]                             = stamp
    cont_all[54]                             = Rmsd
    cont_all[55]                             = 0

    for icon in range(56,66): cont_all[icon] = title

#   calculated internal parameters

    RadMu   = RadFact + RadAdd

    Uang[0] = Uang[0] * math.pi/180.0
    Uang[1] = Uang[1] * math.pi/180.0
    Uang[2] = Uang[2] * math.pi/180.0

#   save parameters

    return RadFact, RadMu, Ncrs,Scrs,Nxyz,Uabc,Uang,cont_all, FileMap, FileAtoms, FileOut

#=====================================
def GetAtomicModel(FileAtoms,Uabc,Uang) :

    print("LOOK FileAtoms", FileAtoms)

    nfatom = open(FileAtoms, 'r')

    ModelValues = []
    ModelTypes  = []

    ScaleB = 1.0 / (8.0 * math.pi**2)

#   read title

    Line = nfatom.readline()
    LineCut  = Line.rstrip('\n')
    LineCont = LineCut.split()
    acell, bcell, ccell = float(LineCont[1]), float(LineCont[2]), float(LineCont[3])
    alpha, beta , gamma = float(LineCont[4]), float(LineCont[5]), float(LineCont[6])

    if Uabc[0] == 0.0 :
       Uabc    = [acell, bcell, ccell]
       Uang[0] = alpha * math.pi/180.0
       Uang[1] = beta  * math.pi/180.0
       Uang[2] = gamma * math.pi/180.0

    while (True) :
        Line = nfatom.readline()
        if not Line : break

        LineCut  = Line.rstrip('\n')
        LineCont = LineCut.split()
        if len(LineCont) == 0 : break

        xat = float(LineCont[6])
        yat = float(LineCont[7])
        zat = float(LineCont[8])
        cat = float(LineCont[9])
        bat = float(LineCont[10])
        rat = float(LineCont[11])
        tat =       LineCont[12]
#
        nat = bat * ScaleB

#       save current line content

        ModelValues.append([xat, yat, zat, cat, nat])
        ModelTypes.append([tat, rat])

    nfatom.close( )

    print("LOOK ModelValues:",ModelValues)
    print("LOOK ModelTypes:",ModelTypes)
    print("LOOK Uabc,Uang:",Uabc,Uang)

    return ModelValues,ModelTypes,Uabc,Uang

#============================
def GetMNKCoefficients(ModelTypes,RadFact,RadMu) :

#   Array ModelTypes consists of a line per atom.
#   This line contains :

#   Types contains labels of atomic types, as they found in the list of atoms
#   ResMin, ResMax contains min and max resolution for the given atomic type
#   TermsDecomp contains consecutively all decomposition coefficients in the order
#        of the Types and in the limits (ResMin, ResMax) for each type;
#        for each atomic type, groups of coefficients are given increasing resolution
#   PosTerms is the starting position in TermsDecomp of the respective group of terms
#   ResTerms is the resolution of this group of terms
#   TypeStart is the position of the respective atomic type in (PosTerms, ResTerms)

#   extract atomic types and resolution limits for each type

    Types,ResMin,ResMax,TermsAtom = FindTypes(ModelTypes)

#   read coefficients for the required atomic types and the resolution ranges

    TermsDecomp,TypeStart,ResTerms,PosTerms = GetDecomposition(Types,ResMin,ResMax)

#   establish references from each atom to respective lines of the global table

    TermsAtom = GetTables(TermsDecomp,TypeStart,ResTerms,PosTerms,RadFact,RadMu,TermsAtom)

    return TermsAtom, TermsDecomp

#============================
def FindTypes(ModelTypes) :

    Natoms = len(ModelTypes)
    Ntypes = 0
    Types, ResMin, ResMax = [], [], []

    TermsAtom = [ [0.0, 0, 0] for ia in range(Natoms)]

#   find atomic types and their resolution bounds

    for iat in range(Natoms) :
        TypeAtom, ResAtom = ModelTypes[iat]
        TermsAtom[iat][0] = ResAtom
        FoundType = False
        for itype in range(Ntypes) :
            if TypeAtom == Types[itype] :
               FoundType = True
               break
        if FoundType :
           if ResMin[itype] > ResAtom : ResMin[itype] = ResAtom
           if ResMax[itype] < ResAtom : ResMax[itype] = ResAtom
           TermsAtom[iat][1] = itype
        else :
           Types.append(TypeAtom)
           ResMin.append(ResAtom)
           ResMax.append(ResAtom)
           TermsAtom[iat][1] = Ntypes
           Ntypes            = Ntypes + 1

    return Types,ResMin,ResMax, TermsAtom

#============================
def GetDecomposition(Types,ResMin,ResMax) :

    Ntypes = len(Types)

    ScaleB = 1.0 / (8.0 * math.pi**2)

#   read file and save terms for the given types and resolution within bounds

    iterm       = 0
    ilist       = 0
    TermsDecomp = []
    ResTerms    = []
    PosTerms    = []
    TypeStart   = []

    for itype in range(Ntypes) :
        TypeMin = ResMin[itype]
        TypeMax = ResMax[itype]
        FileTerms = Types[itype] + '_MNK.tab'

        nfterms  = open(FileTerms, 'r')
        TypeStart.append(ilist)

        while (True) :
           Line = nfterms.readline()
           if not Line : break

           LineClean = Line.rstrip('\n')
           LineCut   = LineClean.replace("="," ")
           LineCont  = LineCut.split()
           if len(LineCont) == 0 : break

           if LineCont[0] == 'Atom' :
              ResTable = float(LineCont[5])
              if ResTable > TypeMax : break
              TakeLine = False
              if ResTable >= TypeMin :
                 TakeLine = True
                 ResTerms.append(ResTable)
                 PosTerms.append(iterm)
                 ilist = ilist + 1
              Line = nfterms.readline()

           else :
              if TakeLine :
                 mu, sigma, kappa = float(LineCont[1]), float(LineCont[2]), float(LineCont[3])
                 nu    = sigma * ScaleB
                 musq  = mu * mu
                 kappi = kappa / math.pi**1.5
                 TermsDecomp.append([mu,nu,kappa,musq,kappi])
                 iterm = iterm + 1

        nfterms.close()

    return TermsDecomp, TypeStart, ResTerms, PosTerms

#============================
def GetTables(TermsDecomp, TypeStart, ResTerms, PosTerms,RadFact,RadMu,TermsAtom) :

    Natoms = len(TermsAtom)
    Ntypes = len(TypeStart)
    Nlists = len(PosTerms)
    Nterms = len(TermsDecomp)

    for iatom in range(Natoms) :
        ResAtom, itype, dummy = TermsAtom[iatom]
        RadAtom, MuMax        = ResAtom * RadFact , ResAtom * RadMu

#       find the part of the list (PosTerms,ResTerms) corresponding to the given type

        k1 = TypeStart[itype]
        if itype == Ntypes-1 :
           k2 = Nlists
        else :
           k2 = TypeStart[itype+1]

#       find the line (PosTerms) for the given type and resolution

        for kres in range(k2-1, k1-1, -1) :
            if ResAtom >= ResTerms[kres] :
               k0 = kres
               break

        if k0 < k2-1 :
           if (ResAtom - ResTerms[k0]) > (ResTerms[k0+1] - ResAtom ) : k0 = k0 + 1

#       select from TermsDecomp all terms within the given maximal distance

        N1 = PosTerms[k0]
        if k0 == Nlists-1 :
           NTMax = Nterms
        else :
           NTMax = PosTerms[k0+1]

        N2 = NTMax
        for it in range (N1,NTMax) :
            if TermsDecomp[it][0] > MuMax :
               N2 = it
               break

        TermsAtom[iatom] = [RadAtom, N1, N2]

    return TermsAtom

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
def CalcOmegaMap(ModelValues,TermsDecomp,TermsAtom,Ncrs, Scrs, Nxyz, unit_cell) :
    print()
    print("Ncrs:", Ncrs)
    print("Scrs:", Scrs)
    print("Nxyz:", Nxyz)
    print("Uabc:", Uabc)
    print("Uang:", Uang)
    print()


    Natoms = len(ModelValues)

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

    for iatom in range(Natoms) :

        RadAtom, Nterms1, Nterms2 = TermsAtom[iatom]
        RadAtom2  = RadAtom   * RadAtom
        RadAtomX  = RadAtom   * RprojX
        RadAtomY  = RadAtom   * RprojY
        RadAtomZ  = RadAtom   * RprojZ

        xat, yat, zat, cat, bat = ModelValues[iatom]

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

        for iterm in range(Nterms1, Nterms2) :
            mu, nu, kappa, musq, kappi = TermsDecomp[iterm]
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
            OmegaMap[iz][iy][ix] = GridValues[ig] * cat

        ix0,iy0,iz0 = AtomGrid0
        if ix0 >= 0 :
           OmegaMap[iz0][iy0][ix0] = GridValue0 * cat

    return OmegaMap

#=====================================
def CalcGradAtom(GradMap,ModelValues,TermsDecomp,TermsAtom,Ncrs, Scrs, Nxyz, unit_cell):

    Natoms = len(ModelValues)

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

    for iatom in range(Natoms) :

        RadAtom, Nterms1, Nterms2 = TermsAtom[iatom]
        RadAtom2  = RadAtom   * RadAtom
        RadAtomX  = RadAtom   * RprojX
        RadAtomY  = RadAtom   * RprojY
        RadAtomZ  = RadAtom   * RprojZ

        xat, yat, zat, cat, bat = ModelValues[iatom]

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

        for iterm in range(Nterms1, Nterms2) :
            mu, nu, kappa, musq, kappi = TermsDecomp[iterm]
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

#=====================================
def Deorth(xat,yat,zat,DeortMatrix) :

    xfrac = xat * DeortMatrix[0] + yat * DeortMatrix[3] + zat * DeortMatrix[6]
    yfrac = xat * DeortMatrix[1] + yat * DeortMatrix[4] + zat * DeortMatrix[7]
    zfrac = xat * DeortMatrix[2] + yat * DeortMatrix[5] + zat * DeortMatrix[8]

    return xfrac, yfrac, zfrac

#=====================================
def GetOrthDeorth(Uabc,Uang) :

    DeortMatrix = [0.0 for i in range(9)]
    OrthMatrix  = [0.0 for i in range(9)]

    acell, bcell, ccell = Uabc
    alpr,  betr,  gamr  = Uang

#   deorthogonalization matrix (Angstroms -> fract.coord)

    calr = math.cos(alpr)
    salr = math.sin(alpr)
    cber = math.cos(betr)
    sber = math.sin(betr)
    cgar = math.cos(gamr)
    sgar = math.sin(gamr)

    cha = (calr-cber*cgar) / (sber*sgar)
    if cha >  1.0 : cha=  1.0
    if cha < -1.0 : cha= -1.0
    sha = math.sqrt(1.0 - cha*cha)

    chb =(cber-calr*cgar) /( salr*sgar)
    if chb >  1.0 : chb =  1.0
    if chb < -1.0 : chb = -1.0
    shb = math.sqrt(1.0 - chb*chb)

    chg = (cgar-calr*cber) / (sgar*sber)
    if chg >  1.0 : chg =  1.0
    if chg< -1.0  : chg = -1.0
    shg = math.sqrt(1.0 - chg*chg)

#   calcualtion of the matrices

    acelli, bcelli, ccelli = 1.0/acell , 1.0/bcell , 1.0/ccell

    DeortMatrix[0] =  acelli
    DeortMatrix[3] = -acelli * cgar       /  sgar
    DeortMatrix[4] =  bcelli              /  sgar
    DeortMatrix[6] = -acelli * salr * chb / (sgar * sber * sha)
    DeortMatrix[7] = -bcelli * cha        / (sgar * sha)
    DeortMatrix[8] =  ccelli              / (sber * sha)

#   orthogonalization matrix (fract.coord -> angstroms)

    OrthMatrix[0] = acell
    OrthMatrix[3] = bcell * cgar
    OrthMatrix[4] = bcell * sgar
    OrthMatrix[6] = ccell * cber
    OrthMatrix[7] = ccell * sber * cha
    OrthMatrix[8] = ccell * sber * sha

    return DeortMatrix, OrthMatrix


####################################################
#
#  MAIN PROGRAM

#==================== input information =====================

time_t01 = time.time()

nflog = open('omega_refine.log', 'w')

RadFact,RadMu,Ncrs,Scrs,Nxyz,Uabc,Uang,cont_all,FileMap,FileAtoms,FileOut = InputData(nflog)

print('Read control / experimental map...')

mm = iotbx.map_manager.map_manager(file_name = FileMap)
Ncrs = mm.map_data().all()
Scrs = [0,0,0]
Nxyz = mm.map_data().all()
unit_cell = mm.crystal_symmetry().unit_cell()
Uabc = unit_cell.parameters()[:3]
Uang = unit_cell.parameters()[3:]
Uang = [it*math.pi/180. for it in Uang]


Mx,My,Mz = Nxyz
ControlMap = [[[ mm.map_data()[ix,iy,iz] for ix in range(Mx) ] for iy in range(My)] for iz in range(Mz)]


print('Read atomic model...')

ModelValues,ModelTypes,Uabc,Uang = GetAtomicModel(FileAtoms,Uabc,Uang)

#=============== principal part ==========

# INPUR DATA REQUIRED :

# ModelValues,ModelTypes     - both contain a line per atom
#                              [x, y, z, c, nu] and [scattering_type, atomic_resolution]
#                              x, y, z - Cartesian coordinates; nu = B / (8*pi*pi)
#                              scattering_type is C, N, O etc (as that in positions 78-79 PDB)
# RadFact,RadMu              - atomic radius is calculated as atomic_resolution * RadFact
#                            - maximal mu values is taken as atomic_resolution * (RadFact+RadMu)
# ControlMap                 - input (control) map; a map is ALWAYS kept as ;
#                              columns (fastest index) by x, sections (slowest index) by z
# Uabc = [acell,bcell,ccell] - unit cell lengths (in Angstorm )
# Uang = [alpha,beta,gamma]  - unit cell angles (in radians)
# Ncrs = [Ncol,Nrow,Nsec]    - number of columns, rows and sections;
# Scrs = [Scol,Srow,Ssec]    - starting values of the indices;
# Nxyz = [Nx,  Ny,  Nz]      - number of grid intervals by x, y, z

# INTERNAL KEY ARRAYS : (eventually, can be calculated once at the beginning)

# TermsDecomp - a table of decomposition parameters for all required types and resolution)
# TermsAtoms  - reference from each atom to the respective part of TermsDecomp

# OUTPUT DATA :

# OmegaMam               - map of a limited (variable) resolution calculated from an atomic model
# GradAtom               - array containing a line per atom of partial derivatives by x,y,z,c,nu ;
#                          to calculate to the derivative by B, dR/dB = (dR/dnu) / (8*pi*pi)

#-------------- calculation --------

TermsAtom, TermsDecomp = GetMNKCoefficients(ModelTypes,RadFact,RadMu)

OmegaMap = CalcOmegaMap(ModelValues,TermsDecomp,TermsAtom,Ncrs,Scrs,Nxyz, unit_cell)

FuncMap  = CalcFuncMap(OmegaMap, ControlMap, Ncrs)
GradMap  = CalcGradMap(OmegaMap, ControlMap, Ncrs)

GradAtom = CalcGradAtom(GradMap,ModelValues,TermsDecomp,TermsAtom,Ncrs,Scrs,Nxyz,unit_cell)

#=============== output results ======================================

Natoms = len(GradAtom)
for iatom in range(Natoms) :
        print('Model',iatom,ModelValues[iatom])

assert approx_equal(GradAtom[0], [0.6210971239199201, -1.2420777490312254, -2.484579744482787, 6.2528792787377645, -7.192527753707328])
assert approx_equal(GradAtom[1], (0,0,0,0,0), 1.e-4 )
assert approx_equal(GradAtom[2], (0,0,0,0,0), 1.e-4 )

####
assert approx_equal(FuncMap, 8.219926081295183)
#
nx,ny,nz = Nxyz
m0 = flex.double(flex.grid(Nxyz))
for iz in range(0, nz):
  for iy in range(0, ny):
    for ix in range(0, nx):
      v = OmegaMap[iz][iy][ix]
      m0[ix,iy,iz] = v

print(Nxyz)
mm = iotbx.map_manager.map_manager(file_name = FileOut)
m1 = mm.map_data()
print(m1.all())

diff = flex.abs(m1-m0)
print("max:", flex.max(diff))
assert flex.max(diff) < 0.005
for d in diff:
  assert d < 0.005


#assert approx_equal(m1, m0, 1.e-3)
#
mm2 = iotbx.map_manager.map_manager(
  map_data                   = m1,
  unit_cell_grid             = m1.all(),
  unit_cell_crystal_symmetry = mm.crystal_symmetry(),
  wrapping                   = True)


time_t08 = time.time()

print('')
print('time total                  ...... ',f'{(time_t08-time_t01)/60.:12.2f} min')
