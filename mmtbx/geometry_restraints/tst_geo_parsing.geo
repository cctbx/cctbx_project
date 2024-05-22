
# Geometry restraints

Bond restraints: 2
bond 38
     39
  ideal  model  delta    sigma   weight residual
  1.522  1.553 -0.030 1.18e-02 7.18e+03 6.53e+00
bond 37
     38
  ideal  model  delta    sigma   weight residual
  1.460  1.485 -0.025 1.17e-02 7.31e+03 4.40e+00

Bond angle restraints: 2
Sorted by residual:
angle 23
      24
      25
    ideal   model   delta    sigma   weight residual
   108.90  113.48   -4.58 1.63e+00 3.76e-01 7.90e+00
angle 37
      38
      39
    ideal   model   delta    sigma   weight residual
   108.02  111.93   -3.91 1.78e+00 3.16e-01 4.84e+00


Dihedral angle restraints: 2
  sinusoidal: 1
    harmonic: 1
Sorted by residual:
dihedral 24
         25
         37
         38
    ideal   model   delta  harmonic     sigma   weight residual
   180.00  166.21   13.79     0      5.00e+00 4.00e-02 7.60e+00
dihedral 58
         59
         60
         61
    ideal   model   delta sinusoidal    sigma   weight residual
     0.00  -72.39   72.39     2      3.00e+01 1.11e-03 4.85e+00


C-Beta improper torsion angle restraints: 2
Sorted by residual:
dihedral 37
         39
         38
         41
    ideal   model   delta  harmonic     sigma   weight residual
   122.80  126.95   -4.15     0      2.50e+00 1.60e-01 2.75e+00
dihedral 56
         54
         55
         58
    ideal   model   delta  harmonic     sigma   weight residual
  -122.60 -126.01    3.41     0      2.50e+00 1.60e-01 1.86e+00


Chirality restraints: 2
Sorted by residual:
chirality 55
          54
          56
          58
  both_signs  ideal   model   delta    sigma   weight residual
    False      2.51    2.39    0.12 2.00e-01 2.50e+01 3.48e-01
chirality 72
          71
          73
          75
  both_signs  ideal   model   delta    sigma   weight residual
    False      2.51    2.62   -0.11 2.00e-01 2.50e+01 2.86e-01

Planarity restraints: 2
Sorted by residual:
             delta    sigma   weight rms_deltas residual
plane 89     0.007 2.00e-02 2.50e+03   1.43e-02 6.15e+00
      90     0.029 2.00e-02 2.50e+03
      91    -0.003 2.00e-02 2.50e+03
      92     0.001 2.00e-02 2.50e+03
      93     0.003 2.00e-02 2.50e+03
      94    -0.001 2.00e-02 2.50e+03
      95    -0.013 2.00e-02 2.50e+03
      96    -0.001 2.00e-02 2.50e+03
      102   -0.029 2.00e-02 2.50e+03
      103   -0.016 2.00e-02 2.50e+03
      104    0.017 2.00e-02 2.50e+03
      105    0.005 2.00e-02 2.50e+03
            delta    sigma   weight rms_deltas residual
plane 13   -0.003 2.00e-02 2.50e+03   1.55e-02 3.62e+00
      14   -0.022 2.00e-02 2.50e+03
      15    0.019 2.00e-02 2.50e+03
      16   -0.002 2.00e-02 2.50e+03
      19   -0.012 2.00e-02 2.50e+03
      20    0.020 2.00e-02 2.50e+03


Nonbonded interactions: 3
Sorted by model distance:
nonbonded 106
          110
   model   vdw sym.op.
   1.719 1.850 -x+1,y-1/2,-z+1
nonbonded 29
          34
   model   vdw sym.op.
   1.859 1.850 x,y+1,z

nonbonded 33
          61
   model   vdw
   1.876 1.850

