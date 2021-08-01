import pandas
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--pkl", type=str, help="input pickle")
parser.add_argument("--erf", type=str, default=None, help="expref file for filtering")
parser.add_argument("--minscore", type=float, default=None, help="minimum prediction offset loss")
args = parser.parse_args()

from simtbx.diffBragg.utils import is_outlier
import sys
from pylab import *

thresh_a = 10
thresh_c = 10
thresh_N = 5
thresh_Na = 6
thresh_Nb = 7
thresh_Nc = 7
thresh_d = 7
thresh_s = 1350
thresh_rot = 45
thresh_rotx = 25
thresh_roty = 25
thresh_rotz = 25

tag = None
#if len(sys.argv)==3:
#    tag = sys.argv[2]
outname = args.pkl.replace(".pkl", "_noout.pkl")
df = pandas.read_pickle(args.pkl)
figure()
cmin = min(df.c.min(), df.c_init.min())
cmax = min(df.c.max(), df.c_init.max())
c_out = is_outlier(df.c, thresh_c)
bins = linspace(cmin, cmax, 125)
df.c.hist(bins=bins, log=True, label="c (after stage 1)")
df.c_init.hist(bins=bins, log=True, histtype='step', label="c (from DIALS)")
legend(prop={"size": 11})
xlabel("$\AA$",fontsize=14)

figure()
amin = min(df.a.min(), df.a_init.min())
amax = min(df.a.max(), df.a_init.max())
bins =  linspace(amin, amax,125)
df.a.hist(bins=bins, log=True, label="a (after stage 1)")
df.a_init.hist(bins=bins, log=True,histtype='step', label="a (from DIALS)")
legend(prop={"size":11})
xlabel("$\AA$",fontsize=14)
a_out = is_outlier(df.a, thresh_a)

if 'detz_shift_mm' in list(df):
    figure()
    df.detz_shift_mm.hist(bins=100, log=True)
    xlabel("detector Z shift (millimeters)", fontsize=14)
    det_out = is_outlier(df.detz_shift_mm, thresh_d)

na,nb,nc = np.vstack([v for v in df['ncells']]).T
#figure()
df['Na'] = na
df['Nb'] = nb
df['Nc'] = nc
df['Na x Nb x Nc'] = na*nb*nc
df[['Na','Nb','Nc', 'Na x Nb x Nc']].hist(bins=100, log=True)
N_out = is_outlier(df['Na x Nb x Nc'], thresh_N)
Na_out = is_outlier(df['Na'], thresh_Na)
Nb_out = is_outlier(df['Nb'], thresh_Nb)
Nc_out = is_outlier(df['Nc'], thresh_Nc)

figure()
s = df.spot_scales.values
bins = logspace(log10(s.min()), log10(s.max()), 100)
df.spot_scales.hist(bins=bins, log=True)
ax = gca()
ax.set_xscale("log")
s_out = is_outlier(df.spot_scales, thresh_s)
xlabel("spot scales", fontsize=14)

rotX = df.rotX.values*180 / pi
rotY = df.rotY.values*180 / pi
rotZ = df.rotZ.values*180 / pi
rotxyz = np.sqrt(rotX**2 + rotY**2 + rotZ**2)
df['rx (deg.)'] = rotX
df['ry (deg.)'] = rotY
df['rz (deg.)'] = rotZ
df['|rx,ry,rz| (deg.)'] = rotxyz
df[['rx (deg.)','ry (deg.)','rz (deg.)', '|rx,ry,rz| (deg.)']].hist(bins=100, log=True)
rot_out = is_outlier(df['|rx,ry,rz| (deg.)'],  thresh_rot)
rx_out = is_outlier(df['rx (deg.)'],  thresh_rotx)
ry_out = is_outlier(df['rx (deg.)'],  thresh_roty)
rz_out = is_outlier(df['rx (deg.)'],  thresh_rotz)
ax = gca()
legend(prop={"size":11})
for i in get_fignums():
    figure(i)
    suptitle(args.pkl)
    axs = gcf().get_axes()
    for ax in axs:
        ax.set_ylabel("# of shots", fontsize=12)

print("Outliers:")
if "detz_shift_mm" in list(df):
    print("det_z: %d" % sum(det_out))
    df_d = df.loc[~det_out]
    print("\t detz shift  range of values afer outlier removal: %f - %f" % (df_d.detz_shift_mm.min(), df_d.detz_shift_mm.max()))
    print("\t detz shift median=%f, mean=%f, std=%f (millimeters)" % (df_d.detz_shift_mm.median(), df_d.detz_shift_mm.mean(), df_d.detz_shift_mm.std()) )
    print("\t suitable beta for detz_shift: %e (meters)" % ((df_d.detz_shift_mm.std()*1e-3)**2))

print("unit cell a: %d" % sum(a_out))
df_a = df.loc[~a_out]
print("\t ucell a range of values afer outlier removal: %f - %f" % (df_a.a.min(), df_a.a.max()))
print("\t ucell a median=%f, mean=%f, std=%f (Angstrom)" % (df_a.a.median(), df_a.a.mean(), df_a.a.std()) )
print("\t suitable beta for ucell a: %e (Angstrom)" % (df_a.a.std()**2))

print("unit cell c: %d" % sum(c_out))
df_c = df.loc[~c_out]
print("\t ucell c range of values afer outlier removal: %f - %f" % (df_c.c.min(), df_c.c.max()))
print("\t ucell c median=%f, mean=%f, std=%f (Angstrom)" % (df_c.c.median(), df_c.c.mean(), df_c.c.std()) )
print("\t suitable beta for ucell c: %e (Angstrom)" % (df_c.c.std()**2))

print("Na: %d" % sum(Na_out))
print("Nb: %d" % sum(Nb_out))
print("Nc: %d" % sum(Nc_out))
df_n = df.loc[~N_out]
df_na = df.loc[~Na_out]
df_nb = df.loc[~Nb_out]
df_nc = df.loc[~Nc_out]
namin = df_na['Na'].min()
nbmin = df_nb['Nb'].min()
ncmin = df_nc['Nc'].min()

namax = df_na['Na'].max()
nbmax = df_nb['Nb'].max()
ncmax = df_nc['Nc'].max()
print("\t Na range of values afer outlier removal: %f - %f" % (namin, namax))
print("\t Nb range of values after outlier removal: %f - %f" % (nbmin, nbmax))
print("\t Nc range of values after outlier removal: %f - %f" % (ncmin, ncmax))
print("\t Na median=%f, mean=%f, std=%f" % (df_na.Na.median(), df_na.Na.mean(), df_na.Na.std() ) )
print("\t Nb median=%f, mean=%f, std=%f" % (df_nb.Nb.median(), df_nb.Nb.mean(), df_nb.Nb.std()) )
print("\t Nc median=%f, mean=%f, std=%f" % (df_nc.Nc.median(), df_nc.Nc.mean(), df_nc.Nc.std()) )
print("\t suitable beta for Na: %e" % (df_na.Na.std()**2))
print("\t suitable beta for Nb: %e" % (df_nb.Nb.std()**2))
print("\t suitable beta for Nc: %e" % (df_nc.Nc.std()**2))
print("Spot scale G: %d" % sum(s_out))
df_s = df.loc[~s_out]
print("\t Spot scale range of values afer outlier removal: %f - %f" % (df_s.spot_scales.min(), df_s.spot_scales.max()))
print("\t Spot scale median=%f, mean=%f, std=%f" % (df_s.spot_scales.median(), df_s.spot_scales.mean(), df_s.spot_scales.std()) )
print("\t suitable beta for spot scale: %e" % (df_s.spot_scales.std()**2))

print("RotXYZ: %d" % sum(rot_out))
print("RotX: %d" % sum(rx_out))
print("RotY: %d" % sum(ry_out))
print("RotZ: %d" % sum(rz_out))
df_rot = df.loc[~rot_out]
df_rotx = df.loc[~rx_out]
df_roty = df.loc[~ry_out]
df_rotz = df.loc[~rz_out]
rxmin = df_rot['rx (deg.)'].min()
rymin = df_rot['ry (deg.)'].min()
rzmin = df_rot['rz (deg.)'].min()

rxmax = df_rot['rx (deg.)'].max()
rymax = df_rot['ry (deg.)'].max()
rzmax = df_rot['rz (deg.)'].max()
print("\t RotX range of values after outlier removal: %f - %f" % (rxmin, rxmax))
print("\t RotY range of values after outlier removal: %f - %f" % (rymin, rymax))
print("\t RotZ range of values after outlier removal: %f) - %f" % (rzmin, rzmax))
rrx = df_rot['rx (deg.)']
rry = df_rot['ry (deg.)']
rrz = df_rot['rz (deg.)']
print("\t RotX median=%f, mean=%f, std=%f" % (rrx.median(), rrx.mean(), rrx.std() ) )
print("\t RotY median=%f, mean=%f, std=%f" % (rry.median(), rry.mean(), rry.std() ) )
print("\t RotZ median=%f, mean=%f, std=%f" % (rrz.median(), rrz.mean(), rrz.std() ) )
print("\t Suitable beta for rotX: %e (rad.)" % ((rrx.std()   )**2))
print("\t Suitable beta for rotY: %e (rad.)" % ((rry.std()   )**2))
print("\t Suitable beta for rotZ: %e (rad.)" % ((rrz.std()   )**2))
total_unique_outliers = np.logical_or(a_out, c_out)
total_unique_outliers = np.logical_or(total_unique_outliers, Na_out)
total_unique_outliers = np.logical_or(total_unique_outliers, Nb_out)
total_unique_outliers = np.logical_or(total_unique_outliers, Nc_out)
total_unique_outliers = np.logical_or(total_unique_outliers, s_out)
total_unique_outliers = np.logical_or(total_unique_outliers, rot_out)
if "detz_shift_mm" in list(df):
    total_unique_outliers = np.logical_or(total_unique_outliers, det_out)

df_no_out = df.loc[~total_unique_outliers]
print("Total unique outliers: %d / %d" % (sum(total_unique_outliers), len(df)))

if args.minscore is not None:
    df_no_out['score'] = [np.median(d_dials) - np.median(d)
                   for d, d_dials in zip(df_no_out.pred_offsets, df_no_out.pred_offsets_dials)]
    above_minscore = df_no_out.score > args.minscore
    df_no_out = df_no_out.loc[above_minscore]
    print("Shots below minscore: %d" % (sum(~above_minscore)))

df_no_out.to_pickle(outname)
print("Wrote outlier filterted pkl %s" % outname)
import os
if tag is not None:
    exp_ref_name = outname.replace(".pkl", "_histfilt.txt")
    o = open(exp_ref_name, "w")
    for e,eopt,s in zip(df_no_out.exp_name, df_no_out.opt_exp_name, df_no_out.spectrum_filename):
        r = eopt.replace(".expt", "_%s.refl" % tag)
        assert os.path.exists(r)
        o.write("%s %s %s\n" % (eopt,r,s))
    o.close()
    print("Wrote exp ref file %s" % exp_ref_name)


if args.erf is not None:
    lines = open(args.erf, "r").readlines()
    good_e = set(df_no_out.exp_name.values)
    new_ers_name = os.path.splitext(args.erf)[0]+ "_noout.txt"
    new_ers = open(new_ers_name, "w")
    for l in lines:
        lsplit = l.strip().split()
        if len(lsplit)==2:
            e,r = lsplit
            s = ""
        else:
            e,r,s = l.strip().split()
        if e not in good_e:
            continue
        new_ers.write(l)
    new_ers.close()
    print("Wrote new file %s" % new_ers_name)

show()

