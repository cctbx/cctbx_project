from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.shoebx

from pylab import *

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("modeler_file", type=str, help="path to a diffBragg modeler file (output from hopper, see the imgs folder in the outdir)")
parser.add_argument("--stacked", action="store_true")
parser.add_argument("--model_clim", nargs=2, default=[0,None], type=float)
parser.add_argument("--data_clim", nargs=2, default=[0,None], type=float)
args = parser.parse_args()


FIG,(ax0,ax1,ax2) = subplots(nrows=1,ncols=3)
FIG.set_size_inches((5,2))


def press(event):
    if event.key == 'right':
        FIG.loop_counter += 1
    elif event.key=="left":
        FIG.loop_counter = FIG.loop_counter -1
    FIG.loop_counter = max(FIG.loop_counter,0)
    FIG.loop_counter = min(FIG.loop_counter, FIG.nspots-1)

    if event.key=="escape":
        FIG.loop_counter = FIG.nspots

def centroid_poly(x,y):
    y1 = y[:-1]
    y2 = y[1:]

    x1 = x[:-1]
    x2 = x[1:]
    xy = (x1*y2 - x2*y1)
    Cx = np.sum((x1+x2)*xy)/6.
    Cy = np.sum((y1+y2)*xy)/6.
    A = 0.5*np.sum(xy)
    return Cx/A, Cy/A


M = np.load(args.modeler_file, allow_pickle=True)[()]
num_spots = len(M.pids)

FIG.loop_counter = 0
FIG.nspots = num_spots
FIG.canvas.mpl_connect('key_press_event', press)

assert len(set(M.roi_id)) == max(M.roi_id)+1
sigma_rdout = M.params.refiner.sigma_r / M.params.refiner.adu_per_photon
#M.hi = (0,0,0)]
#hkls = M.Hi
cmap = 'gray_r'

from cctbx import uctbx
# TODO put unit cell manager in the modeler file
#PSII
uc = uctbx.unit_cell((116.92, 221.635, 307.834, 90, 90, 90))
#lysozyme?
uc = uctbx.unit_cell((79.1,79.1,38.4,90,90,90))
diffs = []
d_max = 0
d_min = 9999

#F2 = plt.figure()
#AX2 = F2.gca()
#AX2.add_patch(Rectangle(xy=(0,0), width=4000, height=4000, fc='none', ec='b'))
#AX2.set_xlim(-10,4010)
#AX2.set_ylim(-10,4010)


if args.stacked:

    if isinstance(M.all_sigma_rdout, np.ndarray):
        data_subimg, model_subimg, trusted_subimg, bragg_subimg, sigma_rdout_subimg = M.get_data_model_pairs()
    else:
        data_subimg, model_subimg, trusted_subimg, bragg_subimg = M.get_data_model_pairs()
        sigma_rdout_subimg = None

    sub_sh = tuple(np.max([im.shape for im in model_subimg], axis=0))
    size_edg = int(np.sqrt(len(data_subimg))) + 1
    full_im = np.zeros((size_edg * sub_sh[0], size_edg * sub_sh[1]))
    full_dat_im = np.zeros((size_edg * sub_sh[0], size_edg * sub_sh[1]))
    for j in range(size_edg):
        for i in range(size_edg):
            mod_idx = j * size_edg + i
            if mod_idx >= len(model_subimg):
                continue
            im = model_subimg[mod_idx].copy()
            dat_im = data_subimg[mod_idx].copy()
            ydim, xdim = im.shape
            if im.shape != sub_sh:
                im = np.pad(im, ((0, sub_sh[0] - ydim), (0, sub_sh[1] - xdim)), mode='constant', constant_values=np.nan)
                dat_im = np.pad(dat_im, ((0, sub_sh[0] - ydim), (0, sub_sh[1] - xdim)), mode='constant',
                                constant_values=np.nan)
            Ysl = slice(j * sub_sh[0], (j + 1) * sub_sh[0], 1)
            Xsl = slice(i * sub_sh[1], (i + 1) * sub_sh[1], 1)
            full_im[Ysl, Xsl] = im
            full_dat_im[Ysl, Xsl] = dat_im
    subplot(121)
    gca().set_title("DATA", fontsize=18)
    imshow(full_dat_im, vmin=args.data_clim[0], vmax=args.data_clim[1], cmap='gray_r')
    gca().set_xticks(np.arange(sub_sh[1], size_edg*sub_sh[1], sub_sh[1]))
    gca().set_yticks(np.arange(sub_sh[0], size_edg*sub_sh[0], sub_sh[0]))
    gca().grid(1, color='k', ls='--')
    gca().set_xticklabels([])
    gca().set_yticklabels([])
    gca().tick_params(length=0)
    #gca().set_xlim(0, size_edg*sub_sh[1])
    #gca().set_ylim(0, size_edg*sub_sh[0])
    #gca().set_yticks([])
    subplot(122)
    gca().set_title("MODEL", fontsize=18)
    imshow(full_im, vmin=args.model_clim[0], vmax=args.model_clim[1], cmap='gnuplot')
    gca().set_xticks(np.arange(sub_sh[1], size_edg*sub_sh[1], sub_sh[1]))
    gca().set_yticks(np.arange(sub_sh[0], size_edg*sub_sh[0], sub_sh[0]))
    gca().grid(1, color='w', ls='--')
    gca().set_xticklabels([])
    gca().set_yticklabels([])
    gca().tick_params(length=0)
    #gca().set_xlim(0, size_edg*sub_sh[1])
    #gca().set_ylim(0, size_edg*sub_sh[0])
    #gca().set_xticks([])
    #gca().set_yticks([])
    subplots_adjust(wspace=0, hspace=0, bottom=0.02, top=0.95, left=0.01, right=0.99)
    show()
    exit()


a = b = None
while FIG.loop_counter < num_spots:
    i_h = FIG.loop_counter
    h = k = l = np.nan
    d = -1 #uc.d((h,k,l))
    #if d > d_max: d_max = d
    #if d < d_min: d_min = d

    # you could put a filter here to choose only those refls in a certain resolution
    sel = M.roi_id == i_h
    x1, x2, y1, y2 = M.rois[i_h]
    sh = y2 - y1, x2 - x1
    data = M.all_data[sel].reshape(sh)
    trusted = M.all_trusted[sel].reshape(sh)
    any_trusted = np.any(trusted)

    bg = M.all_background[sel].reshape(sh)
    bragg = M.best_model[sel].reshape(sh)
    model = bragg + bg
    diff = data.sum() - model.sum()

    Z = (model - data) / np.sqrt(model + sigma_rdout ** 2)
    diffs.extend(Z.ravel())

    vals = Z[trusted]
    data_trust = data[trusted]
    if not any_trusted:
        data_thresh = 0
    else:
        data_thresh = np.percentile(data_trust,97)
    y,x = np.indices(data.shape)

    ycent = (bragg.ravel()*y.ravel()).sum() /bragg.sum()
    xcent = (bragg.ravel()*x.ravel()).sum() / bragg.sum()

    if not any_trusted:
        sigZ_val = np.nan
    else:
        sigZ_val = Z[trusted].std()

    for ax in ax0,ax1,ax2:
        ax.clear()
    im = ax0.imshow(data, cmap=cmap)
    vmin, vmax = im.get_clim()
    ax1.imshow(model,cmap=cmap, vmin=vmin, vmax=vmax)

    ax1.contour(bragg, levels=[1e-4,1e-3,1e-2,1e-1,1,10], cmap='jet')
    C = ax0.contour(bragg, levels=[1e-4,1e-3,1e-2,1e-1,1,10], cmap='jet')
    ax1.plot( xcent, ycent, 'rx')
    #C = ax0.contour(bragg, levels=[0.1*data_thresh, 0.5*data_thresh,data_thresh], cmap='jet')
    #from IPython import embed;embed()

    #try:
    #    segs = C.allsegs[-1][0]
    #    xcent, ycent = centroid_poly(segs[:,0], segs[:,1])
    #    ax0.plot( [xcent], [ycent], 'o', mec='b', ms=10, alpha=0.5)
    #    ax0.plot( segs[:,0], segs[:,1], '.', ms=3, color='k')
    #except IndexError:
    #    pass

    v = np.max(np.abs(Z))
    v = 10
    ax2.imshow(Z, cmap='RdYlBu', vmin=-v, vmax=v)

    for AX in [ax0, ax1, ax2]:
        AX.grid(1, color='#777777', ls='--', lw=0.4)
        AX.set_xticklabels([])
        AX.set_yticklabels([])
    ax0.set_title("data %d,%d,%d,%d"% (x1,x2,y1,y2), pad=0, fontsize=8)
    ax1.set_title("model", pad=0, fontsize=8)
    ax2.set_title("Z (sigZ=%f)" % sigZ_val, pad=0, fontsize=8)

    FIG.suptitle("spot %d" % i_h)
    if a is not None:
        a.remove()
    if b is not None:
        b.remove()
    a = FIG.add_axes([.91, .25, 0.02, .5])
    FIG.colorbar(ax2.images[0], cax=a)
    a.tick_params(length=0, labelsize=8, pad=1)
    a.set_title("Z", pad=0, fontsize=8)

    b = FIG.add_axes([.2, .1, .33, .02])
    FIG.colorbar(ax0.images[0], cax=b, orientation='horizontal')
    b.yaxis.set_ticks_position('left')
    b.tick_params(length=0, labelsize=8, pad=5)
    b.set_title("counts", pad=0, fontsize=8)
    draw()
    waitforbuttonpress()

    i_h += 1

plt.close()
