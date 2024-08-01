from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.shoebx

from pylab import *

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("modeler_file", type=str, help="path to a diffBragg modeler file (output from hopper, see the imgs folder in the outdir)")
parser.add_argument("--scroll", action="store_true", help="if provided, scroll through shoeboxes one-by-one using arrow keys")
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


if not args.scroll:

    if isinstance(M.all_sigma_rdout, np.ndarray):
        data_subimg, model_subimg, trusted_subimg, bragg_subimg, sigma_rdout_subimg, stats = M.get_data_model_pairs(reorder=True, return_stats=True)
    else:
        data_subimg, model_subimg, trusted_subimg, bragg_subimg, stats = M.get_data_model_pairs(reorder=True, return_stats=True)
        sigma_rdout_subimg = None

    sub_sh = tuple(np.max([im.shape for im in model_subimg], axis=0))
    size_edg = int(np.sqrt(len(data_subimg))) + 1
    full_im = np.zeros((size_edg * sub_sh[0], size_edg * sub_sh[1]))
    full_dat_im = np.zeros((size_edg * sub_sh[0], size_edg * sub_sh[1]))
    full_trust_im = np.zeros((size_edg * sub_sh[0], size_edg * sub_sh[1])).astype(bool)
    labs = {"res":[], "hkl":[], "pfs":[], "sigZ":[], "xy":[], "dat_max":[], "mod_max":[]}
    contour_kwargs = []
    for j in range(size_edg):
        for i in range(size_edg):
            mod_idx = j * size_edg + i


            if mod_idx >= len(model_subimg):
                continue
            im = model_subimg[mod_idx].copy()
            dat_im = data_subimg[mod_idx].copy()
            trust_im = trusted_subimg[mod_idx].copy()
            bragg_im = bragg_subimg[mod_idx].copy()
            Z = (im - dat_im) / np.sqrt(im + sigma_rdout ** 2)
            sigmaZ = Z[trust_im].std()

            # make label strings for later use...
            spot_sigZ_lab = "%.1f" % sigmaZ
            spot_res_lab = "%.1f" % (stats["spot_d"][mod_idx])
            spot_hkl_lab = "%d,%d,%d"% (tuple(stats["spot_hkl"][mod_idx]))
            spot_pid_cent_lab = "%d,%d,%d"% (tuple(stats["spot_pid_and_cent"][mod_idx]))
            dat_max_lab = mod_max_lab = "nan"
            if trust_im.sum()>0:
                dat_max_lab = "%d" % dat_im[trust_im].max()
                mod_max_lab = "%d" % im[trust_im].max()

            # loop over sub images
            ydim, xdim = im.shape
            if im.shape != sub_sh:
                im = np.pad(im, ((0, sub_sh[0] - ydim), (0, sub_sh[1] - xdim)), mode='constant', constant_values=np.nan)
                dat_im = np.pad(dat_im, ((0, sub_sh[0] - ydim), (0, sub_sh[1] - xdim)), mode='constant',
                                constant_values=np.nan)
                trust_im = np.pad(trust_im, ((0, sub_sh[0] - ydim), (0, sub_sh[1] - xdim)), mode='constant',
                                constant_values=False)
            Ysl = slice(j * sub_sh[0], (j + 1) * sub_sh[0], 1)
            Xsl = slice(i * sub_sh[1], (i + 1) * sub_sh[1], 1)
            full_im[Ysl, Xsl] = im
            full_dat_im[Ysl, Xsl] = dat_im
            full_trust_im[Ysl,Xsl] = trust_im

            #dat_lvls = [1e-4,1e-3,1e-2,1e-1, 1,10, 100] #np.logspace( log10(bragg_im.min()-1e-6), log10(bragg_im.max()*0.5 +1e-6), 5)
            dat_lvls = np.logspace( log10(dat_im.min()-1e-6), log10(dat_im.max() +1e-6), 5)[:-1]
            mod_lvls = np.logspace( log10(im.min()-1e-6), log10(im.max() +1e-6), 5)[:-1]
            ycoord, xcoord = np.indices(sub_sh)
            ycoord += Ysl.start
            xcoord += Xsl.start
            contour_kwargs.append({'X':xcoord, 'Y':ycoord, 'dat_Z':dat_im, 'mod_Z': im, 'dat_levels': dat_lvls, 'mod_levels': mod_lvls })

            labs["res"].append( spot_res_lab)
            labs["hkl"].append(spot_hkl_lab)
            labs["pfs"].append(spot_pid_cent_lab)
            labs["sigZ"].append(spot_sigZ_lab)
            labs["dat_max"].append(dat_max_lab)
            labs["mod_max"].append(mod_max_lab)
            xlab = (Xsl.start+Xsl.stop)*.5
            ylab = Ysl.start + (Ysl.stop-Ysl.start)*.25
            labs["xy"].append( (xlab, ylab))

    subplot(121)
    gca().set_facecolor('tomato')
    gca().set_title("DATA", fontsize=18)
    masked_dat_vals = full_dat_im[~full_trust_im]
    full_dat_im[~full_trust_im] = np.nan
    imshow(full_dat_im, vmin=args.data_clim[0], vmax=args.data_clim[1], cmap='gray_r')
    xt = np.arange(sub_sh[1], size_edg*sub_sh[1], sub_sh[1])-0.5
    yt = np.arange(sub_sh[0], size_edg * sub_sh[0], sub_sh[0])-0.5
    gca().set_xticks(xt)
    gca().set_yticks(yt)
    dx = abs(xt[1]-xt[0])/2
    dy = abs(yt[1]-yt[0])/2
    mid_x = list((xt[1:]+xt[:-1])//2)
    mid_y = list((yt[1:]+yt[:-1])//2)
    mid_x = [xt[0]-dx] + mid_x + [xt[-1]+dx]
    mid_y = [yt[0]-dy] + mid_y + [yt[-1]+dy]
    gca().set_xticks(mid_x, minor=True)
    gca().set_yticks(mid_y, minor=True)
    gca().grid(1, color='k', ls='--', which='major')
    gca().grid(1, color='k',alpha=0.5,lw=0.5, ls='--', which='minor')
    gca().set_xticklabels([])
    gca().set_yticklabels([])
    gca().tick_params(length=0, which='both')
    dat_ax = gca()
    #gca().set_xlim(0, size_edg*sub_sh[1])
    #gca().set_ylim(0, size_edg*sub_sh[0])
    #gca().set_yticks([])
    subplot(122)
    gca().set_title("MODEL", fontsize=18)
    masked_mod_vals = full_im[~full_trust_im]
    full_im[~full_trust_im] = np.nan
    gca().set_facecolor('tomato')
    imshow(full_im, vmin=args.model_clim[0], vmax=args.model_clim[1], cmap='gray_r')#gnuplot')
    xt = np.arange(sub_sh[1], size_edg*sub_sh[1], sub_sh[1])-.5
    yt = np.arange(sub_sh[0], size_edg*sub_sh[0], sub_sh[0])-.5
    gca().set_xticks(xt)
    gca().set_yticks(yt)
    xt = gca().get_xticks()
    yt = gca().get_yticks()
    dx = abs(xt[1]-xt[0])/2
    dy = abs(yt[1]-yt[0])/2
    mid_x = list((xt[1:]+xt[:-1])//2)
    mid_y = list((yt[1:]+yt[:-1])//2)
    mid_x = [xt[0]-dx] + mid_x + [xt[-1]+dx]
    mid_y = [yt[0]-dy] + mid_y + [yt[-1]+dy]
    gca().set_xticks(mid_x, minor=True)
    gca().set_yticks(mid_y, minor=True)
    gca().grid(1, color='k', ls='--', which='major')
    gca().grid(1, color='k',alpha=0.5,lw=0.5, ls='--', which='minor')

    gca().set_xticklabels([])
    gca().set_yticklabels([])
    gca().tick_params(length=0, which='both')
    mod_ax = gca()
    display_fig = gcf()
    #gca().set_xlim(0, size_edg*sub_sh[1])
    #gca().set_ylim(0, size_edg*sub_sh[0])
    #gca().set_xticks([])
    #gca().set_yticks([])
    subplots_adjust(wspace=0, hspace=0, bottom=0.02, top=0.95, left=0.01, right=0.99)

    for x,y in labs["xy"]:
        mod_ax.text(x=x,y=y,s="", horizontalalignment="center", verticalalignment="center")
    for x,y in labs["xy"]:
        dat_ax.text(x=x,y=y,s="", horizontalalignment="center", verticalalignment="center")

    dat_contour_collections = []
    mod_contour_collections = []
    for kwargs in contour_kwargs:
        C=dat_ax.contour(kwargs['X'], kwargs['Y'], kwargs['dat_Z'], levels=kwargs['dat_levels'], cmap='jet', alpha=0.5)
        dat_contour_collections.append(C)
        C=mod_ax.contour(kwargs['X'], kwargs['Y'], kwargs['mod_Z'], levels=kwargs['mod_levels'], cmap='jet', alpha=0.5)
        mod_contour_collections.append(C)

    # make widgets to control colormap and things
    figure()
    from matplotlib.widgets import RadioButtons, TextBox, CheckButtons
    # control the cmap
    toggle_cmap_ax = plt.axes([0.01, 0.35, 0.15, 0.3])
    toggle_dat_clim_ax = plt.axes([0.4, 0.65, 0.25, 0.1])
    toggle_mod_clim_ax = plt.axes([0.4, 0.45, 0.25, 0.1])
    toggle_fontsize_ax = plt.axes([0.4, 0.25, 0.25, 0.1])
    toggle_label_ax = plt.axes([0.7, 0.35, 0.25, 0.5])
    toggle_contour_ax = plt.axes([0.7, 0.1, 0.25, 0.1])

    # Add a checkbutton for contours
    contour_check = RadioButtons(toggle_contour_ax, ("no contours", "contours") , active=0)
    def overlay_contours(label):
        if label=="no contours":
            for ax in [dat_ax, mod_ax]:
                for c in ax.collections:
                    c.remove()
        elif label=='contours':
            for C in dat_contour_collections:
                for c in C.collections:
                    dat_ax.add_collection(c)
            for C in mod_contour_collections:
                for c in C.collections:
                    mod_ax.add_collection(c)
        display_fig.canvas.draw_idle()
    overlay_contours("no contours")
    contour_check.on_clicked(overlay_contours)

    lab_buttons = RadioButtons(toggle_label_ax, ('none',r'$d$ ($\AA$)','hkl','panel','X (fast dim.)', 'Y (slow dim.)', r'$\sigma_Z$', 'max pixel'), active=0)
    def toggle_label(label):
        for i_txt,( mod_txt, dat_txt) in enumerate(zip(mod_ax.texts, dat_ax.texts)):
            if label=='none':
                s = ''
            elif label=='$d$ ($\AA$)':
                s = labs["res"][i_txt]
            elif label=='hkl':
                s = labs["hkl"][i_txt]
            elif label=="panel":
                s = labs["pfs"][i_txt].split(",")[0]
            elif label.startswith("X"):
                s = labs["pfs"][i_txt].split(",")[1]
            elif label.startswith("Y"):
                s = labs["pfs"][i_txt].split(",")[2]
            elif label=="$\sigma_Z$":
                s = labs["sigZ"][i_txt]
            elif label=="max pixel":
                s = labs["mod_max"][i_txt]
            mod_txt.set_text(s)

            if label=="max pixel" :
                s_dat = labs['dat_max'][i_txt]
                dat_txt.set_text(s_dat)
            else:
                dat_txt.set_text("")

            display_fig.canvas.draw_idle()
    lab_buttons.on_clicked(toggle_label)

    buttons = RadioButtons(toggle_cmap_ax, ('gray_r', 'gnuplot', 'cividis', 'hot', 'viridis'), active=0)
    def toggle_cmap(label):
        mod_ax.images[0].set_cmap(label)
        dat_ax.images[0].set_cmap(label)
        display_fig.canvas.draw_idle()
    buttons.on_clicked(toggle_cmap)


    def update_clim(text, ax):
        if "," in text:
            vmin, vmax = [float(l.strip()) for l in text.split(",")]
        else:
            vmin, vmax = [float(l.strip()) for l in text.split()]
        ax.images[0].set_clim(vmin,vmax)
        display_fig.canvas.draw_idle()

    def update_dat_clim(text):
        update_clim(text, dat_ax)

    def update_mod_clim(text):
        update_clim(text, mod_ax)

    init_dat_clim = "%d,%d" % dat_ax.images[0].get_clim()
    init_mod_clim = "%d,%d" % dat_ax.images[0].get_clim()
    dat_clim_box = TextBox(toggle_dat_clim_ax, 'Data vmin,vmax: ', initial=init_dat_clim)
    dat_clim_box.on_submit(update_dat_clim)
    mod_clim_box = TextBox(toggle_mod_clim_ax, 'Model vmin,vmax: ', initial=init_mod_clim)
    mod_clim_box.on_submit(update_mod_clim)

    def update_fontsize(text):
        new_fs, new_fc = text.split(",")
        new_fs = int(new_fs.strip())
        new_fc = new_fc.strip()
        for dat_txt, mod_txt in zip(dat_ax.texts, mod_ax.texts):
            dat_txt.set_fontsize(new_fs)
            dat_txt.set_color(new_fc)
            mod_txt.set_fontsize(new_fs)
            mod_txt.set_color(new_fc)
        display_fig.canvas.draw_idle()

    fs=mod_ax.texts[-1].get_fontsize()
    fc = mod_ax.texts[-1].get_color()
    fontsize_box = TextBox(toggle_fontsize_ax, 'Font size,color: ', initial="%d,%s"% (fs, fc))
    fontsize_box.on_submit(update_fontsize)

    for side in ['top', 'bottom', 'left', 'right']:
        mod_ax.spines[side].set(lw=2)
        dat_ax.spines[side].set(lw=2)

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
    pid = M.pids[i_h]
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

    try:
        segs = C.allsegs[-1][0]
        xcent, ycent = centroid_poly(segs[:,0], segs[:,1])
        ax0.plot( [xcent], [ycent], 'o', mec='b', ms=10, alpha=0.5)
        ax0.plot( segs[:,0], segs[:,1], '.', ms=3, color='k')
    except IndexError:
        pass

    v = np.max(np.abs(Z))
    v = 10
    ax2.imshow(Z, cmap='RdYlBu', vmin=-v, vmax=v)

    for AX in [ax0, ax1, ax2]:
        AX.grid(1, color='#777777', ls='--', lw=0.4)
        AX.set_xticklabels([])
        AX.set_yticklabels([])
    ax0.set_title("data %d: %d,%d,%d,%d"% (pid, x1,x2,y1,y2), pad=0, fontsize=8)
    ax1.set_title("model", pad=0, fontsize=8)
    ax2.set_title("Z (sigZ=%f)" % sigZ_val, pad=0, fontsize=8)

    FIG.suptitle("spot %d (%d/%d trusted)" % (i_h, trusted.sum(), trusted.size))
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
