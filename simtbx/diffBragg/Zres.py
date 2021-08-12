
from simtbx.diffBragg import utils
from scipy.ndimage import find_objects
import pandas
from pylab import *
Nrows = 7


def plot_Os(O, Nrows=7, prefix="some_plot", title=""):
    Nimg = np.ceil(len(O) / float(Nrows))
    for i_img in range(int(Nimg)):
        subO = O[i_img*Nrows: (i_img+1)*Nrows]
        plot_O(subO, Nrows=Nrows, outname=prefix+"_%d.png" % i_img, title=title)


def plot_O(O, Nrows=7, outname="dumb_plot.png", title="", tx=0.5, ty=0.98):
    close()
    close()
    close()
    fig, axs = subplots(nrows=Nrows, ncols=4, figsize=(8.5,11))
    ydims = [o[1][0].shape[0] for o in O]
    xdims = [o[1][0].shape[1] for o in O]
    print(ydims)
    print(xdims)
    min_ydim = min(ydims)
    min_xdim = min(xdims)
    Y = min_ydim
    X = min_xdim
    for i,(sigZ,b,trusted,d,reso) in enumerate(O):
        print(sigZ)
        A = axs[i]
        b[:,~trusted] = np.nan
        Z,model,data=b
        ydim,xdim = Z.shape

        for i_a, a in enumerate(A[:3]):
            a.clear()
            print(xdim/2.)
            print(ydim/2.)
            a.set_xticks([xdim/2,])
            a.set_yticks([ydim/2.])
            a.set_yticklabels([])
            a.set_xticklabels([])
            a.tick_params(length=0)
            #if i_a==2:
            #    a.grid(1, color='k', alpha=1, ls='--')
            #else:
            a.grid(1, color='#222222', alpha=1, ls='--')
        if i==0:
            #A[3].set_title("Z-score histogram")
            A[0].set_title("data (photons)")
            A[1].set_title("model (photons)")
            A[2].set_title("Z-score")
        #b[:,~c] = np.nan
        #Z,model, data = b
        #print(np.isnan(Z).sum())
        cm = 'RdYlBu'
        #cm2 = 'Greys'
        cm2 = matplotlib.cm.get_cmap('Greys')
        cm2.set_bad(color='red')

        A[0].imshow(b[2], cmap=cm2)
        cl = A[0].images[0].get_clim()
        A[1].imshow(b[1], cmap=cm2)
        A[1].images[0].set_clim(cl)
        A[2].imshow(b[0], cmap=cm)

        for ii,a in enumerate(A[:3]):
            cbar = fig.colorbar(a.images[0], ax=a, pad=0.05)
            cbar.ax.tick_params(labelsize=8, pad=2, length=2.5, direction='in', color="#444444")
            #if ii in [0,1]: cbar.ax.set_ylabel("photons", fontsize=8)
            #else:
            #    cbar.ax.set_ylabel("")
            a.set_ylim((ydim/2.-Y/2.-0.5, ydim/2.+Y/2.-0.5))
            a.set_xlim((xdim/2.-X/2.-0.5, xdim/2.+X/2.-0.5))

        vals = Z[trusted]
        xl = max(np.abs(vals))
        bins = np.linspace(-xl, xl, 50)
        binvals, bins = np.histogram(vals, bins=bins)
        A[3].clear()
        _, histbins, histpatch = A[3].hist(vals, bins=bins, color='#777777') #( bin_cent, bin_vals, '.')
        bin_cent= histbins[:-1]*.5 + histbins[1:] *.5

        gauss = utils.fit_gauss(binvals, bin_cent, amp=max(binvals))
        gauss = None
        if gauss is None:
            muZ = np.mean(vals)
            sigZ = np.std(vals)
            s='Population stats:\n  $\mu_Z=$%.3f\n  $\sigma_Z=$%.1f\nmin: %.1f\nmax: %.1f' % (muZ, sigZ, min(vals), max(vals))
        else:
            A[3].plot(bin_cent, gauss[2], color='tomato', label='$\sigma Z=$%.1f' % sigZ)#, bbox={'color': 'w'})
            ampZ, muZ, sigZ = gauss[0]
            s='Gaussian fit:\n  $\mu_Z=$%.3f\n  $\sigma_Z=$%.1f\nmin: %.1f\nmax: %.1f' % (muZ, sigZ, min(vals), max(vals))
        #fwhm = abs(sigZ)*2.355
        #A[3].text(x=0.85, y=0.75,s='$\fwhm_Z=$%.1f' % sigZ, transform=A[3].transAxes, bbox={'color':'w', 'alpha':0.9})

        A[3].text(x=0.85, y=0.35,s=s, fontsize=8, transform=A[3].transAxes, bbox={'color':'w', 'alpha': 1, 'ec':'k'})
        A[0].text(x=-.30, y=0.75,s='%.2f $\AA$' % (reso), transform=A[0].transAxes, bbox={'color':'w', 'alpha':0.9, 'ec':'k'}, fontsize=10)

        #A[3].legend()
        mx = max(binvals)
        A[3].set_yscale("log")
        A[3].set_yticks([1,10,100])
        A[3].set_yticklabels(['', '$10^1$','$10^2$'])
        #A[3].set_yticks([10,100])
        #A[3].set_yticklabels(['10', '100'])

        A[3].spines['top'].set_visible(False)
        A[3].spines['right'].set_visible(False)
        A[3].tick_params(axis='y', which='both', direction='in', labelsize=8, pad=-20)
        A[3].tick_params(axis='y', which='major', length=4, direction='in')
        A[3].tick_params(axis='y', which='minor', length=2, direction='in')
        A[3].tick_params(direction='in', pad=1,axis='x', length=3, labelsize=8)
        #A[3].set_facecolor('lightgray')
        xl = max(np.abs(vals))
        A[3].set_xlim(-xl, xl)
        cm = plt.cm.get_cmap('RdYlBu')
        col = bin_cent - min(bin_cent)
        col /= max(col)
        for c,p in zip(col, histpatch):
            plt.setp(p, 'facecolor', cm(c))
            plt.setp(p, 'ec', '#444444')
            plt.setp(p, 'lw', 0.5)
        A[2].images[0].set_clim(-xl, xl)
        if i == 0:
            #A[3].set_ylabel("pixel count", labelpad=-31, fontsize=8)
            A[3].set_ylabel("pixel count", labelpad=1, fontsize=8)
        #A[3].set_yscale("log")
        #A[3].grid(1, alpha=0.5)
    axs[-1][3].set_xlabel("Z-score")
    axs[0][3].set_title("Z-score histogram", pad=7)
    for i in range(len(O), len(axs)):
        for a in axs[i]:
            a.remove()
    suptitle(title, x=tx, y=ty,horizontalalignment='left', fontsize=18, fontweight='bold')
    subplots_adjust(left=0.05, right=0.89, bottom=0.04, top=0.9)
    savefig(outname)

    draw()
    pause(0.1)


def main():
    ALL_OUT = []
    nbins = 10
    img_sh = 256, 254, 254

    df = pandas.read_pickle("some_models_from_8_stage2_2_iter1000.pkl")
    res_bins = [s[0] for s in np.array_split(np.sort(df.res), 10)] + [df.res.max()]
    df['res_id'] = np.digitize(df.res, res_bins)
    for i in range(nbins):
        i = i+1
        df_res = df.query("res_id==%d" % i)
        i_fcells = df_res.i_fcell.unique()
        OUT = []
        for iii,i_fcell in enumerate(i_fcells[:10]):
            if i_fcell == 0:
                print("skipping 0")
            img = np.zeros(img_sh, int)
            img2 =np.zeros((3,)+img_sh)
            trus = np.zeros(img_sh, bool)
            df_i = df_res.query("i_fcell==%d" % i_fcell)
            npan = len(df_i.p.unique())
            print("i_f=%d, npan=%d" % (i_fcell, npan))
            img[df_i.p.astype(int), df_i.s.astype(int), df_i.f.astype(int)] = iii+1
            objs = [find_objects(pimg) for pimg in img]
            img2[0,df_i.p, df_i.s, df_i.f] = df_i.Zscore
            img2[1,df_i.p, df_i.s, df_i.f] = df_i.model
            img2[2,df_i.p, df_i.s, df_i.f] = df_i.data
            trus[df_i.p, df_i.s, df_i.f] = df_i.trust
            spot_res = df_i.res.values[0]
            for pid, O in enumerate(objs):
                if not O:
                    continue
                for oo in O:
                    if oo is None:
                        continue
                    S,F = oo
                    Z = img2[:,pid,S,F]
                    T = trus[pid,S,F]
                    sigZ = Z[0][T].std()
                    OUT.append((sigZ, Z,T,i, spot_res))
                    print("resbin %d, i_fcell=%d, obj=%d, sigZ=%f" % (i, i_fcell, pid, sigZ))
        ALL_OUT.append(OUT)


def main2(best=True):
    nresbin = 10
    df = pandas.read_pickle("500pkls.pkl")

    res_bins = [s[0] for s in np.array_split(np.sort(df.res), nresbin)] + [df.res.max()]
    df['res_id'] = np.digitize(df.res, res_bins)

    results = np.load("results.npy")

    p,i_fcell, shot, sigZ, rng  = results.T
    p = p.astype(int)
    i_fcell = i_fcell.astype(int)
    shot = shot.astype(int)

    keys = list(zip(p,i_fcell, shot))
    results_sort = sorted(zip(sigZ, rng, keys))

    res_id_from_i_fcell = {i_fcell:res_id for i_fcell, res_id in zip(df.i_fcell, df.res_id)}
    ALL_O = {i:[] for i in range(nresbin)}
    R = results_sort
    if not best:
        R = R[::-1]
    for sigZ, rng, k in R:
        pid, i_fcell, shot = k
        if not best and pid in [234, 245]:
            continue
        res_id = res_id_from_i_fcell[i_fcell] - 1
        if res_id > nresbin-1:
            continue
        if len(ALL_O[res_id]) == 7:
            continue
        shot = str(shot)
        df_g = df.query("p==%d" % pid).query("i_fcell==%d" % i_fcell).query("shot=='%s'" % shot)
        if (~df_g.trust).sum() > 10:
            continue
        if len(df_g) != 169:
            continue
        if len(ALL_O[res_id]) < 7:
            X = df_g.f - df_g.f.min()
            Y = df_g.s - df_g.s.min()
            xdim = X.max()+1
            ydim = Y.max()+1
            imgs = np.zeros( (3, ydim, xdim))
            imgs[0,Y,X] = df_g.Zscore
            imgs[1,Y,X] = df_g.model
            imgs[2,Y,X] = df_g.data
            trust = np.zeros((ydim,xdim), bool)
            trust[Y,X] = df_g.trust


            ALL_O[res_id].append((sigZ, imgs, trust, k, df_g.res.values[0]))
        if all([len(ALL_O[ii])==7 for ii in range(nresbin)]):
            break
        print(res_id, len(ALL_O[res_id]), sigZ)


    for i_res in range(nresbin):
        name = "Zplots_2/best_resbin%d.png" % i_res
        if not best:
            name = "Zplots_2/worst_resbin%d.png" % i_res
        acc = "high" if best else "low"
        title = r'%s accuracy; resolution range %.2f - %.2f $(\mathrm{\AA})$' % (acc, res_bins[i_res], res_bins[i_res+1])
        plot_O(sorted(ALL_O[i_res], key=lambda x: x[-1]), outname=name, title=title, tx=0.075)


if __name__ == "__main__":
    main2()
    #main2(best=False)
