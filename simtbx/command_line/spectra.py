from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.spectra

#TODO: add text entry boxes for other spectrum filter params, position entry boxes sensibly

from simtbx.diffBragg import hopper_utils
from pylab import *
import dxtbx
from simtbx.diffBragg.phil import philz, hopper_phil
from libtbx.phil import parse
from matplotlib.widgets import TextBox

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("image_file", type=str, help="path to a diffBragg modeler file (output from hopper, see the imgs folder in the outdir)")
parser.add_argument("--filt_freq", default=0.07, type=float)
parser.add_argument("--filt_order", default=3, type=float)
parser.add_argument("--tail", default=50, type=int)
parser.add_argument("--delta_en", default=0.5, type=float)
parser.add_argument("--skip", action="store_true")
args = parser.parse_args()

FIG,ax0 = subplots(nrows=1,ncols=1)
FIG.set_size_inches((5,3))

class P:

    def __init__(self, params, imgset, ax0, FIG):
        self.params = params
        self.imgset = imgset
        self.ax0 = ax0
        self.FIG= FIG

    def entry(self, text):
        delta_en = float(text)
        self.params.downsamp_spec.delta_en = delta_en
        print(delta_en)
        self.update_plot(self.imgset)

    def update_plot(self, i_img):
        raw_spec = self.imgset.get_spectrum(i_img)
        raw_en = raw_spec.get_energies_eV()
        raw_wt = raw_spec.get_weights()

        spec = hopper_utils.downsamp_spec_from_params(self.params, imgset=self.imgset, i_img=i_img)
        en, wt = map(np.array, zip(*spec))
        en = hopper_utils.utils.ENERGY_CONV / en

        self.ax0.clear()
        self.ax0.plot( raw_en, raw_wt, lw=2, label="raw spec (%d chan)" % len(raw_en))
        self.ax0.plot( en, wt, '--', lw=1, label="filt spec (%d chan)" % len(en))
        self.FIG.suptitle("image %d" % (i_img, ), fontsize=14)
        self.ax0.set_xlabel("Energy (eV)", fontsize=12)
        self.ax0.tick_params(labelsize=11)
        self.ax0.legend(prop={"size":11}, loc="upper right")
        draw()


def press(event):
    if event.key == 'right':
        FIG.loop_counter += 1
    elif event.key=="left":
        FIG.loop_counter = FIG.loop_counter -1
    FIG.loop_counter = max(FIG.loop_counter,0)
    FIG.loop_counter = min(FIG.loop_counter, FIG.nspots-1)

    if event.key=="escape":
        FIG.loop_counter = FIG.nspots

FIG.loop_counter = 0
loader = dxtbx.load(args.image_file)
imgset = loader.get_imageset(loader.get_image_file())
FIG.nspots = len(imgset)
FIG.canvas.mpl_connect('key_press_event', press)

phil_scope = parse(philz + hopper_phil)
params = phil_scope.fetch(sources=[phil_scope]).extract()

params.downsamp_spec.filt_freq = args.filt_freq
params.downsamp_spec.filt_order = args.filt_order
params.downsamp_spec.tail = args.tail
params.downsamp_spec.delta_en = args.delta_en
params.downsamp_spec.skip = args.skip

Pinst = P(params, imgset, ax0, FIG)
axbox = plt.axes([0.25, 0.75, 0.1, 0.075])
text_box = TextBox(axbox, 'delta_en', initial="%.2f" % args.delta_en)
text_box.on_submit(Pinst.entry)

while FIG.loop_counter < len(imgset):
    i_img = FIG.loop_counter

    Pinst.update_plot(i_img)

    waitforbuttonpress()
    i_img += 1

plt.close()
