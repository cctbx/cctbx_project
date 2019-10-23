
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    has_mpi = True
except ImportError:
    rank = 0
    size = 1
    has_mpi = False

if rank == 0:
    from numpy import mean, unique, log
    from numpy import all as np_all
    from numpy.linalg import norm
    from sys.stdout import flush as stdout_flush

    from simtbx.diffBragg.refiners import BreakToUseCurvatures
    from scitbx.array_family import flex
    flex_double = flex.double
    from scitbx.matrix import col
    from simtbx.diffBragg.refiners import PixelRefinement

else:
    mean = unique = log = np_all = norm = None
    stdout_flush = None
    BreakToUseCurvatures = None
    flex_double = None
    col = None
    PixelRefinement = None

if has_mpi:
    mean = comm.broadcast(mean, root=0)
    stdout_flush = comm.broadcast(stdout_flush)
    BreakToUseCurvatures = comm.broadcast(BreakToUseCurvatures)
    flex_double = comm.broadcast(flex_double)
    col = comm.broadcast(col)
    PixelRefinement = comm.broadcast(PixelRefinement)

if rank == 0:
    import pylab as plt
    from IPython import embed


class GlobalFat(PixelRefinement):

    def __init__(self, n_global_params, n_local_params, local_idx_start,
                 shot_ucell_managers, shot_rois, shot_nanoBragg_rois,
                 shot_roi_imgs, shot_spectra, shot_crystal_GTs,
                 shot_crystal_models, shot_xrel, shot_yrel, shot_abc_inits,
                 global_param_idx_start,
                 shot_panel_ids, init_gain=1, init_scale=1):
        super(GlobalFat, self).__init__()

        # dictionaries whose keys are the shot indices
        self.UCELL_MAN = shot_ucell_managers
        # cache shot ids and make sure they are identical in all other input dicts
        self.shot_ids = sorted(shot_ucell_managers.keys())
        self.n_shots = len(self.shot_ids)
        # sanity check: no repeats of the same shot
        assert len(self.shot_ids) == len(set(self.shot_ids))

        self.ROIS = self._check_keys(shot_rois)
        self.NANOBRAGG_ROIS = self._check_keys(shot_nanoBragg_rois)
        self.ROI_IMGS = self._check_keys(shot_roi_imgs)
        self.SPECTRA = self._check_keys(shot_spectra)
        self.CRYSTAL_GT = self._check_keys(shot_crystal_GTs)
        self.CRYSTAL_MODELS = self._check_keys(shot_crystal_models)
        self.XREL = self._check_keys(shot_xrel)
        self.YREL = self._check_keys(shot_yrel)
        self.ABC_INIT = self._check_keys(shot_abc_inits)
        self.PANEL_IDS = self._check_keys(shot_panel_ids)

        # Total number of parameters in the MPI world
        self.n_global_params = n_global_params

        # total number of local parameters
        self.n_local_params = n_local_params

        # here are the indices of the local parameters in the global paramter arrays
        self.local_idx_start = local_idx_start

        self.calc_func = True  # NOTE: leave True, debug flag from older code
        self.multi_panel = True  # we are multi panel
        self.f_vals = []  # store the functional over time

        # start with the first shot
        self._i_shot = self.shot_ids[0]

        # These are the per-shot parameters
        self.n_rot_param = 3
        self.n_ucell_param = len(self.UCELL_MAN[self._i_shot].variables)
        self.n_ncells_param = 1
        self.n_spot_scale_param = 1
        self.n_per_shot_params = self.n_rot_param + self.n_ucell_param + self.n_ncells_param + self.n_spot_scale_param

        assert self.n_per_shot_params == self.n_local_params

        self._ncells_id = 9  # diffBragg specific
        self._originZ_id = 10  # diffBragg specific
        self._init_scale = init_scale
        self._init_gain = init_gain
        self.num_positive_curvatures = 0
        self._panel_id = None

        self.pid_from_idx = None
        self.idx_from_pid = None

        # where the global parameters being , initially just gain and detector distance
        self.global_param_idx_start = global_param_idx_start

        self.a = self.b = self.c = None

    def _setup_plots(self):
        if rank == 0:
            if self.plot_images:
                if self.plot_residuals:
                    self.fig = plt.figure()
                    self.ax = self.fig.gca(projection='3d')
                    self.ax.set_yticklabels([])
                    self.ax.set_xticklabels([])
                    self.ax.set_zticklabels([])
                    self.ax.set_zlabel("model residual")
                    self.ax.set_facecolor("gray")
                else:
                    self.fig, (self.ax1, self.ax2) = plt.subplots(nrows=1, ncols=2)
                    self.ax1.imshow([[0, 1, 1], [0, 1, 2]])
                    self.ax2.imshow([[0, 1, 1], [0, 1, 2]])

    @property
    def n(self):
        """LBFGS property"""
        return self.n_global_params

    @property
    def x(self):
        """LBFGS parameter array"""
        return self._x

    @x.setter
    def x(self, val):
        self._x = val

    def _check_keys(self, shot_dict):
        """checks that the dictionary keys are the same"""
        if not sorted(shot_dict.keys()) == self.shot_ids:
            raise KeyError("input data funky, check GlobalFat inputs")

    def _evaluate_averageI(self):
        """model_Lambda means expected intensity in the pixel"""
        self.model_Lambda = \
            self.gain_fac*self.gain_fac*(self.tilt_plane + self.scale_fac * self.scale_fac * self.model_bragg_spots)

    def _evaluate_log_averageI(self):
        # fix log(x<=0)
        try:
            self.log_Lambda = log(self.model_Lambda)
        except FloatingPointError:
            pass
        #if any((self.model_Lambda <= 0).ravel()):
        #    print("\n<><><><><><><><>\n\tWARNING: NEGATIVE INTENSITY IN MODEL!!!!!!!!!\n<><><><><><><><><>\n")
        #    raise ValueError("model of Bragg spots cannot have negative intensities...")
        self.log_Lambda[self.model_Lambda <= 0] = 0  # FIXME with Gaussian model

    def _set_background_plane(self, i_spot):
        xr = self.XREL[self._i_shot][i_spot]
        yr = self.YREL[self._i_shot][i_spot]
        self.tilt_plane = xr*self.a + yr*self.b + self.c

    def _setup(self):
        # Here we go!

        # Make a mapping of panel id to parameter index and backwards
        self.pid_from_idx = {}
        self.idx_from_pid = {}

        # Make te global sized parameter array, though here we only update the local portion
        self.x = flex_double(self.n_global_params)

        self.rotX_xpos = {}
        self.rotY_xpos = {}
        self.rotZ_xpos = {}
        self.ucell_xstart = {}
        self.ncells_xpos = {}
        self.spot_scale_xpos = {}
        self.n_panels = {}

        for i_shot in self.shot_ids:

            self.pid_from_idx[i_shot] = {i: pid for i, pid in enumerate(unique(self.panel_ids))}
            self.idx_from_pid[i_shot] = {pid: i for i, pid in enumerate(unique(self.panel_ids))}
            self.n_panels[i_shot] = len(self.pid_from_idx)

            self.rotX_xpos[i_shot] = self.local_idx_start + i_shot*self.n_per_shot_params
            self.rotY_xpos[i_shot] = self.rotX_xpos[i_shot] + 1
            self.rotZ_xpos[i_shot] = self.rotY_xpos[i_shot] + 1

            self.x[self.rotX_xpos] = 0
            self.x[self.rotY_xpos] = 0
            self.x[self.rotZ_xpos] = 0

            self.ucell_xstart[i_shot] = self.rotZ_xpos[i_shot] + 1
            # populate the x-array with initial values
            for i_uc in range(self.n_ucell_param):
                self.x[self.ucell_xstart[i_shot] + i_uc] = self.UCELL_MAN[i_shot].variables[i_uc]

            # put in Ncells abc estimate
            self.ncells_xpos[i_shot] = self.ucell_xstart[i_shot] + self.n_ucell_param
            self.x[self.ncells_xpos[i_shot]] = self.S.crystal.Ncells_abc[0]  # NOTE: each shot gets own starting Ncells

            self.spot_scale_xpos[i_shot] = self.ncells_xpos[i_shot] + 1
            self.x[self.spot_scale_xpos[i_shot]] = self._init_scale  # TODO: each shot gets own starting scale factor

        self.originZ_xpos = self.global_param_idx_start
        self.gain_xpos = self.n_global_params - 1
        if rank == 0:
            # put in estimates for origin vectors
            # TODO: refine at the different hierarchy
            # get te first Z coordinate for now..
            self.x[self.originZ_xpos] = self.S.detector[0].get_local_origin()[2]  # NOTE maybe just origin instead?
            self.x[self.gain_pos] = self._init_gain  # gain factor

            # TODO per panel gain corretion / pedestal correction?
            #n_panels = len(self.S.detector)
            #self.origin_xstart = self.global_param_idx_start
            #for i_pan in range(n_panels):
            #    pid = self.pid_from_idx[i_pan]
            #    self.x[self.origin_xstart + i_pan] = self.S.detector[pid].get_local_origin()[2]
            # lastly, the panel gain correction factor
            #for i_pan in range(n_panels):
            #    self.x[-1] = self._init_gain
            #self.x[-1] = self._init_scale  # initial scale factor

        # reduce then broadcast self.x
        # TODO do I work properly ?
        self.x = comm.reduce(self.x, MPI.SUM, root=0)
        self.x = comm.broadcast(self.x, root=0)

        # setup the diffBragg instance
        self.D = self.S.D

        if self.refine_Umatrix:
            if self.refine_rotX:
                self.D.refine(0)  # rotX
            if self.refine_rotY:
                self.D.refine(1)  # rotY
            if self.refine_rotZ:
                self.D.refine(2)  # rotZ
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self.D.refine(i + 3)  # unit cell params

        if self.refine_ncells:
            self.D.refine(self._ncells_id)
        if self.refine_detdist:
            self.D.refine(self._originZ_id)
        self.D.initialize_managers()

    def _send_gradients_to_derivative_managers(self):
        """Needs to be called once each time the orientation is updated"""
        for i in range(self.n_ucell_param):
            self.D.set_ucell_derivative_matrix(
                i + 3,
                self.UCELL_MAN[self._i_shot].derivative_matrices[i])
            if self.calc_curvatures:
                self.D.set_ucell_second_derivative_matrix(
                    i + 3, self.UCELL_MAN[self._i_shot].second_derivative_matrices[i])

    def _run_diffBragg_current(self, i_spot):
        """needs to be called each time the ROI is changed"""
        self.D.region_of_interest = self.NANOBRAGG_ROIS[self._i_shot][i_spot]
        self.D.add_diffBragg_spots()

    def _update_rotXYZ(self):
        if self.refine_rotX:
            self.D.set_value(0, self.x[self.rotX_xpos[self._i_shot]])
        if self.refine_rotY:
            self.D.set_value(1, self.x[self.rotY_xpos[self._i_shot]])
        if self.refine_rotZ:
            self.D.set_value(2, self.x[self.rotZ_xpos[self._i_shot]])

    def _update_ncells(self):
        self.D.set_value(self._ncells_id, self.x[self.ncells_xpos[self._i_shot]])

    def _update_dxtbx_detector(self):

        det = self.S.detector
        self.S.panel_id = self._panel_id
        # TODO: select hierarchy level at this point
        # NOTE: what does fast-axis and slow-axis mean
        # for the different hierarchy levels?
        node = det[self._panel_id]
        orig = node.get_local_origin()
        new_originZ = self.x[self.originZ_xpos]
        new_local_origin = orig[0], orig[1], new_originZ
        node.set_local_frame(node.get_local_fast_axis(),
                             node.get_local_slow_axis(),
                             new_local_origin)
        self.S.detector = det  # TODO  update the sim_data detector? maybe not necessary after this point
        self.D.update_dxtbx_geoms(det, self.S.beam.nanoBragg_constructor_beam, self._panel_id)

    def _extract_pixel_data(self):
        self.rot_deriv = [0, 0, 0]
        self.rot_second_deriv = [0, 0, 0]
        if self.refine_Umatrix:
            if self.refine_rotX:
                self.rot_deriv[0] = self.D.get_derivative_pixels(0).as_numpy_array()
                if self.calc_curvatures:
                    self.rot_second_deriv[0] = self.D.get_second_derivative_pixels(0).as_numpy_array()
            if self.refine_rotY:
                self.rot_deriv[1] = self.D.get_derivative_pixels(1).as_numpy_array()
                if self.calc_curvatures:
                    self.rot_second_deriv[1] = self.D.get_second_derivative_pixels(1).as_numpy_array()
            if self.refine_rotZ:
                self.rot_deriv[2] = self.D.get_derivative_pixels(2).as_numpy_array()
                if self.calc_curvatures:
                    self.rot_second_deriv[2] = self.D.get_second_derivative_pixels(2).as_numpy_array()

        self.ucell_derivatives = [0]*self.n_ucell_param
        self.ucell_second_derivatives = [0]*self.n_ucell_param
        if self.refine_Bmatrix:
            for i in range(self.n_ucell_param):
                self.ucell_derivatives[i] = self.D.get_derivative_pixels(3 + i).as_numpy_array()
                if self.calc_curvatures:
                    self.ucell_second_derivatives[i] = self.D.get_second_derivative_pixels(3 + i).as_numpy_array()

        self.ncells_deriv = self.detdist_deriv = 0
        self.ncells_second_deriv = self.detdist_second_deriv = 0
        if self.refine_ncells:
            self.ncells_deriv = self.D.get_derivative_pixels(self._ncells_id).as_numpy_array()
            if self.calc_curvatures:
                self.ncells_second_deriv = self.D.get_second_derivative_pixels(self._ncells_id).as_numpy_array()
        if self.refine_detdist:
            self.detdist_deriv = self.D.get_derivative_pixels(self._originZ_id).as_numpy_array()
            if self.calc_curvatures:
                self.detdist_second_deriv = self.D.get_second_derivative_pixels(self._originZ_id).as_numpy_array()

        self.model_bragg_spots = self.D.raw_pixels_roi.as_numpy_array()

    def _unpack_bgplane_params(self, i_spot):
        self.a, self.b, self.c = self.ABC_INIT[self._i_shot][i_spot]

    def _update_ucell(self):
        _s = slice(self.ucell_xstart[self._i_shot], self.ucell_xstart[self._i_shot] + self.n_ucell_param, 1)
        pars = list(self.x[_s])
        self.ucell_manager.variables = pars
        self._send_gradients_to_derivative_managers()
        self.D.Bmatrix = self.ucell_manager.B_recipspace

    def _update_umatrix(self):
        self.D.Umatrix = self.CRYSTAL_MODELS[self._i_shot].get_U()

    def _update_beams(self):
        # sim_data instance has a nanoBragg beam object, which takes spectra and converts to nanoBragg xray_beams
        self.S.beam.spectra = self.SPECTRA[self._i_shot]
        self.D.xray_beams = self.S.beam.xray_beams

    def compute_functional_and_gradients(self):
        if self.calc_func:
            if self.verbose:
                if self.use_curvatures:
                    print("Compute functional and gradients Iter %d (Using Curvatures)\n<><><><><><><><><><><><><>"
                          % (self.iterations+1))
                else:
                    print("Compute functional and gradients Iter %d PosCurva %d\n<><><><><><><><><><><><><>"
                          % (self.iterations+1, self.num_positive_curvatures))
            f = 0
            g = flex_double(self.n)
            if self.calc_curvatures:
                self.curv = flex_double(self.n)

            self.gain_fac = self.x[self.gain_xpos]
            G2 = self.gain_fac**2

            for self._i_shot in self.shot_ids:
                self.scale_fac = self.x[self.spot_scale_xpos[self._i_shot]]
                S2 = self.scale_fac**2
                #TODO: Omatrix update? All crystal models here should have the same to_primitive operation, ideally
                self._update_beams()
                self._update_umatrix()
                self._update_ucell()
                self._update_ncells()
                self._update_rotXYZ()
                for i_spot in range(self.n_spots):
                    self._panel_id = self.PANEL_IDS[self._i_shot][i_spot]
                    if self.verbose:
                        print "\rdiffBragg: img %d/%d; spot %d/%d; panel %d" \
                              % (self._i_shot+1, self.n_shots, i_spot+1, self.n_spots, self._panel_id),
                        stdout_flush()

                    self._update_dxtbx_detector()
                    self._run_diffBragg_current(i_spot)
                    self._unpack_bgplane_params(i_spot)
                    self._set_background_plane(i_spot)
                    self._extract_pixel_data()
                    self._evaluate_averageI()
                    self._evaluate_log_averageI()

                    Imeas = self.ROI_IMGS[self._i_shot][i_spot]
                    f += (self.model_Lambda - Imeas * self.log_Lambda).sum()
                    one_over_Lambda = 1. / self.model_Lambda
                    one_minus_k_over_Lambda = (1. - Imeas * one_over_Lambda)
                    if self.calc_curvatures:
                        k_over_squared_Lambda = Imeas * one_over_Lambda * one_over_Lambda

                    if self.plot_images and self.iterations % self.plot_stride == 0:
                        xr = self.XREL[self._i_shot][i_spot]  # fast scan pixels
                        yr = self.YREL[self._i_shot][i_spot]  # slow scan pixels
                        if self.plot_residuals:
                            self.ax.clear()
                            residual = self.model_Lambda - Imeas
                            if i_spot == 0:
                                x = residual.max()
                            else:
                                x = mean([x, residual.max()])

                            self.ax.plot_surface(xr, yr, residual, rstride=2, cstride=2, alpha=0.3, cmap='coolwarm')
                            self.ax.contour(xr, yr, residual, zdir='z', offset=-x, cmap='coolwarm')
                            self.ax.set_yticks(range(yr.min(), yr.max()))
                            self.ax.set_xticks(range(xr.min(), xr.max()))
                            self.ax.set_xticklabels([])
                            self.ax.set_yticklabels([])
                            self.ax.set_zlim(-x, x)
                            self.ax.set_title("residual (photons)")
                        else:
                            m = Imeas[Imeas > 1e-9].mean()
                            s = Imeas[Imeas > 1e-9].std()
                            vmax = m+5*s
                            vmin = m-s
                            m2 = self.model_Lambda.mean()
                            s2 = self.model_Lambda.std()
                            self.ax1.images[0].set_data(self.model_Lambda)
                            self.ax1.images[0].set_clim(vmin, vmax)
                            self.ax2.images[0].set_data(Imeas)
                            self.ax2.images[0].set_clim(vmin, vmax)
                        plt.suptitle("Iterations = %d, image %d / %d"
                                     % (self.iterations, i_spot+1, self.n_spots))
                        self.fig.canvas.draw()
                        plt.pause(.02)
                    if self.refine_Umatrix:
                        x_positions = [self.rotX_xpos[self._i_shot],
                                       self.rotY_xpos[self._i_shot],
                                       self.rotZ_xpos[self._i_shot]]
                        for ii, xpos in enumerate(x_positions):
                            d = S2*G2*self.rot_deriv[ii]
                            g[xpos] += (one_minus_k_over_Lambda * d).sum()
                            if self.calc_curvatures:
                                d2 = S2*G2*self.rot_second_deriv[ii]
                                cc = d2*one_minus_k_over_Lambda + d*d*k_over_squared_Lambda
                                self.curv[xpos] += cc.sum()

                    if self.refine_Bmatrix:
                        # unit cell derivative
                        for i_ucell_p in range(self.n_ucell_param):
                            xpos = self.ucell_xstart[self._i_shot] + i_ucell_p
                            d = S2*G2*self.ucell_derivatives[i_ucell_p]
                            g[xpos] += (one_minus_k_over_Lambda * d).sum()

                            if self.calc_curvatures:
                                d2 = S2*G2*self.ucell_second_derivatives[i_ucell_p]
                                cc = d2*one_minus_k_over_Lambda + d*d*k_over_squared_Lambda
                                self.curv[xpos] += cc.sum()

                    if self.refine_ncells:
                        d = self.ncells_deriv
                        xpos = self.ncells_xpos[self._i_shot]
                        g[xpos] += (S2*G2*d*one_minus_k_over_Lambda).sum()
                        if self.calc_curvatures:
                            d2 = S2*G2*self.ncells_second_deriv
                            cc = d2*one_minus_k_over_Lambda + d*d*k_over_squared_Lambda
                            self.curv[xpos] += cc.sum()

                    if self.refine_detdist:
                        raise NotImplementedError("Cannot refined detdist (yet...)")
                        #if self.calc_curvatures:
                        #    raise NotImplementedError("Cannot use curvatures and refine detdist (yet...)")
                        #origin_xpos = self.origin_xstart + self.idx_from_pid[self._panel_id]
                        #g[origin_xpos] += (S2*G2*self.detdist_deriv*one_minus_k_over_Lambda).sum()

                    if self.refine_crystal_scale:
                        d = G2*2*self.scale_fac*self.model_bragg_spots
                        xpos = self.spot_scale_xpos[self._i_shot]
                        g[xpos] += (d*one_minus_k_over_Lambda).sum()
                        if self.calc_curvatures:
                            d2 = d / self.scale_fac
                            self.curv[xpos] += (d2*one_minus_k_over_Lambda + d*d*k_over_squared_Lambda).sum()

                    if self.refine_gain_fac:
                        d = 2*self.gain_fac*(self.tilt_plane + S2*self.model_bragg_spots)
                        g[self.gain_xpos] += (d*one_minus_k_over_Lambda).sum()
                        if self.calc_curvatures:
                            d2 = d / self.gain_fac
                            self.curv[self.gain_xpos] += (d2*one_minus_k_over_Lambda + d*d*k_over_squared_Lambda).sum()

            # reduce the broadcast summed results:
            f = comm.reduce(f, MPI.SUM, root=0)
            f = comm.broadcast(f, root=0)
            g = comm.reduce(g, MPI.SUM, root=0)
            g = comm.broadcast(g, root=0)
            if self.calc_curvatures:
                self.curv = comm.reduce(self.curv, MPI.SUM, root=0)
                self.curv = comm.broadcast(self.curv, root=0)

            if self.calc_curvatures and not self.use_curvatures:
                if np_all(self.curv.as_numpy_array() >= 0):
                    self.num_positive_curvatures += 1
                else:
                    self.num_positive_curvatures = 0

            self._f = f
            self._g = g
            self.D.raw_pixels *= 0
            gnorm = norm(g)
            if self.verbose:
                self.print_step("LBFGS stp", f)
                self.print_step_grads("LBFGS GRADS", gnorm )
            self.iterations += 1
            self.f_vals.append(f)

            if self.calc_curvatures and not self.use_curvatures:
                if self.num_positive_curvatures == self.use_curvatures_threshold:
                    raise BreakToUseCurvatures

        return self._f, self._g

    def _mpi_reduce_broadcast(self, var):
        # TODO: try me out
        var = comm.reduce(var, MPI.SUM, root=0)
        var = comm.broadcast(var, root=0)
        return var

    def print_step(self, message, target):
        names = self.ucell_manager.variable_names
        vals = self.ucell_manager.variables
        ucell_labels = []
        for n, v in zip(names, vals):
            ucell_labels.append('%s=%+2.7g' % (n, v))
        rotX = self.x[self.rotX_xpos]
        rotY = self.x[self.rotY_xpos]
        rotZ = self.x[self.rotZ_xpos]
        rot_labels = ["rotX=%+3.7g" % rotX, "rotY=%+3.7g" % rotY, "rotZ=%+3.4g" % rotZ]
        print("%s: residual=%3.8g, ncells=%f, detdist=%3.8g, gain=%3.4g, scale=%4.8g"
              % (message, target, self.x[self.ncells_xpos], self.x[self.origin_xstart], self.x[-2]**2, self.x[-1]**2))
        print ("Ucell: %s *** Missets: %s" %
               (", ".join(ucell_labels),  ", ".join(rot_labels)))
        print("\n")

    def print_step_grads(self, message, target):
        names = self.ucell_manager.variable_names
        vals = self.ucell_manager.variables
        ucell_labels = []
        for i,(n, v) in enumerate(zip(names, vals)):
            grad = self._g[self.ucell_xstart+i]
            ucell_labels.append('G%s=%+2.7g' % (n, grad))
        rotX = self._g[self.rotX_xpos]
        rotY = self._g[self.rotY_xpos]
        rotZ = self._g[self.rotZ_xpos]
        rot_labels = ["GrotX=%+3.7g" % rotX, "GrotY=%+3.7g" % rotY, "GrotZ=%+3.4g" % rotZ]
        xnorm = norm(self.x)
        print("%s: |x|=%f, |g|=%f, Gncells=%f, Gdetdist=%3.8g, Ggain=%3.4g, Gscale=%4.8g"
              % (message, xnorm, target, self._g[self.ncells_xpos], self._g[self.origin_xstart], self._g[-2], self._g[-1]))
        print ("GUcell: %s *** GMissets: %s" %
               (", ".join(ucell_labels),  ", ".join(rot_labels)))
        print("\n")

    def get_refined_Bmatrix(self, i_shot):
        return self.UCELL_MAN[i_shot].B_recipspace

    def curvatures(self):
        return self.curv

    def get_correction_misset(self, i_shot, as_axis_angle_deg=False, angles=None):
        """
        return the current state of the perturbation matrix
        :return: scitbx.matrix sqr
        """
        if angles is None:
            anglesXYZ = self.x[self.rotX_xpos[i_shot]], self.x[self.rotY_xpos[i_shot]], self.x[self.rotZ_xpos[i_shot]]
        x = col((-1, 0, 0))
        y = col((0, -1, 0))
        z = col((0, 0, -1))
        RX = x.axis_and_angle_as_r3_rotation_matrix(anglesXYZ[0], deg=False)
        RY = y.axis_and_angle_as_r3_rotation_matrix(anglesXYZ[1], deg=False)
        RZ = z.axis_and_angle_as_r3_rotation_matrix(anglesXYZ[2], deg=False)
        M = RX*RY*RZ
        if as_axis_angle_deg:
            q = M.r3_rotation_matrix_as_unit_quaternion()
            rot_ang, rot_ax = q.unit_quaternion_as_axis_and_angle(deg=True)
            return rot_ang, rot_ax
        else:
            return M
