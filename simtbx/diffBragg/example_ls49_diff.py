from simtbx.diffBragg.refine_LS49 import RefineLS49Rot

# data from a single image
data = None # make me not None

# get the scale factor
RLS = RefineLS49Rot(data=data, refine_bg_planes=False, refine_angles=False)
init_scale = RLS.x[-1]
RLS = RefineLS49Rot(data=data, refine_bg_planes=True, refine_angles=True, 
        init_scale=init_scale, refine_scale=True) #NOTE: consider not refine scale again

dxcryst_refined = RLS.update_crystal_model()
print dxcryst_refined.get_A()

