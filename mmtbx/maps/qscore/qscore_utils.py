import math
from collections import defaultdict
from itertools import chain

import numpy as np
import numpy.ma as ma
import cctbx
from cctbx.array_family import flex

from mmtbx.maps.qscore.flex_utils import *


# def sphere_points(ctr, rad, N):
#     if ctr.ndim==1:
#         ctr = ctr[None,:]
#     h = -1.0 + (2.0 * np.arange(N) / float(N-1))[:, np.newaxis]
#     phis = np.arccos(h)
#     thetas = np.zeros_like(phis)
#     thetas[1:-1, :] = (3.6 / np.sqrt(N * (1.0 - h[1:-1]**2))) % (2 * np.pi)
#     thetas = np.cumsum(thetas, axis=0)

#     x = np.sin(phis) * np.cos(thetas)
#     y = np.sin(phis) * np.sin(thetas)
#     z = np.cos(phis)

#     # Stack x, y, z to form points and multiply by rad
#     points = rad * np.stack([x, y, z], axis=-1)

#     # Reshape points to (1, N, 3)
#     points = points.reshape(1, N, 3)

#     # Add center coordinates to all points
#     # ctr shape: (M, 3), points shape: (1, N, 3)
#     # Resultant shape: (M, N, 3)
#     pts = ctr[:, np.newaxis, :] + points


#     return pts
  
# def fisher(r):
#     z = 0.5 * np.log((1 + r) / (1 - r))
#     return z
  
# def points_for_density(radii, density,min_value=1):
#     return max(min_value,np.round(density * 4 * np.pi * radii**2).astype(int))





# def standardize_selection_np(selection,model=None,n_atoms=None):
#     # can get bool,int,or string. If int or string must pass model
#     # Returns numpy bool selection

#     if isinstance(selection,str):
#         assert model is not None,"Must also pass model to use string selection"
#         selection = model.selection(selection) # bool
#     elif isinstance(selection,np.ndarray):
#         if selection.dtype == int:
#             assert model is not None or n_atoms is not None,"Must also pass model to use int selection"
#             if model is not None:
#               n_atoms = model.get_number_of_atoms()
#             selection_bool = np.full(n_atoms,False)
#             selection_bool[selection] = True
#             selection = selection_bool
#         elif selection.dtype == bool:
#             pass
#         else:
#             assert False, "Unable to interpret selection"
#     else:
#         assert False, "Unable to interpret selection"

#     return selection

def trilinear_interpolation(voxel_grid, coords, voxel_size=None, offset=None):
    assert voxel_size is not None, "Provide voxel size as an array or single value"
    
    # Apply offset if provided
    if offset is not None:
        coords = coords - offset

    # Transform coordinates to voxel grid index space
    index_coords = coords / voxel_size

    # Split the index_coords array into three arrays: x, y, and z
    x, y, z = index_coords.T

    # Truncate to integer values
    x0, y0, z0 = np.floor([x, y, z]).astype(int)
    x1, y1, z1 = np.ceil([x, y, z]).astype(int)

    # Ensure indices are within grid boundaries
    x0, y0, z0 = np.clip([x0, y0, z0], 0, voxel_grid.shape[0]-1)
    x1, y1, z1 = np.clip([x1, y1, z1], 0, voxel_grid.shape[0]-1)

    # Compute weights
    xd, yd, zd = [arr - arr.astype(int) for arr in [x, y, z]]

    # Interpolate along x
    c00 = voxel_grid[x0, y0, z0]*(1-xd) + voxel_grid[x1, y0, z0]*xd
    c01 = voxel_grid[x0, y0, z1]*(1-xd) + voxel_grid[x1, y0, z1]*xd
    c10 = voxel_grid[x0, y1, z0]*(1-xd) + voxel_grid[x1, y1, z0]*xd
    c11 = voxel_grid[x0, y1, z1]*(1-xd) + voxel_grid[x1, y1, z1]*xd

    # Interpolate along y
    c0 = c00*(1-yd) + c10*yd
    c1 = c01*(1-yd) + c11*yd

    # Interpolate along z
    c = c0*(1-zd) + c1*zd

    return c

def rowwise_corrcoef(A, B, mask=None):
    assert A.shape == B.shape, f"A and B must have the same shape, got: {A.shape} and {B.shape}"
    
    if mask is not None:
        assert mask.shape == A.shape, "mask must have the same shape as A and B"
        A = ma.masked_array(A, mask=np.logical_not(mask))
        B = ma.masked_array(B, mask=np.logical_not(mask))

    # Calculate means
    A_mean = ma.mean(A, axis=1, keepdims=True)
    B_mean = ma.mean(B, axis=1, keepdims=True)
    
    # Subtract means
    A_centered = A - A_mean
    B_centered = B - B_mean
    
    # Calculate sum of products
    sumprod = ma.sum(A_centered * B_centered, axis=1)
    
    # Calculate square roots of the sum of squares
    sqrt_sos_A = ma.sqrt(ma.sum(A_centered**2, axis=1))
    sqrt_sos_B = ma.sqrt(ma.sum(B_centered**2, axis=1))
    
    # Return correlation coefficients
    cc =  sumprod / (sqrt_sos_A * sqrt_sos_B)
    return cc.data


def sphere_points_np(ctr,rad,N):
    # ctr is shape (N,3)
    #print("sphere_points_ctr_input_shape",ctr.shape)
    h = -1.0 + (2.0 * np.arange(N) / float(N-1))[:, np.newaxis]
    phis = np.arccos(h)
    thetas = np.zeros_like(phis)
    a = (3.6 / np.sqrt(N * (1.0 - h[1:-1]**2)))
    thetas[1:-1, :] = a
    thetas = np.cumsum(thetas,axis=0)
    x = np.sin(phis) * np.cos(thetas)
    y = np.sin(phis) * np.sin(thetas)
    z = np.cos(phis)
    points = rad * np.stack([x, y, z], axis=-1)
    points = points.reshape(1, N, 3)
    ctr = ctr[:, np.newaxis, :]
    points = ctr + points # broadcast add
    return points
    
    
def sphere_points_flex(ctr,rad,N):
    h = -1.0 + (2.0 * flex.double_range(N) / (N-1))
    phis = flex.acos(h)
    thetas = flex.double(len(phis),0.0)
    a = (3.6 / flex.sqrt(N * (1.0 - h[1:-1]**2)))
    thetas = thetas.set_selected(flex.size_t_range(1,N-1),a)
    
    # cumulative sum operation
    def cumsum_flex(arr):
        # should perform as np.cumsum
        result = []
        running_sum = 0.0
        for i,x in enumerate(arr):
            running_sum += x
            result.append(running_sum)
        return flex.double(result)
    
    
    thetas = cumsum_flex(thetas)
    x = flex.sin(phis) * flex.cos(thetas)
    y = flex.sin(phis) * flex.sin(thetas)
    z = flex.cos(phis)
    rad = float(rad)
    points = rad * flex.vec3_double(x,y,z)
    
    # add generated points to center points
    def broadcast_add(ctr, points):

        # add and populate additional dimension. equivalent to:
        # result = ctr[:, np.newaxis, :] + points
        
        N = points.size()
        M = ctr.size()
        # Preallocate an array of shape (M*N, 3)
        result = flex.vec3_double(M*N)

        for i in range(M):
            for j in range(N):
                flat_index = i * N + j
                new_point = tuple(ctr[i:i+1] + points[j:j+1])[0]
                result[flat_index] = new_point


        result = result.as_1d().as_double()
        result.reshape(flex.grid(len(ctr),len(points),3))
        return result
    
    # apply function
    points = broadcast_add(ctr,points)
    return points


def cdist_flex(A,B):
    

    def indices_2d_flex(dimensions):
        N = len(dimensions)
        if N != 2:
            raise ValueError("Only 2D is supported for this implementation.")

        # Create the row indices
        row_idx = flex.size_t(chain.from_iterable([[i] * dimensions[1] for i in range(dimensions[0])]))

        # Create the column indices
        col_idx = flex.size_t(chain.from_iterable([list(range(dimensions[1])) for _ in range(dimensions[0])]))

        return row_idx, col_idx

        
    i_idxs, j_idxs = indices_2d_flex((A.focus()[0],B.focus()[0]))
    
    r = i_idxs
    xi = i_idxs*3
    yi = i_idxs*3 + 1
    zi = i_idxs*3 + 2

    xa = A.select(xi)
    ya = A.select(yi)
    za = A.select(zi)

    xj = j_idxs*3
    yj = j_idxs*3 + 1
    zj = j_idxs*3 + 2


    xb = B.select(xj)
    yb = B.select(yj)
    zb = B.select(zj)

    d = ((xb - xa)**2 + (yb - ya)**2 + (zb - za)**2)**0.5
    d.reshape(flex.grid((A.focus()[0],B.focus()[0])))
    
    return d

def query_atom_neighbors(model,radius=3.5,include_self=True,only_unit=True):
    crystal_symmetry = model.crystal_symmetry()
    hierarchy = model.get_hierarchy()
    sites_cart = hierarchy.atoms().extract_xyz()
    sst = crystal_symmetry.special_position_settings().site_symmetry_table(
    sites_cart = sites_cart)
    conn_asu_mappings = crystal_symmetry.special_position_settings().\
    asu_mappings(buffer_thickness=5)
    conn_asu_mappings.process_sites_cart(
    original_sites      = sites_cart,
    site_symmetry_table = sst)
    conn_pair_asu_table = cctbx.crystal.pair_asu_table(
    asu_mappings=conn_asu_mappings)
    conn_pair_asu_table.add_all_pairs(distance_cutoff=radius)
    pair_generator = cctbx.crystal.neighbors_fast_pair_generator(
    conn_asu_mappings,
    distance_cutoff=radius)
    fm = crystal_symmetry.unit_cell().fractionalization_matrix()
    om = crystal_symmetry.unit_cell().orthogonalization_matrix()


    pairs = list(pair_generator)
    inds = defaultdict(list)
    dists = defaultdict(list)
    
    for pair in pairs:
        i,j = pair.i_seq, pair.j_seq
        rt_mx_i = conn_asu_mappings.get_rt_mx_i(pair)
        rt_mx_j = conn_asu_mappings.get_rt_mx_j(pair)
        rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)

                
        if (only_unit and rt_mx_ji.is_unit_mx()) or (not only_unit):
            d = round(math.sqrt(pair.dist_sq),6)
            inds[i].append(j)
            dists[i].append(d)
            
            # add reverse
            inds[j].append(i)
            dists[j].append(d)
            #print(pair.i_seq, pair.j_seq, rt_mx_ji, math.sqrt(pair.dist_sq), de)
    
    # add self
    if include_self:
        for key,value in list(inds.items()):
            dval = dists[key]
            dists[key]= dval+[0.0]
            inds[key] = value+[key]

    # sort
    for key,value in list(inds.items()):
        dval = dists[key]
        # sort
        sorted_pairs = sorted(set(list(zip(value,dval))))
        value_sorted, dval_sorted = zip(*sorted_pairs)
        inds[key] = value_sorted
        dists[key] = dval_sorted


    return inds,dists



def query_ball_point_flex(tree,tree_xyz,query_xyz,r=None):
    assert r is not None, "provide radius"
    n_atoms,n_probes, _ = query_xyz.focus()
    counts = []

    for atom_i in range(n_atoms):
        probe_range = (n_probes * atom_i * 3, n_probes * (atom_i+1) * 3)
        atom_probes_xyz = query_xyz.select(flex.size_t_range(*probe_range))
        atom_probes_xyz.reshape(flex.grid(n_probes,3))
        nbrs = tree[atom_i]
        n_nbrs = len(nbrs)
        nbrs_xyz = tree_xyz.select(flex.size_t(nbrs)).as_1d().as_double()
        nbrs_xyz.reshape(flex.grid(len(nbrs),3))
        d = cdist_flex(nbrs_xyz,atom_probes_xyz)
        sel = d<r
        count = []
        for nbr_i in range(n_probes):
            nbr_range = (slice(0,n_nbrs),slice(nbr_i,nbr_i+1))
            count_nbr = sel[nbr_range].count(True)
            count.append(count_nbr)

        counts.append(count)

    counts = flex_from_list(counts)
    return counts

