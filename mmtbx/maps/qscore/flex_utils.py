import math
from cctbx.array_family import flex

def flex_from_list(lst,signed_int=False):
    flat_list, shape = flatten_and_shape(lst)
    dtype = get_dtype_of_list(flat_list)
    type_mapper = {int:flex.size_t,
                   float:flex.double,
                   bool:flex.bool}
    if signed_int:
        type_mapper[int] = flex.int16
    
    # make flex array
    assert dtype in type_mapper, f"Unrecognized type: {dtype}"
    flex_func = type_mapper[dtype]
    flex_array = flex_func(flat_list)
    if len(shape)>1:
        flex_array.reshape(flex.grid(*shape))
    return flex_array



def flatten_and_shape(lst):
    """Flatten a nested list and return its shape."""
    def helper(l):
        if not isinstance(l, list):
            return [l], ()
        flat = []
        shapes = []
        for item in l:
            f, s = helper(item)
            flat.extend(f)
            shapes.append(s)
        if len(set(shapes)) != 1:
            raise ValueError("Ragged nested list detected.")
        return flat, (len(l),) + shapes[0]

    flattened, shape = helper(lst)
    return flattened, shape


def get_dtype_of_list(lst):
    dtypes = {type(item) for item in lst}
    
    if len(dtypes) > 1:
        raise ValueError("Multiple data types detected.")
    elif len(dtypes) == 0:
        raise ValueError("Empty list provided.")
    else:
        return dtypes.pop()
    
def optimized_nd_to_1d_indices(i, shape):
    # For fixed input of (None, i, None), we directly compute based on given structure
    result_indices = []
    
    # Pre-compute for 1st dimension which is always a slice
    start1, stop1 = 0, shape[0]
    
    # Pre-compute for 3rd dimension which is always a slice
    start3, stop3 = 0, shape[2]
    stride3 = 1
    
    # Directly compute for 2nd dimension which is variable
    stride2 = shape[2]
    index2 = i * stride2 * shape[0]
    
    for val1 in range(start1, stop1):
        for val3 in range(start3, stop3):
            result_indices.append(val1 * stride2 + index2 + val3 * stride3)
            
    return result_indices

def nd_to_1d_indices(indices, shape):
    # Normalize indices to always use slice objects
    normalized_indices = []
    for dim, idx in enumerate(indices):
        if idx is None:
            normalized_indices.append(slice(0,shape[dim]))
        else:
            normalized_indices.append(idx)
    
    # If any index is a slice, recursively call function for each value in slice
    for dim, (i, s) in enumerate(zip(normalized_indices, shape)):
        if isinstance(i, slice):
            result_indices = []
            start, stop, step = i.indices(s)
            for j in range(start, stop, step):
                new_indices = list(normalized_indices)
                new_indices[dim] = j
                result_indices.extend(nd_to_1d_indices(new_indices, shape))
            return result_indices
    
    # If no slices, calculate single 1D index
    index = 0
    stride = 1
    for i, dim in reversed(list(zip(normalized_indices, shape))):
        index += i * stride
        stride *= dim
    return [index]


def flex_std(flex_array):
    n = flex_array.size()
    if n <= 1:
        raise ValueError("Sample size must be greater than 1")
    
    # Compute the mean
    mean_value = flex.mean(flex_array)
    
    # Compute the sum of squared deviations
    squared_deviations = (flex_array - mean_value) ** 2
    sum_squared_deviations = flex.sum(squared_deviations)
    
    # Compute the standard deviation
    std_dev = (sum_squared_deviations / (n - 1)) ** 0.5
    return std_dev