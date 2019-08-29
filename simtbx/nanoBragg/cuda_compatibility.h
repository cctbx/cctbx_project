
#ifndef SIMTBX_CUDA_COMPAT_H
#define SIMTBX_CUDA_COMPAT_H

// Runtime error encountered, __ldg doesn't exist for compute capability < 3.5
// if you tell nvcc your cimpiler is > 3.5 , it will compile just fine
// but then the code will fail during runtime if your GPU is too old..
// https://stackoverflow.com/a/27302007/2077270
template<typename T>
__device__ __forceinline__ T __ldg(const T* ptr) {
#if __CUDA_ARCH__ >= 350
    return __ldg(ptr);
#else
    return *ptr;
#endif
}

#endif // SIMTBX_CUDA_COMPAT_H

