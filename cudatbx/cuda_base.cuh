#ifndef CUDA_BASE_H
#define CUDA_BASE_H

#include <iostream>

// CUDA includes
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <math_constants.h>
#include <math_functions.h>

#define cudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
  if( cudaSuccess != err) {
    std::cerr << "File:  " << file << "\n"
              << "Line:  " << line << "\n"
              << "Error: " << cudaGetErrorString(err) << "\n";
    exit(-1);
  }
}

#endif // CUDA_BASE_H
