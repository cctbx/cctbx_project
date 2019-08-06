#ifndef CUDA_STRUCT_H
#define CUDA_STRUCT_H

/*
  Header for defining a struct to house pointers on GPU
*/

#include <simtbx/nanoBragg/nanotypes.h>

#ifndef CUDAREAL
#define CUDAREAL double
#endif

#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

struct hklParams {
        int hkls;
        int h_min;
        int h_max;
        int h_range;
        int k_min;
        int k_max;
        int k_range;
        int l_min;
        int l_max;
        int l_range;
};

#endif
