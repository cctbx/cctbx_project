#include "nvToolsExt.h"

#include <simtbx/nanoBragg/nanoBragg.h>



namespace simtbx {
namespace nanoBragg {

void nanoBragg::add_nanoBragg_spots_nvtx(){

    nvtxRangePushA("add_nanoBragg_spots");
    add_nanoBragg_spots();
    nvtxRangePop();

}


void nanoBragg::add_nanoBragg_spots_cuda_nvtx(){

    nvtxRangePushA("add_nanoBragg_spots_cuda");
    add_nanoBragg_spots_cuda();
    nvtxRangePop();

}

} // namespace nanoBragg
} // namespace simtbx
