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



void nanoBragg::add_noise_nvtx(){

    nvtxRangePushA("add_noise");
    add_noise();
    nvtxRangePop();

}



void nanoBragg::to_smv_format_nvtx(std::string const & fileout, double
                                   intfile_scale, int debug_x, int debug_y){

    nvtxRangePushA("to_smv_format");
    to_smv_format(fileout, intfile_scale, debug_x, debug_y);
    nvtxRangePop();

}



} // namespace nanoBragg
} // namespace simtbx
