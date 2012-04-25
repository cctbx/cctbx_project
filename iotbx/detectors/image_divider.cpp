#include <iotbx/detectors/image_divider.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>

Distl::image_divider::image_divider(scitbx::af::flex_int d, const int& n):
  data(d),nullvalue(n),slow_tiles(interval_list()),
                       fast_tiles(interval_list()){

  int slow_size = data.accessor().focus()[0];
  int fast_size = data.accessor().focus()[1];
  scitbx::af::shared<bool> slow_flag(slow_size,true);
  scitbx::af::shared<bool> fast_flag(fast_size,true);
  int* icount = data.begin();
  for (int islow = 0; islow<slow_size; ++islow){
    for (int ifast = 0; ifast<fast_size; ++ifast){
      if (*icount++!=nullvalue){
        slow_flag[islow]=false;
        fast_flag[ifast]=false;
      }
    }
  }

  //analysis of slow modules
  bool bswitch = true;
  int begin,end,islow,ifast;
  for (islow= 0; islow < slow_size; ++islow){
    if (bswitch && !slow_flag[islow]){
      begin = islow;
      bswitch = !bswitch;
    }
    if (!bswitch && slow_flag[islow]){
      end = islow-1;
      bswitch = !bswitch;
      slow_tiles.push_back(interval(begin,end));
    }
  }
  if (!bswitch) {
    end = islow - 1;
    slow_tiles.push_back(interval(begin,end));
  }

  //analysis of fast modules
  bswitch = true;
  for (ifast= 0; ifast < fast_size; ++ifast){
    if (bswitch && !fast_flag[ifast]){
      begin = ifast;
      bswitch = !bswitch;
    }
    if (!bswitch && fast_flag[ifast]){
      end = ifast-1;
      bswitch = !bswitch;
      fast_tiles.push_back(interval(begin,end));
    }
  }
  if (!bswitch) {
    end = ifast - 1;
    fast_tiles.push_back(interval(begin,end));
  }
}

int
Distl::image_divider::module_count()const{
  return slow_tiles.size() * fast_tiles.size();
}

scitbx::af::flex_int
Distl::image_divider::tile_data(const int& module_number)const{
  int n_slow = slow_tiles.size();
  int n_fast = fast_tiles.size();

  int idx_slow = module_number / n_fast;
  int idx_fast = module_number % n_slow;

  int target_size_slow = slow_tiles[idx_slow].size();
  int target_size_fast = fast_tiles[idx_fast].size();

  scitbx::af::flex_int z(scitbx::af::flex_grid<>(target_size_slow,target_size_fast));
  int* begin = z.begin();
  const int* icount = data.begin();

  int slow_size = data.accessor().focus()[0];
  int fast_size = data.accessor().focus()[1];

  icount += slow_tiles[idx_slow].first * fast_size
         + fast_tiles[idx_fast].first;
  for (int islow = 0; islow<target_size_slow; ++islow){
    for (int ifast = 0; ifast<target_size_fast; ++ifast){
      *begin++ = *icount++;
    }
    icount += fast_size - target_size_fast;
  }

  return z;
}

Distl::interval
Distl::image_divider::tile_slow_interval(const int& module_number)const{
  int n_fast = fast_tiles.size();

  int idx_slow = module_number / n_fast;
  SCITBX_ASSERT(idx_slow < slow_tiles.size());

  return slow_tiles[idx_slow];
}

Distl::interval
Distl::image_divider::tile_fast_interval(const int& module_number)const{
  int n_fast = fast_tiles.size();

  int idx_fast = module_number % n_fast;
  SCITBX_ASSERT(idx_fast < fast_tiles.size());

  return fast_tiles[idx_fast];
}
