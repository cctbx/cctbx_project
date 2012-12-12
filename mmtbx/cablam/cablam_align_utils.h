// (jEdit options) :folding=indent:collapseFolds=1:
#include <boost/python.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <assert.h>
#include <cmath>


namespace mmtbx { namespace cablam {

  struct index_mean
  {
    // This is just a little object to hold inecies for simialr regions
    // found by get_similar_regions, which returns a list of these objects
    index_mean(
      int i_1,
      int i_2,
      double mean,
      int window_length)
    :
      i_1(i_1),
      i_2(i_2),
      mean(mean),
      window_length(window_length)
    {}
    int i_1;
    int i_2;
    double mean;
    int window_length;
  };

  boost::python::list get_similar_regions(
    boost::python::list list_1,
    boost::python::list list_2,
    double threshold,
    int window_len)
  {
    // This function returns a a list of index_mean objects,
    // if im = index_mean object then im has three attributes, im.i_1 is the
    // index of list 1, im.i_2 is the index of list 2, and im.mean is the mean
    // of the diff_list as explained below.
    // Only diff_lists with their mean <= threshold will be returned in the
    // return list. IF threshold !> 0 then returns the fragments with the
    // lowest diff_list mean
    boost::python::ssize_t i_1, i_2, i_w;
    boost::python::list return_list;
    double smallest_mean = 9999;
    //
    for(i_1=0;i_1<boost::python::len(list_1) - window_len + 1;i_1++) {
      // get window
      boost::python::list window;// ,w;
      bool none_in = false;
      for(i_w=i_1;i_w<window_len + i_1;i_w++) {
        boost::python::extract<double> value(list_1[i_w]);
        if( !value.check() ) none_in = true;
        window.append(list_1[i_w]);
      }
      if(none_in) continue;
      // slide the window across list_2 and check for similarities
      for(i_2=0;i_2<boost::python::len(list_2) - window_len + 1;i_2++) {
        // get target
        boost::python::list target;
        bool none_in = false;
        for(i_w=i_2;i_w<window_len + i_2;i_w++) {
          boost::python::extract<double> value(list_2[i_w]);
          if( !value.check() ) none_in = true;
          target.append(list_2[i_w]);
        }
        if(none_in) continue;
        assert(boost::python::len(target) == boost::python::len(window));
        // get diff list
        boost::python::list diff_list;
        for(i_w=0;i_w<window_len;i_w++) {
          double v1 = boost::python::extract<double>(window[i_w]);
          double v2 = boost::python::extract<double>(target[i_w]);
          double diff = std::abs(v1 - v2);
          if(diff > 180) diff = 360 - diff;
          diff_list.append(diff);
        }
        assert(boost::python::len(window) == boost::python::len(diff_list));
        // get the average of the diff_list
        double sum = 0;
        for(int i_d=0;i_d<boost::python::len(diff_list);i_d++)
          sum = sum + boost::python::extract<double>(diff_list[i_d]);
        double average;
        average = sum/boost::python::len(diff_list);
        if((threshold > 0) && (average <= threshold))
          return_list.append(index_mean(i_1, i_2, average, window_len));
        else if(threshold <= 0) {
          // return_list.append(index_mean(i_1, i_2, average));
          // return_list.append(average);
          if(average < smallest_mean) {
            smallest_mean = average;
            if(boost::python::len(return_list) > 0) return_list.pop(0);
            return_list.append(index_mean(i_1, i_2, average, window_len));
          }
        }
      }
    }
    return return_list;
  }


}}
