#ifndef SCITBX_GRAPHICS_UTILS_COLOR_H
#define SCITBX_GRAPHICS_UTILS_COLOR_H

#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/error.h>

#include <cmath>

namespace scitbx { namespace graphics_utils {

  namespace af = scitbx::af;

  // Hue, Saturation, Value --> Red, Green, Blue
  inline
  scitbx::vec3<double>
  hsv2rgb (double h, double s, double v)
  {
    if (s == 0) {
      return scitbx::vec3<double>(v, v, v);
    }
    h /= 60.0;
    v *= 255.0;
    int i = std::floor(h);
    double f = h - i;
    double p = v * (1 - s);
    double q = v * (1 - (s * f));
    double t = v * (1 - (s * (1 - f)));
    switch (i) {
      case 0  : return scitbx::vec3<double>(v, t, p) / 255.0; break;
      case 1  : return scitbx::vec3<double>(q, v, p) / 255.0; break;
      case 2  : return scitbx::vec3<double>(p, v, t) / 255.0; break;
      case 3  : return scitbx::vec3<double>(p, q, v) / 255.0; break;
      case 4  : return scitbx::vec3<double>(t, p, v) / 255.0; break;
      default : break;
    }
    return scitbx::vec3<double>(v, p, q) / 255.0;
  }

  // this function may be superfluous here, but could be useful elsewhere
  af::shared< scitbx::vec3<double> >
  make_rainbow_gradient (unsigned nbins)
  {
    SCITBX_ASSERT(nbins > 1);
    af::shared< scitbx::vec3<double> > color_gradient(nbins);
    double f_nbins(nbins);
    for (unsigned i = 0; i < nbins; i++) {
      double gradient_ratio = i / (f_nbins-1);
      color_gradient[i] = hsv2rgb(240.0 - (240 * gradient_ratio), 1., 1.);
    }
    return color_gradient;
  }

  // gradient generation is inlined for speed.
  af::shared< scitbx::vec3<double> >
  color_rainbow (
    af::const_ref< bool > const& selection,
    bool color_all=false)
  {
    size_t n_selected = 0;
    for (unsigned i_seq = 0; i_seq < selection.size(); i_seq++) {
      if (selection[i_seq]) n_selected++;
    }
    SCITBX_ASSERT(n_selected > 1);
    af::shared< scitbx::vec3<double> > colors(selection.size());
    unsigned j = 0;
    double f_bins = (double) n_selected;
    if (color_all) f_bins = (double) selection.size();
    for (unsigned i_seq = 0; i_seq < selection.size(); i_seq++) {
      if ((selection[i_seq]) || (color_all)) {
        double gradient_ratio = j / (f_bins-1);
        colors[i_seq] = hsv2rgb(240.0 - (240 * gradient_ratio), 1., 1.);
        j++;
      } else {
        colors[i_seq] = scitbx::vec3<double>(0.0,0.0,0.0);
      }
    }
    return colors;
  }

  af::shared< scitbx::vec3<double> >
  scale_selected_colors (
    af::const_ref< scitbx::vec3<double> > const& input_colors,
    af::const_ref< bool > const& selection,
    double scale=0.5)
  {
    SCITBX_ASSERT(input_colors.size() == selection.size());
    SCITBX_ASSERT(scale >= 0);
    af::shared< scitbx::vec3<double> > atom_colors(input_colors.size());
    for (unsigned i_seq = 0; i_seq < input_colors.size(); i_seq++) {
      scitbx::vec3<double> c = input_colors[i_seq];
      if (selection[i_seq]) {
        c[0] *= scale;
        c[1] *= scale;
        c[2] *= scale;
      }
      atom_colors[i_seq] = c;
    }
    return atom_colors;
  }

  af::shared< scitbx::vec3<double> >
  color_by_property (
    af::const_ref< double > const& properties,
    af::const_ref< bool > const& selection,
    bool color_all=false,
    unsigned gradient_type=0,
    double min_value=0.1)
  {
    SCITBX_ASSERT(properties.size() > 0);
    SCITBX_ASSERT(gradient_type <= 2);
    af::shared <scitbx::vec3<double> > colors(properties.size());
    double vmax = properties[0];
    double vmin = properties[0];
    for (unsigned i_seq = 0; i_seq < properties.size(); i_seq++) {
      if ((! color_all) && (! selection[i_seq])) continue;
      if (properties[i_seq] > vmax) vmax = properties[i_seq];
      if (properties[i_seq] < vmin) vmin = properties[i_seq];
    }
    for (unsigned i_seq = 0; i_seq < properties.size(); i_seq++) {
      double gradient_ratio = (properties[i_seq]-vmin) / (vmax-vmin);
      if ((! color_all) && (! selection[i_seq])) { // black
        colors[i_seq] = scitbx::vec3<double>(0.0,0.0,0.0);
      } else if (gradient_type == 1) { // red-blue
        colors[i_seq] = hsv2rgb(240.0 + (120 * gradient_ratio), 1., 1.);
      } else if (gradient_type == 2) { // heatmap
        double h = 0.;
        double s = 1.;
        double v = 1.;
        if (gradient_ratio < 0.35) {
          double ratio_norm = gradient_ratio / 0.35;
          v = min_value + (1. - min_value) * (ratio_norm * ratio_norm);
          s = gradient_ratio / 0.35;
        } else if (gradient_ratio < 0.75) {
          h = 60. - (60 * (0.75 - gradient_ratio) / 0.4);
        } else {
          h = 60.;
          s = 1. - (gradient_ratio - 0.75) / 0.25;
        }
        colors[i_seq] = hsv2rgb(h, s, v);
      } else { // rainbow
        colors[i_seq] = hsv2rgb(240.0 - (240 * gradient_ratio), 1., 1.);
      }
    }
    return colors;
  }

  af::shared< scitbx::vec3<double> >
  grayscale_by_property (
    af::const_ref< double > const& properties,
    af::const_ref< bool > const& selection,
    bool shade_all=false,
    bool invert=false,
    double max_value=0.95,
    double max_value_inverted=0.1)
  {
    SCITBX_ASSERT(properties.size() > 0);
    af::shared <scitbx::vec3<double> > colors(properties.size());
    double vmax = properties[0];
    double vmin = properties[0];
    for (unsigned i_seq = 0; i_seq < properties.size(); i_seq++) {
      if ((! shade_all) && (! selection[i_seq])) continue;
      if (properties[i_seq] > vmax) vmax = properties[i_seq];
      if (properties[i_seq] < vmin) vmin = properties[i_seq];
    }
    if (vmax == vmin) {
      vmax = 1.0;
      vmin = 0.0;
    }
    for (unsigned i_seq = 0; i_seq < properties.size(); i_seq++) {
      double gradient_ratio = (properties[i_seq]-vmin) / (vmax-vmin);
      if ((! shade_all) && (! selection[i_seq])) {
        if (invert) {
          colors[i_seq] = scitbx::vec3<double>(0.0,0.0,0.0);
        } else {
          colors[i_seq] = scitbx::vec3<double>(1.0,1.0,1.0);
        }
      } else if (invert) {
        double value = max_value_inverted + (gradient_ratio *
          (1.0-max_value_inverted));
        colors[i_seq] = scitbx::vec3<double>(value, value, value);
      } else {
        double value = max_value - (max_value * gradient_ratio);
        colors[i_seq] = scitbx::vec3<double>(value, value, value);
      }
    }
    return colors;
  }

}} // namespace scitbx::graphics_utils

#endif
