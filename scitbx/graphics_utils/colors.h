#ifndef SCITBX_GRAPHICS_UTILS_COLOR_H
#define SCITBX_GRAPHICS_UTILS_COLOR_H

//#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <cctbx/hendrickson_lattman.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/error.h>
#include <scitbx/constants.h>

#include <cmath>
#include <cstdio>
#include <boost/math/special_functions/fpclassify.hpp> // provides isfinite()

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

  inline
  scitbx::vec3<double>
  get_heatmap_color (
    double gradient_ratio,
    double min_value=0.1)
  {
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
    return hsv2rgb(h, s, v);
  }

  inline
  scitbx::vec3<double>
    get_Phi_FOM_colour(double phi, double fom = 1.0)
  {
    /* return circular rainbow colour indicating phase between [0;2Pi] but greyed according to fom between [0;1]
    inspired by https://en.wikipedia.org/wiki/HSL_and_HSV and python code:
     def HSV_vivid(h):
        h %= 1
        h *= 6
        c = 1
        x = 1 - abs((h % 2) - 1)
        rgb = ()
        if   h < 1: rgb = (c, x, 0)
        elif h < 2: rgb = (x, c, 0)
        elif h < 3: rgb = (0, c, x)
        elif h < 4: rgb = (0, x, c)
        elif h < 5: rgb = (x, 0, c)
        else:       rgb = (c, 0, x)
        return rgb[0], rgb[1], rgb[2]
  */
    if (boost::math::isfinite(phi) == false || boost::math::isfinite(fom) == false)
      return scitbx::vec3<double>(0.5, 0.5, 0.5);

    double h = fmod(phi, 2.0*scitbx::constants::pi)/(2.0*scitbx::constants::pi);

    while (h < 0.0)
      h++;
    h *= 6;
    double mhfom = 0.5 - 0.5*fom;
    double c = (1.0 + fom)/2.0;
    double x = (1.0 - fabs(fmod(h, 2.0) - 1.0)) * fom + mhfom;
    double r, g, b;
    if (h < 1.0) {
      r = c; g = x; b = mhfom;
    }
    else if (h < 2.0) {
      r = x; g = c; b = mhfom;
    }
    else if (h < 3.0) {
      r = mhfom; g = c; b = x;
    }
    else if (h < 4.0) {
      r = mhfom; g = x; b = c;
    }
    else if (h < 5.0) {
      r = x; g = mhfom; b = c;
    }
    else {
      r = c; g = mhfom; b = x;
    }
    return scitbx::vec3<double>(r, g, b);
  }


  af::shared< scitbx::vec3<double> >
    NoNansvec3(af::const_ref< scitbx::vec3<double> > const& vecs,
      double defx = 0.0, double defy = 0.0, double defz = 0.0)
  {
    af::shared <scitbx::vec3<double> > nonanvecs(vecs.size());
    for (unsigned i_seq = 0; i_seq < vecs.size(); i_seq++)
    {
      if (boost::math::isfinite(vecs(i_seq)[0] + vecs(i_seq)[1] + vecs(i_seq)[2]))
        nonanvecs[i_seq] = vecs(i_seq);
      else
        nonanvecs[i_seq] = scitbx::vec3<double>(defx, defy, defz);
    }
    return nonanvecs;
  }


  af::shared< cctbx::hendrickson_lattman<> >
    NoNansHL(af::const_ref< cctbx::hendrickson_lattman<> > const& HL,
      double A = 0.0, double B = 0.0, double C = 0.0, double D = 0.0)
  {
    af::shared <cctbx::hendrickson_lattman<> > nonanshl(HL.size());
    for (unsigned i = 0; i < HL.size(); i++)
    {
      if (boost::math::isfinite(HL[i].a() + HL[i].b() + HL[i].c() + HL[i].d()))
        nonanshl[i] = cctbx::hendrickson_lattman<>(HL[i].a(), HL[i].b(), HL[i].c(), HL[i].d());
      else
        nonanshl[i] = cctbx::hendrickson_lattman<>(A, B, C, D);
    }
    return nonanshl;
  }


  af::shared< double>
  NoNans(af::const_ref< double > const& arr, double def = 0.0)
  {
    af::shared <double> nonanarr(arr.size());
    for (unsigned i_seq = 0; i_seq < arr.size(); i_seq++)
    {
      if (boost::math::isfinite(arr(i_seq)))
        nonanarr[i_seq] = arr(i_seq);
      else
        nonanarr[i_seq] = def;
    }
    return nonanarr;
  }


  af::shared< bool>
  IsNans(af::const_ref< double > const& arr)
  {
    af::shared <bool> nonanarr(arr.size());
    for (unsigned i_seq = 0; i_seq < arr.size(); i_seq++)
    {
      if (boost::math::isfinite(arr(i_seq)))
        nonanarr[i_seq] = false;
      else
        nonanarr[i_seq] = true;
    }
    return nonanarr;
  }


  af::shared< bool >
  IsNansvec3(af::const_ref< scitbx::vec3<double> > const& vecs)
  {
    af::shared <bool > nonanarr(vecs.size());
    for (unsigned i_seq = 0; i_seq < vecs.size(); i_seq++)
    {
      if (boost::math::isfinite(vecs(i_seq)[0] + vecs(i_seq)[1] + vecs(i_seq)[2]))
        nonanarr[i_seq] = false;
      else
        nonanarr[i_seq] = true;
    }
    return nonanarr;
  }


  double round2(double const &val, int const& precision)
  {
    int d = 0;
    if ((val * pow(10.0, precision + 1)) - (floor(val * pow(10.0, precision))) > 4)
      d = 1;
    return (floor(val * pow(10.0, precision)) + d) / pow(10.0, precision);
  }


  double flt_roundoff(double const& val, int const& precision)
  {
    // fast version of libtbx.math_utils.roundoff() for a single double value only
    if (!boost::math::isfinite(val))
      return FP_NAN;

    if (fabs(val) < pow(10.0, -precision))
    {
      const int ndigits = 50;
      char fstr[ndigits];
      strcpy(fstr, "%");
      char fstr1[ndigits];
      snprintf(fstr1, sizeof(fstr1), "%d.", precision);
      strcat(fstr, fstr1);
      snprintf(fstr1, sizeof(fstr1), "%d", precision);
      strcat(fstr, fstr1);
      strcat(fstr, "e");

      char fstr2[ndigits];
      snprintf(fstr2, sizeof(fstr2), fstr, val);
      return atof(fstr2);
    }
    return round2(val, precision);
  }


  af::shared< scitbx::vec3<double> >
  color_by_phi_fom(
    af::const_ref< double > const& phases,
    af::const_ref< double > const& foms
  )
  {
    SCITBX_ASSERT(phases.size() == foms.size());
    af::shared <scitbx::vec3<double> > colors(phases.size());

    for (unsigned i_seq = 0; i_seq < phases.size(); i_seq++)
      colors[i_seq] = get_Phi_FOM_colour(phases[i_seq], foms[i_seq]);

    return colors;
  }


  // this function may be superfluous here, but could be useful elsewhere
  af::shared< scitbx::vec3<double> >
  make_rainbow_gradient (unsigned nbins)
  {
    SCITBX_ASSERT(nbins > 0);
    af::shared< scitbx::vec3<double> > color_gradient(nbins);
    double f_nbins(nbins);
    for (unsigned i = 0; i < nbins; i++) {
      double gradient_ratio = 0;
      if (nbins > 1) {
        gradient_ratio = i / (f_nbins-1);
      }
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
        double gradient_ratio = 0;
        if (n_selected > 0) {
          gradient_ratio = j / (f_bins-1);
        }
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
    color_by_property(
      af::const_ref< double > const& properties,
      af::const_ref< bool > const& selection,
      bool color_all = false,
      unsigned gradient_type = 0,
      double min_value = 0.1)
  {
    SCITBX_ASSERT(properties.size() > 0);
    SCITBX_ASSERT(gradient_type <= 2);
    af::shared <scitbx::vec3<double> > colors(properties.size());
    double vmax = -9e99;
    double vmin = 9e99;
    for (unsigned i_seq = 0; i_seq < properties.size(); i_seq++) {
      if ((!color_all) && (!selection[i_seq])) continue;
      if (!boost::math::isfinite(properties[i_seq])) continue;
      if (properties[i_seq] > vmax) vmax = properties[i_seq];
      if (properties[i_seq] < vmin) vmin = properties[i_seq];
    }
    if (vmax == vmin) {
      vmax = 1.0;
      vmin = 0.0;
    }
    for (unsigned i_seq = 0; i_seq < properties.size(); i_seq++) {
      double gradient_ratio = (properties[i_seq] - vmin) / (vmax - vmin);
      if ((!color_all) && (!selection[i_seq])) { // black
        colors[i_seq] = scitbx::vec3<double>(0.0, 0.0, 0.0);
      }
      else if (gradient_type == 0) { // rainbow
        colors[i_seq] = hsv2rgb(240.0 - (240 * gradient_ratio), 1., 1.);
      }
      else if (gradient_type == 1) { // red-blue
        colors[i_seq] = hsv2rgb(240.0 + (120 * gradient_ratio), 1., 1.);
      }
      else if (gradient_type == 2) { // heatmap
        colors[i_seq] = get_heatmap_color(gradient_ratio, min_value);
      }
    }
    return colors;
  }

  af::shared< scitbx::vec3<double> >
    map_to_rgb_colourmap(
      af::const_ref< double > const& data_for_colors,
      af::shared< scitbx::vec3<double> > const& colourmap,
      af::const_ref< bool > const& selection,
      af::const_ref< double > const& attenuation,
      double powscale = 1.0,
      bool map_directly = false,
      bool color_all = false
      )
  {
    /*
    Map data_for_colors to the colours given in colourmap
    data_for_colors could be a large array, say ~ number of reflections in a data set.
    colourmap is a smallish array of say 200 different colours, like a smooth rgb gradient
    obtained from https://matplotlib.org/examples/color/colormaps_reference.html
    */
    SCITBX_ASSERT(data_for_colors.size() > 0 && data_for_colors.size() == attenuation.size());
    af::shared <scitbx::vec3<double> > colors(data_for_colors.size());
    double vmax = -9e99;
    double vmin = 9e99;
    int nrgbcolours = colourmap.size();
    for (unsigned i_seq = 0; i_seq < data_for_colors.size(); i_seq++)
    {
      if ((!color_all) && (!selection[i_seq])) continue;
      if (!boost::math::isfinite(data_for_colors[i_seq])) continue;
      if (data_for_colors[i_seq] > vmax) vmax = data_for_colors[i_seq];
      if (data_for_colors[i_seq] < vmin) vmin = data_for_colors[i_seq];
    }
    if (vmax == vmin) {
      vmax = 1.0;
      vmin = 0.0;
    }
    // do the table lookup of rgb colours by mapping data values to indices of colourmap array
    // and attenutate colours according to attenuation values
    if (map_directly)
    {
      for (unsigned i_seq = 0; i_seq < data_for_colors.size(); i_seq++)
      {
        if (boost::math::isfinite(data_for_colors[i_seq]))
        {
          int indx = int(data_for_colors[i_seq]) % nrgbcolours;
          colors[i_seq][0] = colourmap[indx][0] * attenuation[i_seq] + 0.5 * (1.0 - attenuation[i_seq]);
          colors[i_seq][1] = colourmap[indx][1] * attenuation[i_seq] + 0.5 * (1.0 - attenuation[i_seq]);
          colors[i_seq][2] = colourmap[indx][2] * attenuation[i_seq] + 0.5 * (1.0 - attenuation[i_seq]);
        }
        else
        {
          colors[i_seq][0] = FP_NAN;
          colors[i_seq][1] = FP_NAN;
          colors[i_seq][2] = FP_NAN;
        }
      }
    }
    else
    {
      for (unsigned i_seq = 0; i_seq < data_for_colors.size(); i_seq++)
      {
        if (boost::math::isfinite(data_for_colors[i_seq]))
        {
          // applying a power scaling skews the colour mapping towards the smaller numbers for powscale > 1
          // but skews it towards the larger numbers for 0 < powscale < 1
          double colindex = (nrgbcolours - 1) * pow(((data_for_colors[i_seq] - vmin) / (vmax - vmin)), powscale);
          int indx = int(colindex);
          SCITBX_ASSERT(indx >= 0 && indx < nrgbcolours);
          colors[i_seq][0] = colourmap[indx][0] * attenuation[i_seq] + 0.5 * (1.0 - attenuation[i_seq]);
          colors[i_seq][1] = colourmap[indx][1] * attenuation[i_seq] + 0.5 * (1.0 - attenuation[i_seq]);
          colors[i_seq][2] = colourmap[indx][2] * attenuation[i_seq] + 0.5 * (1.0 - attenuation[i_seq]);
        }
        else
        {
          colors[i_seq][0] = FP_NAN;
          colors[i_seq][1] = FP_NAN;
          colors[i_seq][2] = FP_NAN;
        }
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
    double vmax = -9e99;
    double vmin = 9e99;
    for (unsigned i_seq = 0; i_seq < properties.size(); i_seq++) {
      if ((! shade_all) && (! selection[i_seq])) continue;
      if ( !boost::math::isfinite(properties[i_seq]) ) continue;
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
