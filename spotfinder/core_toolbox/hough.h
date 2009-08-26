/* Prototypes */

#ifndef HOUGH_H
#define HOUGH_H

#include <string>
#include <iostream>
#include <algorithm>
#include <scitbx/vec3.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace spotfinder {

  class hough {

  public:
    void importData(scitbx::af::flex_int,double);    // import raw data
    void exportData(scitbx::af::flex_int);           // export image data
    void setGeometry(double,double,double,double);   // set detector geometry
    scitbx::af::shared<double> getGeometry();        // get detector geometry
    void gaussianBlur(int);                          // apply Gaussian blur
    void sobelEdge();                                // apply Sobel operator
    void cannyEdge(int,double,double);               // Canny edge detection
    void findEllipse(int,double,double);             // Hough transformation
    double getDistance(double,double,double,double,double,double);
    scitbx::af::shared<double> getRings();           // returns found rings

    /*
      the input for image filters (gaussianBlur, sobelEdge, cannyEdge)
      is always rawData, and their output is always newData
    */

  private:
    scitbx::af::shared<int> origData;            // storage for original data
    scitbx::af::shared<int> rawData;             // storage for raw data
    scitbx::af::shared<int> newData;             // storage for filtered image
    scitbx::af::shared<int> tempData;            // temporary storage

    int row, col;                                // x,y dimensions of image
    double beamX, beamY;                         // center of beam (mm)
    double distance;                             // distance to detector
    double twoTheta;                             // two theta angle
    double pixelSize;                            // pixel size

    int width;                                   // width of rings
    double threshold;                            // threshold cutoff for hough
    double nStdDev;                              // standard deviation cutoff

    scitbx::af::shared<scitbx::af::shared<int> > sampleXY;// x,y of sample points
    scitbx::af::shared<scitbx::af::shared<double> > rings;// rings (7 col)
    scitbx::af::shared<scitbx::af::shared<int> > // sample points in each ring
      ringPoints;
  };
}

#endif  // hough header
