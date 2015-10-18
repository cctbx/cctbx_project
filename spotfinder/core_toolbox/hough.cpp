/*
  Use Hough transform to detect ice rings (ellipses)

  The Hough transform is broken down into the following series of steps:
    1) Gaussian Blur - removes noise from the image
    2) Sobel Filter - finds edges in the image based on the gradient of pixel
                      intensities (large changes from light to dark indicate
                      an edge)
    3) Remove weak edges (gradient is below a certain threshold)
    4) Check that each pixel in an edge has one or two  neighbors that are also
       edges (real edges are continuous and move in one direction)
    5) The remaining edge pixels are set to a constant and all other pixels are
       set to zero.
    (Steps 1-5 is also known as Canny Edge detection.)
    5) Generate a series of ellipses.  For ice ring detection, the ellipses are
       concentric, centered at the beam center, and parameterized by only one
       parameter, the radius of the ice ring when the two theta angle is zero.
       Each ellipse is actually a ring with a finite width
    6) For each elliptical ring, count the number of edge pixels within it.
       The ellipses with the highest counts represent the ice rings.  The
       radius is converted to a resolution.

  The Canny Edge detection algorithm should mark the pixels that belong to an
  ice ring.  The rest of the Hough transforms finds the ellipses that best fits
  those pixels.
*/

#include "spotfinder/core_toolbox/hough.h"
using namespace spotfinder;
  /* ==========================================================================
     store data in array
       data: flex_int array containing image data
       r: int containing number of rows in the image
       c: int containing number of columns in the image
       p: double containing the pixel size
     --------------------------------------------------------------------------
  */
  void hough::importData(scitbx::af::flex_int data,double p) {
    // initialize arrays and image size parameters
    origData = data.deep_copy();
    rawData = data.deep_copy();
    newData = scitbx::af::shared<int>(origData.size(),0);
    tempData = scitbx::af::shared<int>(origData.size(),0);
    row = data.accessor().focus()[1];
    col = data.accessor().focus()[0];
    pixelSize = p;
  }

  /* ==========================================================================
     export data
       data: flex_int array to be exported to
     --------------------------------------------------------------------------
  */
  void hough::exportData(scitbx::af::flex_int data) {
    for (int i=0; i<newData.size(); i++) {
      data[i] = newData[i];
    }
  }

  /* ==========================================================================
     sets the geometry of the detector
       x: double containing the x coordinate of the center (in mm)
       y: double containing the y coordinate of the center (in mm)
       dist: double containing the detector distance (in mm)
       t: double containing the 2 theta angle (in degrees)
     --------------------------------------------------------------------------
  */
  void hough::setGeometry(double x,double y,double dist,double t) {

    distance = dist;
    twoTheta = t;

    // beam center specific to detector
    beamX = y;
    beamY = ( (double(row)*pixelSize - x) -
              dist * std::tan(scitbx::deg_as_rad(twoTheta)) );
  }

  /* ==========================================================================
     returns the detector geometry as a shared array
     --------------------------------------------------------------------------
  */
  scitbx::af::shared<double> hough::getGeometry() {

    scitbx::af::shared<double> geometry(4);
    geometry[0] = ( double(row)*pixelSize - beamY -
                    distance * std::tan(scitbx::deg_as_rad(twoTheta)) );
    geometry[1] = beamX;
    geometry[2] = distance;
    geometry[3] = twoTheta;

    return geometry;
  }

  /* ==========================================================================
     apply Gaussian blur to image
       pixelRadius: int containing the radius for the Gaussian blur
     --------------------------------------------------------------------------
  */
  void hough::gaussianBlur(int pixelRadius) {

    // calculate Gaussian filter using Pascal's triangle
    scitbx::af::shared<int> gaussian;
    int filterLength = 2*pixelRadius + 1;
    for (int y=0; y<filterLength; y++) {
      int c = 1;
      for (int x=0; x<=y; x++) {
        if (y == filterLength-1) {
          gaussian.push_back(c);
        }
        c = c * (y-x) / (x+1);
      }
    }

    // calculate sum of gaussian filter for later rescaling
    int sum = 0;
    for (int i=0; i<gaussian.size(); i++) {
      sum = sum + gaussian[i];
    }

    // loop over image
    for (int i=pixelRadius; i<row-pixelRadius; i++) {
      for (int j=pixelRadius; j<col-pixelRadius; j++) {

        // loop over Gaussian filter (in x direction)
        for (int k=-pixelRadius; k<pixelRadius+1; k++) {
          tempData[i*row + j] = tempData[i*row + j] +
            rawData[i*row + (j+k)] * gaussian[k+pixelRadius];
        }

        // rescale pixel value
        tempData[i*row + j] = tempData[i*row + j] / sum;
      }
    }

    // repeat in y direction
    for (int j=pixelRadius; j<col-pixelRadius; j++) {
      for (int i=pixelRadius; i<row-pixelRadius; i++) {
        for (int k=-pixelRadius; k<pixelRadius+1; k++) {
          newData[i*row + j] = newData[i*row + j] +
            tempData[(i+k)*row + j] * gaussian[k+pixelRadius];
        }
        newData[i*row + j] = newData[i*row + j] / sum;
      }
    }
  }

  /* ==========================================================================
     find edges using Sobel operator
       Applies the x and y operators to the image and stores the edge gradient
       strength in newData and the edge direction in tempData.
     --------------------------------------------------------------------------
  */
  void hough::sobelEdge() {

    // Sobel operator
    int sobelX[3][3] = { { 1, 0,-1},
                         { 2, 0,-2},
                         { 1, 0,-1} };
    int sobelY[3][3] = { { 1, 2, 1},
                         { 0, 0, 0},
                         {-1,-2,-1} };

    // loop over image
    for (int i=3; i<row-3; i++) {
      for (int j=3; j<col-3; j++) {

        // loop over x mask
        for (int ii=-1; ii<2; ii++) {
          for (int jj=-1; jj<2; jj++) {
            tempData[i*row + j] = tempData[i*row + j] +
              rawData[(i+ii)*row + (j+jj)] * sobelX[ii+1][jj+1];
          }
        }

        // loop over y mask
        for (int ii=-1; ii<2; ii++) {
          for (int jj=-1; jj<2; jj++) {
            newData[i*row + j] = newData[i*row + j] +
              rawData[(i+ii)*row + (j+jj)] * sobelY[ii+1][jj+1];
          }
        }

        // classify edge
        int theta;
        double tempTheta;

        if (tempData[i*row + j] == 0) {
          if (newData[i*row + j] == 0) {
            theta = 0;
          }
          else {
            theta = 90;
          }
        }
        else {
          tempTheta = std::atan(double(newData[i*row + j]) /
                           double(tempData[i*row + j]));
          tempTheta = scitbx::rad_as_deg(tempTheta);
          if (tempTheta <= 22.5 && tempTheta > -22.5) {
            theta = 0;
          }
          else if (tempTheta > 22.5 && tempTheta <= 67.5) {
            theta = 45;
          }
          else if (tempTheta > 67.5 || tempTheta <= -67.5) {
            theta = 90;
          }
          else {
            theta = 135;
          }
        }

        /* combine results
           newData stores edge gradient strength
           tempData store edge direction
        */
        newData[i*row + j] = ( std::abs(newData[i*row + j]) +
                               std::abs(tempData[i*row + j]) );
        tempData[i*row + j] = theta;
      }
    }
  }

  /* ==========================================================================
     edge detection using Canny method
       1) reduce noise by applying Gaussian filter
       2) find edges using Sobel operator
       3) highlight edges with non-maximum suppression
       4) clean up edges using thresholds (hysteresis)
       5) clean up noise by removing single pixels and clusters of pixels

       pixelRadius: the radius for the Gaussian blur
       f1: the lower fraction of pixels to be excluded (<1.0)
       f2: the higher fraction of pixels to be excluded (<1.0, f2 > f1)
     --------------------------------------------------------------------------
  */
  void hough::cannyEdge(int pixelRadius,double f1,double f2) {

    // apply Gaussian filter
    hough::gaussianBlur(pixelRadius);

    // replace old image with new image and reset newData and tempData
    for (int i=0; i<rawData.size(); i++) {
      rawData[i] = newData[i];
      newData[i] = 0;
      tempData[i] = 0;
    }

    // find edges using Sobel operator
    hough::sobelEdge();

    // replace old image with new image
    for (int i=0; i<rawData.size(); i++) {
      rawData[i] = newData[i];
      newData[i] = 0;
    }

    // calculate image statistics for thresholding
    scitbx::af::shared<int> sortMe(rawData.size(),0);
    for (int i=0; i<sortMe.size(); i++) {
      sortMe[i] = rawData[i];
    }
    std::sort(sortMe.begin(),sortMe.end());
    int min = sortMe.front();
    int max = sortMe.back();
    int t1 = sortMe[int(f1*sortMe.size())];
    int t2 = sortMe[int(f2*sortMe.size())];

    // remove highest peaks
    for (int i=0; i<rawData.size(); i++) {
      if (rawData[i] > t2) {
        rawData[i] = 0;
      }
    }

    // highlight edges with non-maximum suppression
    int filterRadius = pixelRadius + 1;
    for (int i=filterRadius; i<row-filterRadius; i++) {
      for (int j=filterRadius; j<col-filterRadius; j++) {

        // check gradient strength along edge direction
        if (rawData[i*row + j] < t1) { // remove if less than t1
          newData[i*row + j] = 0;
        }
        else {
          if (tempData[i*row + j] == 0) { // check horizontal line
            if (rawData[i*row + j] > rawData[i*row + (j+1)] &&
                rawData[i*row + j] > rawData[i*row + (j-1)]) {
              newData[i*row + j] = max;
            }
          }
          else if (tempData[i*row + j] == 45) { // check diagonal line (45)
            if (rawData[i*row + j] > rawData[(i+1)*row + (j+1)] &&
                rawData[i*row + j] > rawData[(i-1)*row + (j-1)]) {
              newData[i*row + j] = max;
            }
          }
          else if (tempData[i*row + j] == 90) { // check vertical line
            if (rawData[i*row + j] > rawData[(i+1)*row + j] &&
                rawData[i*row + j] > rawData[(i-1)*row + j]) {
              newData[i*row + j] = max;
            }
          }
          else if (tempData[i*row + j] == 135) { // check diagonal line (135)
            if (rawData[i*row + j] > rawData[(i-1)*row + (j+1)] &&
                rawData[i*row + j] > rawData[(i+1)*row + (j-1)]) {
              newData[i*row + j] = max;
            }
          }
          else {  // use hysteresis to determine if edge is real
            if (rawData[i*row + j] > t2) {  // make edge if larger than t2
              newData[i*row + j] = max;
            }
            else {  // check if neighbors are edges if in between
              if (rawData[(i-1)*row + (j-1)] > t2 ||
                  rawData[(i-1)*row + (j+0)] > t2 ||
                  rawData[(i-1)*row + (j+1)] > t2 ||
                  rawData[(i+0)*row + (j-1)] > t2 ||
                  rawData[(i+0)*row + (j+1)] > t2 ||
                  rawData[(i+1)*row + (j-1)] > t2 ||
                  rawData[(i+1)*row + (j+0)] > t2 ||
                  rawData[(i+1)*row + (j+1)] > t2) {
                newData[i*row + j] = max;
              }
              else {
                newData[i*row + j] = 0;
              }
            }
          }
        }
      }
    }

    /*
      remove excessively connected points (3 points within sqrt(2) distance)
      and points with no connections (no neighbors within sqrt(2) distance)
    */
    for (int i=filterRadius; i<row-filterRadius; i++) {
      for (int j=filterRadius; j<col-filterRadius; j++) {
        int count = 0;
        // count neighbors around pixel
        if (newData[i*row + (j-1)] != 0) {
          count = count + 1;
        }
        else if (newData[i*row + (j+1)] != 0) {
          count = count + 1;
        }
        else if (newData[(i-1)*row + j] != 0) {
          count = count + 1;
        }
        else if (newData[(i-1)*row + j] != 0) {
          count = count + 1;
        }
        // check number of neighbors
        if (count == 0 || count > 2) {
          newData[i*row + j] = 0;
        }
      }
    }

    // store non-zero points separately
    for (int i=0; i<row; i++) {
      for (int j=0; j<col; j++) {
        if (newData[i*row + j] != 0) {
          scitbx::af::shared<int> p(2);
          p[0] = j;
          p[1] = i;
          sampleXY.push_back(p);
        }
      }
    }
  }

  /* ==========================================================================
     find ellipses using Hough transformation

       The equation for an ellipse is given by,

            [( x cos(theta) + y sin(theta)) - x0]^2
            --------------------------------------- +
                              a^2

            [(-x sin(theta) + y cos(theta)) - y0]^2
            --------------------------------------- = 1
                              b^2

       where (x0,y0) denotes the center of the ellipse, (x,y) denotes a point
       on the ellipse, theta is the skew, and a and b are the major/minor axes
       of the ellipse.  Since the center is known (the beam center) and there
       is no skew, there are two remaining parameters for the Hough transform,
       a and b.  Furthermore, since the diffraction image at a non-zero two
       theta angle is the projection of the image at zero onto the detector at
       non-zero two theta, a and b are not independent, which reduces the
       number of parameters to one, the radius of the circle prior to
       projection.

       Hough transform procedure:
       1) generate ellipses for a set of a and b.
       2) for each sample point, find all ellipses that contain it
       3) pick ellipses with highest numbers of points

     --------------------------------------------------------------------------
  */
  void hough::findEllipse(int width,double threshold,double nStdDev) {

    rings.resize(0);
    ringPoints.resize(0);

    // calculate beam center
    int x0 = int(beamX / pixelSize);
    int y0 = row - int(beamY / pixelSize);

    /*
      calculate projection of a circle with radius, r, onto plane offset by
      twoTheta
    */

    scitbx::af::shared<double> r(row/width + 1);  // radius (hough space)
    scitbx::af::shared<double> b(r.size());       // transformed r
    scitbx::af::shared<int> c(r.size());          // transformed center of ellipse
    scitbx::af::shared<int> houghSpace(r.size(),0);
    SCITBX_ASSERT(r.size() > 0);
    /*
      initialize values
           r (set of radii)
           b (transformed major axis)
           c (transformed center of ellipse)
    */
    double h1, h2, alpha;

    /* build concentric circles */
    r[0] = width;
    for (int i=1; i<r.size(); i++) {
      r[i] = r[i-1] + width;
    }
    scitbx::af::shared<double> r0 = r.deep_copy();

    for (int i=0; i<r.size(); i++) {
      // calculate major axis of ellipse
      alpha = std::atan((r[i]*pixelSize)/distance);
      h1 = ( distance * std::sin(alpha) /
             ( std::sin(scitbx::deg_as_rad(90.0 - twoTheta)) *
               std::cos(scitbx::deg_as_rad(twoTheta) - alpha) ) );
      h2 = ( distance * std::cos(scitbx::deg_as_rad(90.0) - alpha) /
             ( std::cos(scitbx::deg_as_rad(twoTheta)) *
               std::cos(scitbx::deg_as_rad(twoTheta) + alpha) ) );
      b[i] = (h1 + h2) / 2;

      // calculate minor axis of ellipse
      r[i] =(r[i]/distance)*(distance/std::cos(scitbx::deg_as_rad(twoTheta)) +
                             (h2-b[i])*std::sin(scitbx::deg_as_rad(twoTheta)));

      // calculate center of ellipse
      c[i] = row - int((beamY - (h2-b[i])) / pixelSize);

      // rescale b[i] to pixels
      b[i] = b[i] / pixelSize;
    }

    /*
      determine which points belong to which ellipse
      each circle of radius r, describes a unique ellipse after projection.
      points between ellipses of radius r[i] and r[i+1] belong to the ellipse
      of radius r[i]+width/2.
    */
    double x2, y2, a20, b20, a2, b2;
    for (int i=0; i<r.size()-1; i++) {

      scitbx::af::shared<int> points;

      // loop over sample points
      for (int j=0; j<sampleXY.size(); j++) {
        x2 = sampleXY[j][0] - x0;
        x2 = x2 * x2;
        y2 = sampleXY[j][1] - c[i];
        y2 = y2 * y2;
        a20 = r[i] * r[i];
        b20 = b[i] * b[i];
        a2 = r[i+1] * r[i+1];
        b2 = b[i+1] * b[i+1];
        if ( (x2/a20 + y2/b20) >= 1.0 &&  // larger or equal to r[i],b[i]
             (x2/a2 + y2/b2) < 1.0) {     // smaller than r[i+1],b[i+1]
          houghSpace[i] = houghSpace[i] + 1;
          points.push_back(sampleXY[j][0]);
          points.push_back(sampleXY[j][1]);
        }
      }
      ringPoints.push_back(points);
    }

    /*
      bin maximums
      the image is broken into radial bins. maximums higher than a certain
      standard deviation (nStdDev) above the bin average or above some
      magnitude (threshold) are considered rings
      nStdDev is used to find weak rings and threshold is used to find rings
      in noisy regions
    */
    int bins = 10;
    int binSize = r.size()/bins;
    scitbx::af::shared<double> average(bins,0.0);
    scitbx::af::shared<double> stdDev(bins,0.0);

    for (int i=0; i<bins; i++) {
      // calculate average
      for (int j=0; j<binSize; j++) {
        average[i] = average[i] + houghSpace[i*binSize + j];
      }
      average[i] = average[i] / binSize;

      // calculate standard deviation
      for (int j=0; j<binSize; j++) {
        double temp = houghSpace[i*binSize + j] - average[i];
        stdDev[i] = stdDev[i] + temp * temp;
      }
      stdDev[i] = std::sqrt(stdDev[i] / binSize);

      // compare values in current bin with bin average
      for (int j=0; j<binSize; j++) {
        if (houghSpace[i*binSize + j] > threshold ||
            (houghSpace[i*binSize + j] > (nStdDev*stdDev[i] + average[i]) &&
             houghSpace[i*binSize + j] > 100) ) {
          scitbx::af::shared<double> points(7);
          points[0] = i*binSize + j;
          points[1] = houghSpace[i*binSize + j];
          points[2] = x0;
          points[3] = c[i*binSize + j];
          points[4] = r[i*binSize + j] + width/2;
          points[5] = b[i*binSize + j] + width/2;
          points[6] = r0[i*binSize + j] + width/2;
          rings.push_back(points);
        }
      }
    }

    /*
      if the width of the ice ring is larger than the width of the elliptical
      rings used in the Hough transform, there will be several elliptical rings
      that belong to one ice ring.
    */
    if (rings.size() > 1) {
      // combine rings
      scitbx::af::shared<scitbx::af::shared<scitbx::af::shared<double> > >
        group;
      scitbx::af::shared<scitbx::af::shared<double> > manyRings;
      for (int i=0; i<rings.size()-1; i++) {
        scitbx::af::shared<double> points(7);
        for (int j=0; j<points.size(); j++) {
          points[j] = rings[i][j];
        }
        manyRings.push_back(points);
        // if next ring is far, make a new group
        if ((rings[i+1][0] - rings[i][0]) > 5) {
          scitbx::af::shared<scitbx::af::shared<double> > tRings;
          for (int k=0; k<manyRings.size(); k++) {
            scitbx::af::shared<double> tPoints(7);
            for (int l=0; l<tPoints.size(); l++) {
              tPoints[l] = manyRings[k][l];
            }
            tRings.push_back(tPoints);
          }
          group.push_back(tRings);
          manyRings.resize(0);
        }
      }
      // check if all rings belong to same group
      if (manyRings.size() != 0) {
        group.push_back(manyRings);
      }
      // check last ring
      if ((rings[rings.size()-1][0] - rings[rings.size()-2][0]) > 5) {
        scitbx::af::shared<scitbx::af::shared<double> > tRings;
        scitbx::af::shared<double> tPoints(7);
        for (int i=0; i<tPoints.size(); i++) {
          tPoints[i] = rings[rings.size()-1][i];
        }
        tRings.push_back(tPoints);
        group.push_back(tRings);
      }
      else { // merge with previous group
        scitbx::af::shared<double> tPoints(7);
        for (int i=0; i<tPoints.size(); i++) {
          tPoints[i] = rings[rings.size()-1][i];
        }
        group.back().push_back(tPoints);
      }
      // calculate weighted average based on score
      rings.resize(0);
      for (int i=0; i<group.size(); i++) {
        scitbx::af::shared<double> points(7,0.0);
        int sum = 0;
        for (int j=0; j<group[i].size(); j++) {
          sum = sum + group[i][j][1];
          for (int k=0; k<group[i][j].size(); k++) {
            points[k] = points[k] + group[i][j][1]*group[i][j][k];
          }
        }
        for (int j=0; j<points.size(); j++) {
          points[j] = points[j]/sum;
        }
        points[1] = sum/group[i].size();
        rings.push_back(points);
      }
    }
  }

  /* ==========================================================================
     given a point,p(u,v), and an ellipse,e(x0,y0,a,b), calculate the shortest
     distance between the point and the ellipse
     --------------------------------------------------------------------------
  */
  double hough::getDistance(double u,double v,double x0,
                                     double y0,double a,double b) {
    /*
      the elllipse is parameterized as
         x = a cos(theta) + x0
         y = b cos(theta) + y0
         r = (x,y)

      the shortest distance between r and p,|p - r|, must satisfy the condition
      that |p - r| be perpendicular to the tangent of r,r'
         (p - r).r' = 0
         f(theta) = (a^2 - b^2) sin(theta) cos(theta) -
                    (u - x0) a sin(theta) +
                    (v - y0) b cos(theta) = 0
         f'(theta) = (a^2 - b^2)(cos^2(theta) - sin^2(theta)) -
                     (u - x0) a cos(theta) -
                     (v - y0) b sin(theta)
    */

    double a2_b2 = a*a - b*b;
    double u_x0 = u - x0;
    double v_y0 = v - y0;
    double theta = std::atan(v_y0/u_x0);
    if (u_x0 < 0.0 && v_y0 >= 0.0) {
      theta = theta + scitbx::constants::pi;
    }
    if (u_x0 < 0.0 && v_y0 < 0.0) {
      theta = theta - scitbx::constants::pi;
    }

    // use Newton's method to find root
    double f = ( a2_b2 * std::sin(theta) * std::cos(theta) -
                 u_x0 * a * std::sin(theta) +
                 v_y0 * b * std::cos(theta) );
    double df = ( a2_b2 * (std::cos(theta)*std::cos(theta) -
                           std::sin(theta)*std::sin(theta)) -
                  u_x0 * a * std::cos(theta) -
                  v_y0 * b * std::sin(theta) );

    for (int i=0; i<5; i++) {
      theta = theta - f/df;
      f = ( a2_b2 * std::sin(theta) * std::cos(theta) -
            u_x0 * a * std::sin(theta) +
            v_y0 * b * std::cos(theta) );
      df = ( a2_b2 * (std::cos(theta)*std::cos(theta) -
                      std::sin(theta)*std::sin(theta)) -
             u_x0 * a * std::cos(theta) -
             v_y0 * b * std::sin(theta) );
    }

    // calculate distance, |p-r|
    double dx = u - a * std::cos(theta) - x0;
    double dy = v - b * std::sin(theta) - y0;
    double distance = std::sqrt(dx*dx + dy*dy);
    return distance;
  }

  /* ==========================================================================
     returns a vector containing the rings
     --------------------------------------------------------------------------
  */

  scitbx::af::shared<double> hough::getRings() {
    int row = rings.size();
    int col = 7;
    scitbx::af::shared<double> rings_1d(row*col);
    for (int i=0; i<row; i++) {
      for (int j=0; j<col; j++) {
        rings_1d[i*col+j] = rings[i][j];
      }
    }
    return rings_1d;
  }
