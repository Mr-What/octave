/*
Sub-pixel edge detection, based on some of the ideas
in the Canny edge detector.

Follows ridges of the gradient^2 image, in the direction
of the perpendicular to the gradient.

It was found that the smoothing provided by the Sobel
gradient operator was desirable for making smooth edges.
Extracting gradients from image splines proved to
be problematic due to false edges and ridges in the
cubic splines.

*/

typedef struct {float x,y;} EdgePoint;
typedef struct {
  int nLines;
  EdgePoint **line;
  int *lineLen;
} EdgeLines;

/** Call this to free resources in an EdgeLines struct */
void freeEdgeLines(EdgeLines *);

/** Get list of points forming edge lines in an image */
EdgeLines subCannyEdge(
     const float *img, /**< Image buffer, raster order */
     const int nx, /**< Number of pixels per scan */
     const int ny, /**< Number of raster scans */
     const float step, /**< Target spacing, in pixels, for edge line samples */
     const float res, /**< Resolution of edge points (pixels) (<<step) */
     const float hiThresh, /**< Minimum normalized strength of gradient^2 to be used for a seed in an edge line */
     const float loThresh); /**< Minimem normalized strength of gradient^2 to be considered part of an edge line */

#ifdef _DEBUG
int bufCleanup();  /* debug garbage clean-up */
#endif

/***********************************************************************
2015-06-09  ab   octave port
2007-10-25  ab   passed icpDemo.m test.
2007-10-23  ab
*/
