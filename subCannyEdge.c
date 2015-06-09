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

#include "subCannyEdge.h"
#include <stdlib.h> /* calloc */
#include <stdio.h>
#include <string.h> /* memset */
#include <math.h>

#ifdef _DEBUG
const static int BUF_CHUNK = 1024;
typedef struct {
  void **p;
  int maxLen,n;
} BufferList;
static BufferList bufList = {NULL,0,0};

/*static void bufReport(){}*/

static void bufPush(void *p)
{
  if (bufList.n >= bufList.maxLen)
    {
      void *q;
      if (bufList.maxLen <= 0) bufList.maxLen  = 2*BUF_CHUNK;
      else                     bufList.maxLen +=   BUF_CHUNK;
      q = bufList.p;
      bufList.p = (void **)realloc(q,bufList.maxLen*sizeof(void *));
      if (p==NULL)
        {
          fprintf(stderr,"Could not grow pointer buffer to %d.\n",
                  bufList.maxLen);
          exit(3);
        }
      fprintf(stderr,"buffer list moved from %p to %p and grown to %d\n",
              q,bufList.p,bufList.maxLen);
    }
  bufList.p[bufList.n++] = p;
}

static void bufPop(void *p)
{
  int k;
  for (k=0;k<bufList.n;k++)
    {
      if (p == bufList.p[k])
        {
          bufList.p[k]=bufList.p[--bufList.n];
          return;
        }
    }
  fprintf(stderr,"Attempt to free unlisted buffer at %p\n",p);
}

/** return number of stranded buffers */
int bufCleanup()
{
  int n = bufList.n;
  fputs("Garbage clean-up:\n",stderr);
  while (bufList.n>0)
    {
      fprintf(stderr,"\t%d) %p\n",bufList.n,bufList.p[n-1]);
      free(bufList.p[--bufList.n]);
    }
  free(bufList.p);
  fprintf(stderr,"Garbage list at %p released\n",bufList.p);
  bufList.p=NULL;
  bufList.maxLen=0;
  return(n);
}
#endif


static void myFree(void *p)
{
#ifdef _DEBUG_VERBOSE
  fprintf(stderr,"Releasing %p\n",p);fflush(stderr);
  bufPop(p);
#endif
  free(p);
}

static void *myCalloc(int nItems, int itemSize)
{
  void *p;
  p = calloc(nItems,itemSize);
#ifdef _DEBUG_VERBOSE
  bufPush(p);
  fprintf(stderr,"calloc(%d,%d)=%p\n",nItems,itemSize,p);
#endif
  if (p==NULL)
    fprintf(stderr,"Memory error on calloc(%d,%d)\n",nItems,itemSize);
  return(p);
}

static void *myRealloc(void *p, int nBytes)
{
  void *q;
  q = realloc(p,nBytes);
#ifdef _DEBUG_VERBOSE
  bufPop(p);
  bufPush(q);
  fprintf(stderr,"realloc(%p,%d)=%p\n",p,nBytes,q);
#endif
  if (q==NULL)
    fprintf(stderr,"Memory error on realloc(%p,%d)\n",p,nBytes);
  return(q);
}


/*-------------------------------------------- EdgeLines struct stuff */

/* grow number of edge lines buffer by this many */
const static int EDGE_LINE_CHUNK=100;

/* Call this to free resources in an EdgeLines struct */
void freeEdgeLines(EdgeLines *e)
{
  int k;

  k = e->nLines-1;
  /*k--;  debug only...  test mem err */
  for (;k>=0;k--)
    {
#ifdef _DEBUG_VERBOSE
      fprintf(stderr,"line %d, %d points\n",k+1,e->lineLen[k]);fflush(stderr);
#endif
      myFree(e->line[k]);
      e->line[k]=NULL;
    }
#ifdef _DEBUG_VERBOSE
  fputs("line pointers: ",stderr);fflush(stderr);
#endif
  myFree(e->line);  e->line=NULL;
#ifdef _DEBUG_VERBOSE
  fputs("line lengths: ",stderr);fflush(stderr);
#endif
  myFree(e->lineLen);  e->lineLen=NULL;
  e->nLines=0;
}

void initEdgeLines(EdgeLines *el)
{
  int k;

  el->nLines=0;
  el->line = (EdgePoint **)myCalloc(EDGE_LINE_CHUNK,sizeof(EdgePoint *));
  if (el->line == NULL)
    {
      fputs("Out of memory allocating line pointer buffer.\n",stderr);
      return;
    }
  el->lineLen = (int *)myCalloc(EDGE_LINE_CHUNK,sizeof(int));
  if (el->lineLen == NULL)
    {
      fputs("Out of memory allocating line length buffer.\n",stderr);
      return;
    }
  for(k=0;k<EDGE_LINE_CHUNK;k++)
    {
      el->line[k]=NULL;
      el->lineLen[k]=0;
    }
}

/* re-alloc line buffers if vectors are full
   return 0 or error code */
int growEdgeLines(EdgeLines *el)
{
  int m;

  if (el->nLines < EDGE_LINE_CHUNK) return(0);
  m = el->nLines % EDGE_LINE_CHUNK;
  if (m > 0) return(0);
  m = el->nLines / EDGE_LINE_CHUNK;
  m++;
  m *= EDGE_LINE_CHUNK;
  el->lineLen = (int *)myRealloc(el->lineLen,m*sizeof(int));
  el->line = (EdgePoint **)myRealloc(el->line,m*sizeof(EdgePoint *));
  if ((el->lineLen == NULL) | (el->line == NULL))
    {
      fputs("Out of memory re-allocating line buffers",stderr);
      return(-1);
    }
  return(0);
}

/*-------------------------------------------- Internal sub-functions */

void write_float(const float *x, const int n, const char *fNam)
{
  FILE *f;
  f = fopen(fNam,"wb");
  fwrite(x,sizeof(float),n,f);
  fclose(f);
}

/** Low-level, generic, sobel gradient filter.
    Also returns greatest abs(gx,gy)

    return max gradient, or negative error code
*/
float sobelGradient(const float *img,
               const int nx, const int ny,
               float *gx,
               float *gy)
{
  int y,yHi,xHi;
  float hi;
  const float *r0,*rp;

  /* just set edges to zero */
  for (y=0; y<ny; y++)
    {
      int y0,y1;

      y0=y*nx;y1=y0+nx-1;
      gx[y0]=gx[y1]=gy[y0]=gy[y1]=0.0f;
    }
  memset(gx,0,nx*sizeof(float));
  memset(gy,0,nx*sizeof(float));
  memset(gx+(y-1)*nx,0,nx*sizeof(float));
  memset(gy+(y-1)*nx,0,nx*sizeof(float));

  hi=0;
  yHi=ny-1;
  xHi=nx-1;
  r0=img;rp=r0+nx;
  for (y=1;y<yHi;y++)
    {
      const float *rm;
      float *rx,*ry;
      int x;

      rx = gx + y*nx;
      ry = gy + y*nx;
      rm = r0;
      r0 = rp;
      rp += nx;
      for (x=1;x<xHi;x++)
        {
          float g;
          g = (float)(   ((double)rp[x+1]) + ((double)rm[x+1]) +
                      2*(((double)r0[x+1]) - ((double)r0[x-1]))
                      -  ((double)rp[x-1]) - ((double)rm[x-1]) );
          rx[x]=g;
          if (g<0) g=-g;
          if (g>hi) hi=g;
          g = (float)(   ((double)rp[x-1]) + ((double)rp[x+1]) +
                      2*(((double)rp[x  ]) - ((double)rm[x  ]))
                      -  ((double)rm[x-1]) - ((double)rm[x+1]) );
          ry[x]=g;
          if (g<0) g=-g;
          if (g>hi) hi=g;
        }
    }
#ifdef _DEBUG
  write_float(gx,256*256,"gx.float");
  write_float(gy,256*256,"gy.float");
#endif
  return(hi);
}

void biCubicCoefs(double z[16],double c[16])
{
  /*
% Given a function, f(y,x) known at all points x,y within [-1..2],
% compute coefficients C(j,i) where i,j within [0..3] s.t.
% f(y,x) where x,y within [0,1] can be estimated as:
% f=0;
% for i=0:3
%    for j=0:3
%       f = f + C(j,i) * x^i * y^j;
%    end
% end
%
% I believe this is the traditional bicubic spline, which (i hope)
% is contiguous at least to the 1st derivative on the edges
% of the central unit cell, even as it is re-computed from cell to cell
function C = biCubicCoefs(f)
global unitBiCubicGenMat;
if (prod(size(unitBiCubicGenMat)) ~= 16*16)
   %% ran this once for answer, then just hard-code answer
   %A = zeros(16);
   %xx=A;yy=A;
   %k=1;
   %for x=-1:2
   %   for y=-1:2
   %      m = 1;
   %      for i = 0:3
   %         xi=x^i;
   %         for j = 0:3
   %            A(k,m) = (y^j) * xi;
   %            xx(k,m) = x;  yy(k,m) = y;
   %            fprintf(1,'A(%2d,%2d) = %2d^%d * %2d^%d = %d\n',k,m,x,i,y,j,A(k,m));
   %            m=m+1;
   %         end
   %      end
   %      k=k+1;
   %   end
   %end
   %A
   %unitBiCubicGenMat = inv(A);
  */
  static double genMat[16][16] = {
    {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,-0.33333333333333333,-0.5,1,-0.16666666666666667,0,0,0,0,0,0,0,0},
    {0,0,0,0,0.5,-1,0.5,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,-0.16666666666666667,0.5,-0.5,0.16666666666666667,0,0,0,0,0,0,0,0},
    {0,-0.33333333333333333,0,0,0,-0.5,0,0,0,1,0,0,0,-0.166666666666666667,0,0},
    {0.111111111111111111,0.1666666666666666667,-0.333333333333333333,0.0555555555555555555,0.16666666666666666667,0.25,-0.5,0.083333333333333333333,-0.3333333333333333333,-0.5,1,-0.166666666666667,0.0555555555555555,0.0833333333333333,-0.166666666666666667,0.0277777777777777778},
    {-0.16666666666666667,0.333333333333333333,-0.1666666666666666667,0,-0.25,0.5,-0.25,0,0.5,-1,0.5,0,-0.083333333333333333333,0.166666666666666666667,-0.0833333333333333,0},
    {0.0555555555555555555556,-0.1666666666666666667,0.1666666666666666667,-0.05555555555555555556,0.083333333333333333333,-0.25,0.25,-0.08333333333333333333,-0.1666666666666666667,0.5,-0.5,0.1666666666666666667,0.02777777777777777778,-0.08333333333333333333,0.08333333333333333333,-0.02777777777777777778},
    {0,0.5,0,0,0,-1,0,0,0,0.5,0,0,0,0,0,0},
    {-0.1666666666666666667,-0.25,0.5,-0.0833333333333333333,0.333333333333333333,0.5,-1,0.1666666666666666667,-0.1666666666666666667,-0.25,0.5,-0.083333333333333333,0,0,0,0},
    {0.25,-0.5,0.25,0,-0.5,1,-0.5,0,0.25,-0.5,0.25,0,0,0,0,0},
    {-0.083333333333333333,0.25,-0.25,0.08333333333333333333,0.166666666666666667,-0.5,0.5,-0.166666666666666667,-0.0833333333333333333,0.25,-0.25,0.083333333333333333,0,0,0,0},
    {0,-0.166666666666666667,0,0,0,0.5,0,0,0,-0.5,0,0,0,0.1666666666666666667,0,0},
    {0.05555555555555555556,0.083333333333333333333,-0.16666666666666666667,0.0277777777777777778,-0.166666666666666667,-0.25,0.5,-0.0833333333333333333,0.166666666666666667,0.25,-0.5,0.08333333333333333333,-0.05555555555555555556,-0.083333333333333333,0.1666666666666666667,-0.02777777777777777778},
    {-0.083333333333333333333,0.16666666666666666667,-0.083333333333333333333,0,0.25,-0.5,0.25,0,-0.25,0.5,-0.25,0,0.083333333333333333333,-0.16666666666666666667,0.083333333333333333333,0},
    {0.027777777777777777778,-0.083333333333333333333,0.083333333333333333333,-0.02777777777777777778,-0.0833333333333333333,0.25,-0.25,0.08333333333333333333,0.08333333333333333333,-0.25,0.25,-0.083333333333333333333,-0.0277777777777777778,0.0833333333333333333,-0.0833333333333333333,0.0277777777777777778}};
  int k;

  for (k=0;k<16;k++)
    {
      int m;
      double zz;
      zz=0.0;
      for (m=0;m<16;m++)
        {
          zz += genMat[k][m] * z[m];
        }
      c[k]=zz;
    }
}

double biCubicValue(double c[16], double x, double y)
{
  double x2,x3,y2,y3,z;

  x2=x*x;x3=x2*x;
  y2=y*y;y3=y2*y;
  z =   c[0]  +    y*c[1]  +    y2*c[2]  +    y3*c[3] +
     x *c[4]  + x *y*c[5]  + x *y2*c[6]  + x *y3*c[7] + 
     x2*c[8]  + x2*y*c[9]  + x2*y2*c[10] + x2*y3*c[11] +
     x3*c[12] + x3*y*c[13] + x3*y2*c[14] + x3*y3*c[15];
  return(z);
}

/*--------- Functions using Internal structure for passing image data */
typedef struct {
  float *gx,*gy;
  float gMax; /* normilization factor for above gradients */
  unsigned char *used;
  int nx,ny;
  float step,res;
} ImgInfo;

static int initImgInfo(const float *img0, const int nx, const int ny,
                ImgInfo *img)
{
  int k;
  float hi;

  img->nx = img->ny = 0;
  img->gMax=0;
  img->step = img->res = 0.0f;

  img->gx = (float *)myCalloc(nx*ny*2,sizeof(float));
  img->used = (unsigned char *)myCalloc(nx*ny,sizeof(unsigned char));
  if ((img->gx == NULL) | (img->used == NULL))
    {
      fputs("Out of memory allocating image data buffers\n",stderr);
      return(-1);
    }
  img->nx = nx;
  img->ny = ny;
  img->gy = img->gx + nx*ny;
  memset(img->used,0,nx*ny);

  hi = sobelGradient(img0,nx,ny,img->gx,img->gy);
  img->gMax = hi/8;  /* normalize sobel filter for unity gain */

  /* Normalize gradients */
  for (k=0; k<ny*nx; k++)
    {
      img->gx[k] /= hi;
      img->gy[k] /= hi;
    }
  return(0);
}
          
static int freeImgInfo(ImgInfo *img)
{
  img->nx = img->ny = 0;
  img->gMax=0;
  myFree(img->used);
  img->used=NULL;
  myFree(img->gx);
  img->gx=NULL;
  img->gy=NULL;
  return(0);
}

static void extractNeigborhood4(const float *img,
      const int nx, const int x, const int y,
      double z[16])
{
  static int idx[16]={0,4, 8,12,
                      1,5, 9,13,
                      2,6,10,14,
                      3,7,11,15};
  const float *r;
  int k,n;

  r=img+((y-1)*nx + x-1);
  n=0;
  for (k=0;k<4;k++)
    {
      int m;
      for (m=0;m<4;m++) z[idx[n++]]=r[m];
      r += nx;
    }
}

/** Compute gradient of image at given point
        returned gradient is normalized.  g^2 is returned */
static void imgGradient(const ImgInfo *img,
    const EdgePoint p,
    EdgePoint *g, float *g2)
{
  double c[16],z[16],dx,dy,gg,gx,gy;
  int x,y;

  x = (int)p.x;
  y = (int)p.y;
  if(x<1)x=1;   /* just scoot to avoid edge issue crash */
  if(y<1)y=1;
  if(x>img->nx-3)x=img->nx-3;
  if(y>img->ny-3)y=img->ny-3;
  dx = p.x - x;
  dy = p.y - y;

  extractNeigborhood4(img->gx,img->nx,x,y,z);
  biCubicCoefs(z,c);
  gx = biCubicValue(c,dx,dy);

  extractNeigborhood4(img->gy,img->nx,x,y,z);
  biCubicCoefs(z,c);
  gy = biCubicValue(c,dx,dy);

  gg = gx*gx + gy*gy;
  if (gg > 0.0)
    {
      *g2 = (float)gg;
      gg = sqrt(gg);
      g->x = (float)(gx/gg);
      g->y = (float)(gy/gg);
    }
  else
    *g2=g->x=g->y=0.0f;
}

/* Follow gradient up to nearest ridge crest */
static EdgePoint refineRidgeCrest(const ImgInfo *img,
     const float x,
     const float y,
     float step,
     EdgePoint *g0,  /* return gradient at peak */
     float *g2hi)  /* return value of g2 at this point */
{
  EdgePoint p,g,pm,pp,gm,gp;
  float g2,g2m,g2p;
  int i,nShrink,dir,nFlip;

  p.x=x;p.y=y;

  imgGradient(img,p,&g,&g2);
  if (g2 <= 0.0f)
    { /* the starting point was in a flat area! */
      g0->x=g0->y=*g2hi=0.0f;
      return(p);
    }
  g.x *= step;   g.y *= step;
  pm.x=p.x-g.x; pm.y=p.y-g.y;
  pp.x=p.x+g.x; pp.y=p.y+g.y;
  imgGradient(img,pm,&gm,&g2m);
  imgGradient(img,pp,&gp,&g2p);
  gp.x *= step;   gp.y *= step;
  gm.x *= step;   gm.y *= step;

  nShrink=dir=nFlip=0;
  i=1;
  for(;;)
    {
      if(((g2>=g2m) && (g2>=g2p)) || (i>1000) || (nShrink>7))
        {
          float d;
          *g2hi = g2;
          d = sqrt(g.x*g.x + g.y*g.y);
          g0->x=g.x/d;
          g0->y=g.y/d;
          if (g2 <= 0.0f)
            g0->x=g0->y=0.0f;
          return(p);
        }

      if (nFlip>1)
        {
          nFlip = 0;
          nShrink++;
          step *= 0.9;
        }

      if (g2m>g2)
        {
          pp=p;gp=g;g2p=g2;
          p=pm;g=gm;g2=g2m;
          pm.x = p.x - g.x;
          pm.y = p.y - g.y;
          imgGradient(img,pm,&gm,&g2m);
          gm.x *= step;   gm.y *= step;
          if (dir != -1) nFlip++;
          dir=-1;
        }
      else
        {
          pm=p;gm=g;g2m=g2;
          p=pp;g=gp;g2=g2p;
          pp.x = p.x + g.x;
          pp.y = p.y + g.y;
          imgGradient(img,pp,&gp,&g2p);
          gp.x *= step;   gp.y *= step;
          if (dir != 1) nFlip++;
          dir=1;
        }

      if (i++>800)
        {
          fprintf(stderr,"p=[%f,%f] g=[%f,%f] %.6f\n",p.x,p.y,g.x,g.y,g2);
        }
    }
}

/* sort point locations by given value
   v array is destroyed */
static void sortPoints(EdgePoint *ep, float *v, const int n)
{
  int m;

  /* Simple selection sort for now.  speed here not an issue */
  for (m=0;m<n-1;m++)
    {
      int k,kHi;
      float hi;
      EdgePoint tp;

      /* find next highest value */
      kHi=m;
      hi=v[kHi];
      for (k=m+1;k<n;k++)
        {
          if (v[k]>hi)
            {
              kHi=k;
              hi=v[kHi];
            }
        }

      tp = ep[m];
      ep[m]=ep[kHi];
      ep[kHi] = tp;
      v[kHi] = v[m];
    }
}

/* find ridge points in image above given thresh, and sort by edge strength */
static EdgePoint *highPoints(const ImgInfo *img, const float thresh,
                             int *nHi)
{
  EdgePoint *ep;
  float *g2;
  int y,nSkip;
  const static int HI_CHUNK = 1024;

  *nHi=0;
  ep = (EdgePoint *)myCalloc(HI_CHUNK,sizeof(EdgePoint));
  g2 = (float *)myCalloc(HI_CHUNK,sizeof(float));
  if ((ep==NULL)|(g2==NULL))
    {
      fputs("Out of memory allocating hi-point buffers\n",stderr);
      myFree(g2);
      myFree(ep);
      return(NULL);
    }

  nSkip=0;
  for(y=1;y<img->ny-1;y++)
    {
      const float *rx,*ry;
      int x;

      rx = img->gx + y*(img->nx);
      ry = img->gy + y*(img->nx);
      for(x=1;x<img->nx-1;x++)
        {
          float gx,gy,gg;

          gx=rx[x];
          gy=ry[x];
          gg=gx*gx+gy*gy;
          if (gg>thresh)
            {
              EdgePoint g,p;
              float g2pk;

              if ((*nHi+1) % HI_CHUNK == 0)
                {
                  ep = (EdgePoint *)myRealloc(ep,(*nHi+1+HI_CHUNK)*sizeof(EdgePoint));
                  g2 = (float *)myRealloc(g2,(*nHi+1+HI_CHUNK)*sizeof(float));
                  if ((ep==NULL)|(g2==NULL))
                    {
                      fputs("Out of memory re-allocating hi-point buffers\n",stderr);
                      myFree(g2);
                      myFree(ep);
                      return(NULL);
                    }
                }
              /* Do a rough ridge-ascent, just enough to guarantee
                 removal when we check previous edge lines */
              p = refineRidgeCrest(img,x,y,img->step/4,&g,&g2pk);
              if (g2pk>thresh)
                {
                  g.x = p.x-x;
                  g.y = p.y-y;
                  gx = g.x*g.x + g.y*g.y;
                  if (gx < 2.01f)
                    {
                      ep[*nHi]=p;
                      g2[*nHi]=g2pk;
                      *nHi += 1;
                    }
                  else
                    {
#ifdef _DEBUG
                      fprintf(stderr,"%d,%d ascended to %.2f,%.2f and was ignored as seed\n",x,y,p.x,p.y);
#endif
                      nSkip++;
                    }
                }
              else
                {
                  fputs("Should never get here.  g2 went down in refine\n",stderr);
                  nSkip++;
                }
            }
        }
    }
  #ifdef _DEBUG
  fprintf(stderr,"%d seeds, %d skipped\n",*nHi,nSkip);
  #endif
  /* sort by g2 value */
  sortPoints(ep,g2,*nHi);
  myFree(g2);
  return(ep);
}

/* Mark points as 'used' from given starting point,
   following gradient to a valley */
static void clearGradientToMinima(ImgInfo *img, EdgePoint p, int dir)
{
  float ds,g2,gPrev;
  EdgePoint g;
  int nx;

  ds = 0.6*dir;  /* step size */
  g2 = 9e9;
  nx = img->nx;
  imgGradient(img,p,&g,&g2);
  gPrev=2*g2;
  while(g2 < gPrev)
    {
      int x,y;

      x = (int)(p.x + 0.5);
      y = (int)(p.y + 0.5);
      img->used[y*nx+x] = (unsigned char)1;

      gPrev = g2;
      imgGradient(img,p,&g,&g2);
      p.x += g.x * ds;
      p.y += g.y * ds;
    }
}

/* Walk along a ridge, starting in the indicated direction */
static EdgePoint *traceRigde(ImgInfo *img,
       EdgePoint pk, EdgePoint dv,
       float step, float res, float thresh,
       int *np)
{
  const static int CHUNK_SIZ=1024;
  EdgePoint *ep,dvPrev,p;
  int k,lag,maxOverlap,nOverlap,nx;
  float z,g2,maxStep;

  *np=1;
  ep = (EdgePoint *)myCalloc(CHUNK_SIZ,sizeof(EdgePoint));

  k = 0;
  ep[k++]=pk;
  z = (float)sqrt((double)(dv.x*dv.x + dv.y*dv.y));
  z = step/z;
  dv.x *= z;  dv.y *= z;
  dvPrev=dv;
  p = pk;
  g2 = 2*thresh;
  lag = (int)(2.0/step+0.999); /* mark pixels as used this many samples back */
  maxOverlap = lag; /* max overlapping line points to allow */
  nOverlap=0;
  nx = img->nx;
  maxStep=2*1.4*step;
  if (maxStep < 1.5) maxStep=1.5;

  while(g2 > thresh)
    {
      float s;
      EdgePoint d,g;
      int x,y;

      p.x += dv.x;
      p.y += dv.y;

      if ((k % CHUNK_SIZ == 0) && (k > 0))
        {
          ep = (EdgePoint *)myRealloc(ep,(k+CHUNK_SIZ)*sizeof(EdgePoint));
          if (ep == NULL)
            {
              fputs("Out of memory growing buffer for edge points\n",stderr);
              return(ep);
            }
        }
      *np = k;

      p = refineRidgeCrest(img,p.x,p.y,step/4,&g,&g2);
      s = step/10;
      while (res < s)
        {
          p = refineRidgeCrest(img,p.x,p.y,s,&g,&g2);
          s /= 10;
        }
      if (res < s*10.01)
          p = refineRidgeCrest(img,p.x,p.y,res,&g,&g2);

      d.x = p.x - ep[k-1].x;
      d.y = p.y - ep[k-1].y;
      s = (float)sqrt(d.x*d.x + d.y*d.y);
      if (s > maxStep) /* skipped to a different ridge, end this line */
          return(ep);

      /* check if we are headed off of image */
      if ((p.x<=1) || (p.y<=1) || (p.x >= nx-2) || (p.y > img->ny-2))
          return(ep);

      x = (int)(p.x+0.5);
      y = (int)(p.y+0.5);
      if (img->used[y*nx+x]>0)
        {
          if (++nOverlap >= maxOverlap) /* we are re-tracing a line */
            {
              if (nOverlap >= k) *np=0; /* no new edges found */
              return(ep);
            }
        }
      else nOverlap=0;

      ep[k++] = p;

      dvPrev = dv;
      dv.x = -g.y;
      dv.y =  g.x;
      s = dv.x*dvPrev.x + dv.y*dvPrev.y;
      if (s < 0) /* make sure in same general direction */
        {
          dv.x *= -1;
          dv.y *= -1;
        }
      dv.x *= step;
      dv.y *= step;

      if (k > lag) /* mark points as used which are well in interior of edge */
        {
          clearGradientToMinima(img,ep[k-lag-1], 1);
          clearGradientToMinima(img,ep[k-lag-1],-1);
        }
    }
  return(ep);
}


static EdgePoint *followEdgeRidge(const EdgePoint start,
                                  const float step, const float res,
                                  const float thresh,
                                  ImgInfo *img, int *nLine)
{
  EdgePoint pk,g,dv,*e1,*e2,gg,*ep;
  float g2;
  int ne1,ne2;

  /* refine peak location and gradient */
  *nLine=0;
  pk = refineRidgeCrest(img,start.x,start.y,step/10,&g,&g2);
  if (g2 <= 0.0f) return(NULL);
  if (res < step/10)
    pk = refineRidgeCrest(img,pk.x,pk.y,res,&g,&g2);
  if (g2 <= 0.0f) return(NULL);

  gg.x = pk.x-start.x;
  gg.y = pk.y-start.y;
  g2 = (float)sqrt((double)(gg.x*gg.x + gg.y*gg.y));
  if (g2 > step*1.5f) /* wandered to a different edge.  ignore this peak */
    return(NULL);

  /* Follow edge in tangent to gradient direction */
  dv.x=-g.y;  dv.y=g.x;

  e1 = traceRigde(img,pk,dv,step,res,thresh,&ne1);
  dv.x *= -1;  dv.y *= -1;
  e2 = traceRigde(img,pk,dv,step,res,thresh,&ne2);

  if((ne1<=0)||(ne2<=0))
    { /* only good in one direction.  return this trace. */
      EdgePoint *e;
      
      if (ne1<=0)
        {
          if (ne2<=0)
            { /* no good lines */
              *nLine=0;
              e=NULL;
            }
          else
            {
              *nLine=ne2;
              e=e2;
            }
        }
      else
        {
          *nLine=ne1;
          e=e1;
        }
      if (*nLine>0)
        {
          int k;
          ep = (EdgePoint *)myCalloc(*nLine,sizeof(EdgePoint));
          if (ep==NULL){fputs("memory error on line buffer\n",stderr);*nLine=0;}
          for(k=0;k<*nLine;k++)ep[k]=e[k];
        }
    }
  else
    { /* combine traces into one line */
      int i,k;

      *nLine = ne1+ne2-1;
      ep = (EdgePoint *)myCalloc(*nLine,sizeof(EdgePoint));
      i=0;
      for(k=ne1-1;k>0;k--) ep[i++] = e1[k];
      for(k=0;k<ne2;k++)   ep[i++] = e2[k];
    }
  myFree(e1);
  myFree(e2);
  return(ep);
}

/** Prune points on used ridges */
static int removeTraversedSeeds(ImgInfo *img, EdgePoint *seeds,
                                const int nSeeds)
{
  int k,n,nx;
  unsigned char *u;

  nx = img->nx;
  u = img->used;
  n=0;
  for (k=1;k<nSeeds;k++)
    {
      int x,y;
      float dx,dy;

      x = (int)(seeds[k].x + 0.5f);
      y = (int)(seeds[k].y + 0.5f);
      if(u[y*nx+x]) continue;

      /* Add a little fuzz for points near pixel edges */
      dx = seeds[k].x-x;
      dy = seeds[k].y-y;
      if ((dx >  0.3) && (u[y*nx+x+1])) continue;
      if ((dx < -0.3) && (u[y*nx+x-1])) continue;
      if ((dy >  0.3) && (u[y*nx+nx+x])) continue;
      if ((dy < -0.3) && (u[y*nx-nx+x])) continue;

      /* If we got here, this string edge point must be in a
         region which is still active */
      seeds[n++] = seeds[k];
    }
  return(n);
}

/** Return index of closest point on given line to point */
static int closestPointOnLine(const EdgePoint *line, const int nLine,
			      const EdgePoint p,
			      float *dist)
{
  double d0;
  int k,k0;

  k0=0;
  d0=9e99;
  for (k=0;k<nLine;k++)
    {
      double dx,dy,d;
      dx = (double)(p.x) - (double)(line[k].x);
      dy = (double)(p.y) - (double)(line[k].y);
      d = dx*dx + dy*dy;
      if (d < d0)
	{
	  k0=k;
	  d0=d;
	}
    }
  *dist = (float)sqrt(d0);
  return(k0);
}

/** find the line segments surrounding the closest point on line
        to given point
*/
static const EdgePoint *closestLineSegment(
        const EdgePoint *line, const int nLine,
        const EdgePoint p,
	float *dist,
	int *np)
{
  int k;

  k = closestPointOnLine(line,nLine,p,dist);
  *np = (nLine < 3) ? nLine : 3;
  if (k<=0)
    { /* closest was start of line */
      return(line);
    }
  if (k>=nLine-2)
    { /* closest to end of line */
      return(line+(nLine-*np));
    }

  /* typical, mid-line case */
  return(line+(k-1));
}

/* find closest point on piecewise-linear line segments to point */
static float distanceToSegment(const EdgePoint *q,const int nq,
			       const EdgePoint p)
{
  double d2;
  int k;

  /* check endpoints first */
  d2 = 9e99;
  for (k=0;k<nq;k++)
    {
      double dx,dy,d;
      dx=(double)p.x - (double)q[k].x;
      dy=(double)p.y - (double)q[k].y;
      d = dx*dx + dy*dy;
      if (d<d2)
	d2=d;
    }

  /* check segments */
  for (k=0;k<nq-1;k++)
    {
      double px,py,x0,y0,x1,y1,d1,c;

      x0=q[k].x;
      y0=q[k].y;
      px=(double)p.x - x0;
      py=(double)p.y - y0;
      x1 = (double)q[k+1].x - x0;
      y1 = (double)q[k+1].y - y0;
      d1 = sqrt(x1*x1 + y1*y1);
      x1/=d1;y1/=d1;
      c = px*x1 + py*y1;  /* project p onto line */
      if ((c>0) && (c<d1))
	{ /* projection IS on segment */
	  x1*=c; y1*=c;  /* move x1,y1 to closest point on segment */
	  d1 = x1*x1 + y1*y1;
          if (d1<d2)
	    d2=d1;  /* new closest point */
	}
    }
  return((float)sqrt(d2));
}

/** check for any remaining seeds near line */
static int removeSeedsNearLine(const EdgePoint *line, const int nLine,
			       const float step,
			       EdgePoint *seeds, int *nSeeds)
{
  int m,k;
  float smallStep, neigh;

  smallStep = 0.5f * step;
  neigh     = 5.0f * step;

  m=0; /* No. of retained seeds */
  for (k=0;k<*nSeeds;k++)
    {
      const EdgePoint *p;
      float d;
      int np;

      p = closestLineSegment(line,nLine,seeds[k],&d,&np);
      if (d < smallStep)
	continue;  /* close enough.  discard this seed */
      if (d < neigh)
	{ /* pretty close.  take a closer look */
	  d = distanceToSegment(p,np,seeds[k]);
	}
      if (d < smallStep)
	continue;  /* close enough.  discard this seed */

      /* not too close.  keep it around */
      seeds[m++]=seeds[k];
    }
  *nSeeds=m;
  return(0);
}

/*======================================================================
    Get list of points forming edge lines in an image
*/
EdgeLines subCannyEdge(
     const float *img0, /* Image buffer, raster order */
     const int nx, /* Number of pixels per scan */
     const int ny, /* Number of raster scans */
     const float step, /* Target resolution, in pixels, for edge lines */
     const float res, /* Resolution of edge points (pixels) */
     const float hiThresh, /* Minimum normalized strength of gradient^2 to be used for a seed in an edge line */
     const float loThresh) /* Minimem normalized strength of gradient^2 to be considered part of an edge line */
{
  EdgeLines el;
  ImgInfo img;
  EdgePoint *seeds;
  int nSeeds,i;

  initEdgeLines(&el);
  initImgInfo(img0,nx,ny,&img);
  img.step=step;
  img.res=res;

  seeds = highPoints(&img,hiThresh,&nSeeds);

  i=0;
  while(nSeeds>0)
    {
      EdgePoint *line;
      int nLine;

      /* Trace ridge associated with next strongest edge seed */
      line = followEdgeRidge(seeds[0],step,res,loThresh,&img,&nLine);
      #ifdef _DEBUG
      fprintf(stderr,"Line %d, %d points (%d seeds)\n",i+1,nLine,nSeeds);fflush(stderr);
      #endif
      if (nLine > 0)
        {
          growEdgeLines(&el);
          el.line[i] = line;
          el.lineLen[i] = nLine;
          el.nLines=++i;
        }
      #ifdef _DEBUG
      else fputs("Skipping NULL line\n",stderr);
      #endif

      /* prune line seeds which are already used in tracing lines */
      nSeeds = removeTraversedSeeds(&img,seeds,nSeeds);

      /* check for any remaining seeds near this new line */
      removeSeedsNearLine(line,nLine,step,seeds,&nSeeds);
    }

  /* check edge line ends for overlap and prune here (TBD) */
#ifdef _DEBUG_VERBOSE
  fputs("seed buffer: ",stderr);fflush(stderr);
#endif
  myFree(seeds);
  freeImgInfo(&img);
  return(el);
}

/***********************************************************************
2015-06-09   ab   octave port
2008-02-06   ab   prune seeds points as we go
...
2007/10/23   ab
*/
