/*

MEX front end for matlab to call subCannyEdge

Can be compiled with:
mex subCannyEdgeMEX.c subCannyEdge.c

For optimization, I recommend compiling the main code by hand

gcc -O -fPIC -Wall -c subCannyEdge.c
mex subCannyEdgeMEX.c subCannyEdge.o

Detailed options to try:
     i386 and x86-64 Options:
         -mtune=x86-64 -march=x86-64 -mmmx  -msse  -msse2 -msse3 -m64 
     General opt:
       -fprefetch-loop-arrays -ffinite-math-only -fno-trapping-math
       -foptimize-reg
       -ffast-math 
       -floop-optimize -finline-functions -finline-functions-called-once
       -fkeep-inline-functions -fkeep-static-consts
       -fmerge-constants  -fmerge-all-constants
       -floop-optimize2 -fmove-loop-invariants
       -fno-math-errno
       -fno-peephole2 -funsafe-math-optimizations
       -frerun-loop-opt -frounding-math
       -fsingle-precision-constant
       -funroll-loops

For debug:
mex -g -D_DEBUG subCannyEdgeMEX.c subCannyEdge.c

*/
#include "subCannyEdge.h"
#include "mex.h"

static float *parseInputs(const int nrhs, const mxArray *prhs[],
       int *nx, int *ny,
       float *step, float *res,
       float *hi, float *lo)
{
  float *img;
  if (nrhs != 5)
    {
      puts("incorrect no. of inputs : subCannyEdgeMEX(single(img),step,res,hi,lo)");
      //mxErrMsgTxt("incorrect no. of inputs : subCannyEdgeMEX(single(img),step,res,hi,lo)");
      return(NULL);
    }

  img = (float *)mxGetPr(prhs[0]);
  *ny   =  mxGetN(prhs[0]);  /* cols */
  *nx   =  mxGetM(prhs[0]);  /* rows, contigous in matlab */
  *step = (float)mxGetScalar(prhs[1]);
  *res  = (float)mxGetScalar(prhs[2]);
  *hi   = (float)mxGetScalar(prhs[3]);
  *lo   = (float)mxGetScalar(prhs[4]);
  return(img);
}

static mxArray  *unpackEdgeLines(const EdgeLines *el)
{
  const static char *fieldName="p";
  mxArray *es;
  int k;

  es = mxCreateStructMatrix(1,el->nLines,1,&fieldName);
  for (k=0; k< el->nLines; k++)
    {
      mxArray *v;
      double *p;
      int m;

      v = mxCreateNumericMatrix(2,el->lineLen[k],mxDOUBLE_CLASS,mxREAL);
      p = (double *)mxGetPr(v);

      for (m=0;m<el->lineLen[k];m++)
        {
          EdgePoint ep;

          ep = el->line[k][m];
          p[2*m  ] = ep.x+1;  /* make match matlab conventions */
          p[2*m+1] = ep.y+1;
        }
      mxSetField(es,k,fieldName,v);
    }
  return(es);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* inputs */
  float *img,step,res,hiThresh,loThresh;
  int nx,ny;

  /* output */
  double *pStat,dummyStatus;
  /*mxArray *es;   Structure of edge-line arrays */

  /* Local */
  EdgeLines el;

  if (nlhs>1)
    { /* enable status return */
      plhs[1] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
      pStat = (double *)mxGetPr(plhs[1]);
    }
  else pStat = &dummyStatus;

  img = parseInputs(nrhs,prhs,&nx,&ny,&step,&res,&hiThresh,&loThresh);

  if (img)
    {
      el = subCannyEdge(img,nx,ny,step,res,hiThresh,loThresh);
      plhs[0] = unpackEdgeLines(&el);
      freeEdgeLines(&el);
    }
  else
    { plhs[0] = NULL; }
}

/***********************************************************************
2008-02-06   ab  transposed dims
2007/11/02   ab  removed image size params.  get them from image matrix
2007/10/23   ab
*/
