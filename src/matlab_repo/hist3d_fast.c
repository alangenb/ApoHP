/*    hist3d_fast.c

      Mike Lawrence 2012-11-17

      2020-07-26 added optional last argument "output_type" which can be:
        1   double (default)
        2   uint8 (no bin's counts go over 255)
        3   uint16 (no bin's counts go over 65535)
        4   uint32 (no bin's counts go over 4.2 billion)

      % to compile: (from Matlab prompt)
        mex hist3d_fast.c

*/

#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
  long xlen, ylen, zlen, len;
  long firstx, lastx, firsty, lasty, firstz, lastz, nxvals, nyvals, nzvals;
  int dims[3];
  double *x, *y, *z, *xvals, *yvals, *zvals;
  long i, xi, yi, zi, idx;

  int output_type = 1;
  double *h_double;
  unsigned char *h_uint8;
  unsigned short int *h_uint16;
  unsigned int *h_uint32;
  long n_overflow = 0;
  
  if (nrhs!=3 && nrhs!=4 && nrhs!=9 && nrhs!=10) mexErrMsgTxt("wrong number of inputs, was expecting:  x, y [, firstx, lastx, firsty, lasty, firstz, lastz] [, output_type]");
  if (nlhs>1) mexErrMsgTxt("too many outputs requested");

  xlen = mxGetN(prhs[0]); if (xlen==1) xlen = mxGetM(prhs[0]);
  ylen = mxGetN(prhs[1]); if (ylen==1) ylen = mxGetM(prhs[1]);
  zlen = mxGetN(prhs[2]); if (zlen==1) zlen = mxGetM(prhs[2]);
  if (xlen!=ylen || xlen!=zlen) mexErrMsgTxt("error: x,y,z are of different lengths");
  len = xlen;
  x = mxGetPr(prhs[0]);
  y = mxGetPr(prhs[1]);
  z = mxGetPr(prhs[2]);

  if (nrhs>=9) {
    if (nrhs==10) output_type = (long)mxGetScalar(prhs[9]);
    /* user supplied ranges */
    firstx = (long)mxGetScalar(prhs[3]);
    lastx = (long)mxGetScalar(prhs[4]);
    firsty = (long)mxGetScalar(prhs[5]);
    lasty = (long)mxGetScalar(prhs[6]);
    firstz = (long)mxGetScalar(prhs[7]);
    lastz = (long)mxGetScalar(prhs[8]);
  } else {
    if (nrhs==4) output_type = (long)mxGetScalar(prhs[3]);
    /* compute ranges as min and max of x,y,z */
    firstx = *x;
    lastx = *x;
    firsty = *y;
    lasty = *y;
    firstz = *z;
    lastz = *z;
    for (i=1;i<len;i++) {
      if ((*(x+i))<firstx) firstx = (*(x+i));
      if ((*(x+i))>lastx) lastx = (*(x+i));
      if ((*(y+i))<firsty) firsty = (*(y+i));
      if ((*(y+i))>lasty) lasty = (*(y+i));
      if ((*(z+i))<firstz) firstz = (*(z+i));
      if ((*(z+i))>lastz) lastz = (*(z+i));
    }
    printf("x range = %d to %d;  y range = %d to %d;  z range = %d to %d\n",firstx,lastx,firsty,lasty,firstz,lastz);
  }

  nxvals = lastx-firstx+1;
  nyvals = lasty-firsty+1;
  nzvals = lastz-firstz+1;
  if (nxvals<1 || nyvals<1 || nzvals<1) mexErrMsgTxt("cannot have firstx>lastx, firsty>lasty, or firstz>lastz");

  dims[0] = nxvals;
  dims[1] = nyvals;
  dims[2] = nzvals;  

  /*    x = rows        */
  /*    y = columns     */
  /*    z = pages       */

  switch(output_type) {
  case 1:
    plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    h_double = mxGetPr(plhs[0]);
    for(i=0;i<len;i++) {
      xi = (long)(*(x+i))-firstx;
      yi = (long)(*(y+i))-firsty;
      zi = (long)(*(z+i))-firstz;
      if (xi<0 || xi>=nxvals || yi<0 || yi>=nyvals || zi<0 || zi>=nzvals) continue;  /* ignore out-of-range values */
      idx = (nyvals*nxvals*zi)+(nxvals*yi)+xi;
      /* no bin overflow check needed for double output */
      (*(h_double+idx))++;
    }
    break;
  case 2:
    plhs[0] = mxCreateNumericArray(3,dims,mxUINT8_CLASS,mxREAL);
    h_uint8 = (unsigned char *)mxGetData(plhs[0]);
    for(i=0;i<len;i++) {
      xi = (long)(*(x+i))-firstx;
      yi = (long)(*(y+i))-firsty;
      zi = (long)(*(z+i))-firstz;
      if (xi<0 || xi>=nxvals || yi<0 || yi>=nyvals || zi<0 || zi>=nzvals) continue;  /* ignore out-of-range values */
      idx = (nyvals*nxvals*zi)+(nxvals*yi)+xi;
      if (*(h_uint8+idx)==255) n_overflow++;
      else (*(h_uint8+idx))++;
    }
    break;
  case 3:
    plhs[0] = mxCreateNumericArray(3,dims,mxUINT16_CLASS,mxREAL);
    h_uint16 = (unsigned short int *)mxGetData(plhs[0]);
    for(i=0;i<len;i++) {
      xi = (long)(*(x+i))-firstx;
      yi = (long)(*(y+i))-firsty;
      zi = (long)(*(z+i))-firstz;
      if (xi<0 || xi>=nxvals || yi<0 || yi>=nyvals || zi<0 || zi>=nzvals) continue;  /* ignore out-of-range values */
      idx = (nyvals*nxvals*zi)+(nxvals*yi)+xi;
      if (*(h_uint16+idx)==65535) n_overflow++;
      else (*(h_uint16+idx))++;
    }
    break;
  case 4:
    plhs[0] = mxCreateNumericArray(3,dims,mxUINT32_CLASS,mxREAL);
    h_uint32 = (unsigned int *)mxGetData(plhs[0]);
    for(i=0;i<len;i++) {
      xi = (long)(*(x+i))-firstx;
      yi = (long)(*(y+i))-firsty;
      zi = (long)(*(z+i))-firstz;
      if (xi<0 || xi>=nxvals || yi<0 || yi>=nyvals || zi<0 || zi>=nzvals) continue;  /* ignore out-of-range values */
      idx = (nyvals*nxvals*zi)+(nxvals*yi)+xi;
      if (*(h_uint32+idx)==4294967295) n_overflow++;
      else (*(h_uint32+idx))++;
    }
    break;
  default: mexErrMsgTxt("invalid output_type");
  }

  if (n_overflow) printf("WARNING: bin overflow led to loss of %ld item(s).\n",n_overflow);
  
    
}


