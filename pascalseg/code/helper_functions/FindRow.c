#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *S, *s;
  mwSize n1, n2, i1, i2, k;
    S  = mxGetPr(prhs[0]);
    n1 = mxGetM(prhs[0]);
    n2 = mxGetN(prhs[0]);
    s  = mxGetPr(prhs[1]);  
    for (i1 = 0; i1 < n1; i1++) 
    {
       k = i1;
       for (i2 = 0; i2 < n2; i2++) 
       {
          if (S[k] != s[i2]) 
          {
             break;
          }
          k += n1;
       }
       if (i2 == n2) {
          plhs[0] = mxCreateDoubleScalar((double) (i1 + 1));
          return;
       }
    }
    
    plhs[0] = mxCreateDoubleScalar(mxGetNaN());
}
