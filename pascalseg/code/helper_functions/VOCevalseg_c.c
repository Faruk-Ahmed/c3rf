/* A C version of the PASCAL VOCdevkit version of the
 * semantic segmentation evaluation in MATLAB.
 * Being in C, this runs faster.
 * (Faruk Ahmed, 2014)
 */

#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variables */
    int i, j, k, nclasses, num, npredsegs, ngtsegs, count, ROWSGT, ROWSRESIM, COLSGT, COLSRESIM, *sumim, *hs;    
    double *accuracies, *confcounts, *gtim, *resim, gtj, resj, gtjresj;
    
    /* Input */
    nclasses = mxGetPr(prhs[0])[0];

    /* Sizes */
    npredsegs = mxGetDimensions(prhs[1])[1];
    ngtsegs   = mxGetDimensions(prhs[2])[1];
    
    if (npredsegs != ngtsegs) mexErrMsgTxt("The number of predsegs and gtsegs do not match.");
    
    /* Output */
    plhs[0] = mxCreateDoubleMatrix(nclasses+1,1,mxREAL);
    
    /* Set up things */
    accuracies = mxGetPr(plhs[0]);
    count = 0;
    num = nclasses + 1;
    
    confcounts = mxCalloc(num * num, sizeof(double));
    
    for (i = 0; i < npredsegs; i++)
    {
        resim = mxGetPr(mxGetCell(prhs[1], i));
        gtim  = mxGetPr(mxGetCell(prhs[2], i));
        
        /* Get sizes of image */
        ROWSRESIM = mxGetDimensions(mxGetCell(prhs[1], i))[0];
        COLSRESIM = mxGetDimensions(mxGetCell(prhs[1], i))[1];

        ROWSGT = mxGetDimensions(mxGetCell(prhs[2], i))[0];
        COLSGT = mxGetDimensions(mxGetCell(prhs[2], i))[1];

        /* Sizes of gtim and resim should match. */
        if ((ROWSRESIM != ROWSGT) || (COLSRESIM != COLSGT))
            mexErrMsgTxt("Ground truth dimensions do not match predicted image dimensions. What are you doing?");
        
        sumim      = mxMalloc(ROWSGT * COLSGT * sizeof(int));
        hs         = mxCalloc(num * num, sizeof(int));

        /* Construct joint histogram: */
        for (j = 0; j < ROWSGT; j++)
        {
            for (k = 0; k < COLSGT; k++)
            {
                sumim[j*COLSGT+k] = 1 + gtim[j*COLSGT+k] + resim[j*COLSGT+k]*num;
                if (resim[j*COLSGT+k] > nclasses) mexErrMsgTxt("The results image has an out-of-range value.");
            }
        }

        for (j = 0; j < ROWSGT; j++)
        {
            for (k = 0; k < COLSGT; k++)
            {
                if (gtim[j*COLSGT+k] != 255)
                {
                    hs[sumim[j*COLSGT+k]-1] = hs[sumim[j*COLSGT+k]-1] + 1;
                    count += 1;
                }
            }
        }

        for (j = 0; j < num*num; j++)
            confcounts[j] += hs[j];
        
        mxFree(sumim);
        mxFree(hs);
    }

    for (j = 0; j < num; j++)
    {
        gtj  = 0;
        resj = 0;
        for (k = 0; k < num; k++)
        {
            gtj += confcounts[j*num+k];
            resj += confcounts[k*num+j];
        }
        gtjresj = confcounts[j*num+j];
        accuracies[j] = 100*gtjresj/(gtj + resj - gtjresj);
    }   

    mxFree(confcounts);
    return;
}
