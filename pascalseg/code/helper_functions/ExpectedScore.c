#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variables */
    int i, j, k, ROWS, COLS; 
    double *predseg, *expectedscore, dummy, *marginalsForK, sumExpectedUnion, sumExpectedIntersect, nlabels;

    /* Input */
    predseg       = mxGetPr(prhs[0]);
    nlabels       = mxGetPr(prhs[2])[0];
    
    /* Compute Sizes */
    ROWS = mxGetDimensions(prhs[0])[0];
    COLS = mxGetDimensions(prhs[0])[1];
    
    /* Output */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    /* Set up things */
    expectedscore = mxGetPr(plhs[0]);
    dummy = 0;
    
    /* Do what you have to do */    
    for (i = 0; i < nlabels; i++)
    {
        marginalsForK         = mxGetPr(mxGetCell(prhs[1], i));
        sumExpectedUnion      = 0.0;
        sumExpectedIntersect  = 0.0;

        for (j = 0; j < ROWS; j++)
        {
            for (k = 0; k < COLS; k++)
            {
                if (predseg[j*COLS+k] == i)
                {
                    sumExpectedUnion     = sumExpectedUnion + 1;
                    sumExpectedIntersect = sumExpectedIntersect + marginalsForK[j*COLS+k];
                }
                else
                {
                    sumExpectedUnion     = sumExpectedUnion + marginalsForK[j*COLS+k];
                }
            }
        }

        dummy = dummy + (1/nlabels)*sumExpectedIntersect/sumExpectedUnion;
    }
    expectedscore[0] = dummy;
}

