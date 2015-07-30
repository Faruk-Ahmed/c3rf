#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"
#include <time.h>

static int cmpfnc (const void * a, const void * b)
{
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variables */
    int ROWS, COLS, LLen0, LLen1, i, j, n;
    double *spids, *uniqueSpids, *L, *splabels, *seg, numUniques;
    
    /* Input */
    L           = mxGetPr(prhs[0]);
    splabels    = mxGetPr(prhs[1]);
    uniqueSpids = mxGetPr(prhs[2]);

    /* Compute Sizes */
    ROWS       = mxGetDimensions(prhs[1])[0];
    COLS       = mxGetDimensions(prhs[1])[1];
    numUniques = mxGetDimensions(prhs[2])[0];

    /* Output */
    plhs[0] = mxCreateDoubleMatrix(ROWS, COLS, mxREAL);
    seg     = mxGetPr(plhs[0]);

    for (i = 0; i < ROWS; i++)
    {
        for (j = 0; j < COLS; j++)
	{
	    seg[i+ROWS*j] = 0;
	    for (n = 0; n < numUniques; n++)
	    {
 	        if (splabels[i+ROWS*j] == uniqueSpids[n])
		{
		    seg[i+ROWS*j] = L[n];
                    break;
		}
            }
	}
    }
}

