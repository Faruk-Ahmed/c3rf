#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

static int cmpfnc (const void * a, const void * b)
{
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variables */
    int ROWS, COLS, LLen0, LLen1, i, j, n, numUniques;
    double *spids, *uniqueSpids, *splabels, *numUniquesTotal;

    /* Input */
    splabels = mxGetPr(prhs[0]);

    /* Compute Sizes */
    ROWS    = mxGetDimensions(prhs[0])[0];
    COLS    = mxGetDimensions(prhs[0])[1];

    /* Output */
    plhs[0]          = mxCreateDoubleMatrix(ROWS*COLS, 1, mxREAL);
    plhs[1]          = mxCreateDoubleMatrix(1, 1, mxREAL);
    uniqueSpids      = mxGetPr(plhs[0]);
    numUniquesTotal  = mxGetPr(plhs[1]);

    /* Get unique spids: */
    spids = (double *)malloc(ROWS*COLS*sizeof(double));

    for (i = 0; i < ROWS; i++)
        for (j = 0; j < COLS; j++)
	    spids[i+ROWS*j] = splabels[i+ROWS*j];
    qsort(spids, ROWS*COLS, sizeof(spids[0]), cmpfnc);

    numUniques = 0;
    for (i = 0; i < ROWS*COLS; i++)
    {
        if (i == 0)   
	{
	    uniqueSpids[numUniques] = spids[i];
	    numUniques += 1;
	}
	else
	{
	    if (spids[i] != uniqueSpids[numUniques-1])
	    {
	        uniqueSpids[numUniques] = spids[i];
		numUniques += 1;
            }
	}
    }
    numUniquesTotal[0] = numUniques;
    free(spids);
}
