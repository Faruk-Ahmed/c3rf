/* Loopy Belief Propagation in log-space.
 * This is strongly based on Mark Schmidt's version of Loopy BP available here at
 * http://www.cs.ubc.ca/~schmidtm/Software/UGM.html
 * (Faruk Ahmed, 2014)
 */

#include <math.h>
#include "mex.h"
#include "UGM_common.h"

#define EPSILON log(1e-16)
#define MINUSINF -100000

double myLog(double);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Variables */
    int n, s, e, e2, n1, n2, neigh, Vind, Vind2, s1, s2, di, L, foundIt, 
        nNodes, nEdges, maxState, dims[3],
        iter,maxIter,nNbrs,
        *edgeEnds, *nStates, *V, *E,*y,
        calcMargs, useDamping, *nNeighbors;
    
    double *msgFromNode, *edgePot, *nodeBel, *edgeBel, *logZ,
           z,energy1,energy2,entropy1,entropy2,*prodMsgs,*oldMsgs,*newMsgs,*tmp,
           *updatedMsg, *nodeBelUptree, *msgFromNodeForUptree, *prodMsgForUptree, *prodMsgForFactorMargs, z1, z2, z3, *nodePot, *msgDownTree, *factorMarginals,
           nodeTerm, factorTerm, zf, *messagesFromLastLBP, alpha, M1, temp;
    
   /* Input */
    msgFromNode = mxGetPr(prhs[0]);
    msgFromNodeForUptree = mxGetPr(prhs[1]);
    nodePot = mxGetPr(prhs[2]);
    msgDownTree = mxGetPr(prhs[3]);
    nNeighbors = (int*)mxGetPr(prhs[4]);
    edgePot = mxGetPr(prhs[5]);
    edgeEnds = (int*)mxGetPr(prhs[6]);
    nStates = (int*)mxGetPr(prhs[7]);
    V = (int*)mxGetPr(prhs[8]);
    E = (int*)mxGetPr(prhs[9]);
    maxIter = ((int*)mxGetPr(prhs[10]))[0];
    calcMargs = ((int*)mxGetPr(prhs[11]))[0];
    useDamping = ((int*)mxGetPr(prhs[12]))[0];
    messagesFromLastLBP = mxGetPr(prhs[13]);
    alpha = ((double*)mxGetPr(prhs[14]))[0];
   
    if (!mxIsClass(prhs[6],"int32")||!mxIsClass(prhs[7],"int32")||!mxIsClass(prhs[8],"int32")||!mxIsClass(prhs[9],"int32")||!mxIsClass(prhs[10],"int32"))
		mexErrMsgTxt("edgeEnds, nStates, V, E, maxIter must be int32");
	
    /* Compute Sizes */

    nNodes = mxGetDimensions(prhs[0])[0];
    maxState = mxGetDimensions(prhs[0])[1];
    nEdges = mxGetDimensions(prhs[6])[0];
    
    /* Output */
    plhs[0] = mxCreateDoubleMatrix(nNodes,maxState,mxREAL);
    dims[0] = maxState;
    dims[1] = maxState;
    dims[2] = nEdges;
    plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nEdges*maxState*2,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nNodes, maxState, mxREAL);
    
    nodeBel = mxGetPr(plhs[0]);
    edgeBel = mxGetPr(plhs[1]);
    logZ = mxGetPr(plhs[2]);
    updatedMsg = mxGetPr(plhs[3]);
    nodeBelUptree = mxGetPr(plhs[4]);
    
    prodMsgs = mxCalloc(maxState*nNodes,sizeof(double));
    prodMsgForUptree = mxCalloc(maxState*nNodes,sizeof(double));
    prodMsgForFactorMargs = mxCalloc(maxState*nNodes,sizeof(double)); 
    oldMsgs = mxCalloc(maxState*nEdges*2,sizeof(double));
    newMsgs = mxCalloc(maxState*nEdges*2,sizeof(double));
    tmp = mxCalloc(maxState,sizeof(double));
    factorMarginals = mxCalloc(maxState*nNodes,sizeof(double));
    
    /* Initialize */
    for(e = 0; e < nEdges; e++)
    {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        for(s = 0; s < nStates[n2]; s++)
            newMsgs[s+maxState*e] = myLog(1./nStates[n2]);
        for(s = 0; s < nStates[n1]; s++)
            newMsgs[s+maxState*(e+nEdges)] = myLog(1./nStates[n1]);
    }
    
    for(iter = 0; iter < maxIter; iter++)
    {        
        for(n = 0; n < nNodes; n++)
        {            
            /* Update Messages */
            for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
            {
                e = E[Vind]-1;
                n1 = edgeEnds[e]-1;
                n2 = edgeEnds[e+nEdges]-1;
                
                /* First part of message is msgFromNode*/
                for(s = 0; s < nStates[n]; s++)
                    tmp[s] = msgFromNode[n + nNodes*s];
                
                /* Multiply by messages from neighbors except j */
                for(Vind2 = V[n]-1; Vind2 < V[n+1]-1; Vind2++)
                {
                    e2 = E[Vind2]-1;
                    if (e != e2)
                    {
                        if (n == edgeEnds[e2+nEdges]-1)
                        {
                            for(s = 0; s < nStates[n]; s++)
                            {
                                tmp[s] += newMsgs[s+maxState*e2];
                            }
                        }
                        else
                        {
                            for(s = 0; s < nStates[n]; s++)
                            {
                                tmp[s] += newMsgs[s+maxState*(e2+nEdges)];
                            }
                        }
                    }
                }
               
                /* Now multiply by edge potential to get new message */
                
                if (n == n2)
                {
                    for(s1 = 0; s1 < nStates[n1]; s1++)
                    {
                        temp = 0.0;
                        M1 = tmp[0] + edgePot[s1+maxState*(0+maxState*e)];

                        for(s2 = 0; s2 < nStates[n2]; s2++)
                        {
                            if ((tmp[s2] + edgePot[s1+maxState*(s2+maxState*e)]) > M1)
                                M1 = tmp[s2] + edgePot[s1+maxState*(s2+maxState*e)];
                        }
                        
                        for(s2 = 0; s2 < nStates[n2]; s2++)
                        {
                            temp += exp(tmp[s2] + edgePot[s1+maxState*(s2+maxState*e)] - M1);
                        }
                        newMsgs[s1+maxState*(e+nEdges)] = M1 + myLog(temp);
                    }
                    
                    /* Normalize */
                    M1 = newMsgs[0+maxState*(e+nEdges)];
                    for(s = 0; s < nStates[n1]; s++)
                        if (newMsgs[s+maxState*(e+nEdges)] > M1)
                            M1 = newMsgs[s+maxState*(e+nEdges)];
                    temp = 0.0;

                    z = M1;
                    for(s = 0; s < nStates[n1]; s++)
                    {
                        temp += exp(newMsgs[s+maxState*(e+nEdges)] - M1);
                    }
                    z += myLog(temp);
                        
                    for(s = 0; s < nStates[n1]; s++)
                    {
                        newMsgs[s+maxState*(e+nEdges)] -= z;
                    }
                        
                    if(useDamping == 1)
                    {
                        for(s = 0; s < nStates[n1]; s++)
                        {
                            newMsgs[s+maxState*(e+nEdges)] = myLog(alpha*exp(newMsgs[s+maxState*(e+nEdges)]) + (1-alpha)*exp(messagesFromLastLBP[s+maxState*(e+nEdges)]));
                        }
                    }
                }
                else
                {
                    for(s2 = 0; s2 < nStates[n2]; s2++)
                    {
                        temp = 0.0;
                        M1 = tmp[0] + edgePot[0+maxState*(s2+maxState*e)];
                        for(s1 = 0; s1 < nStates[n1]; s1++)
                            if ((tmp[s1] + edgePot[s1+maxState*(s2+maxState*e)]) > M1)
                                M1 = tmp[s1] + edgePot[s1+maxState*(s2+maxState*e)];

                        for(s1 = 0; s1 < nStates[n1]; s1++)
                            temp += exp(tmp[s1] + edgePot[s1+maxState*(s2+maxState*e)] - M1);
                        newMsgs[s2+maxState*e] = M1 + myLog(temp);
                    }

                    /* Normalize */
                    M1 = newMsgs[0+maxState*e];
                    for (s = 0; s < nStates[n2]; s++)
                        if (newMsgs[s+maxState*e] > M1)
                            M1 = newMsgs[s+maxState*e];

                    z = M1;
                    temp = 0.0;
                    for (s = 0; s < nStates[n2]; s++)
                        temp += exp(newMsgs[s+maxState*e] - M1);
                    z += myLog(temp);
                    
                    for (s = 0; s < nStates[n2]; s++)
                        newMsgs[s+maxState*e] -= z;
                        
                    if(useDamping == 1)
                    {
                        for(s = 0; s < nStates[n2]; s++)
                        {
                            newMsgs[s+maxState*e] = myLog(alpha*exp(newMsgs[s+maxState*e]) + (1-alpha)*exp(messagesFromLastLBP[s+maxState*e]));
                        }
                    }
                }   
            }            
        }
        

        /* oldMsgs = newMsgs */
        z = 0;
        for(s=0;s<maxState;s++)
        {
            for(e=0;e<nEdges*2;e++)
            {
                z += absDif(newMsgs[s+maxState*e],oldMsgs[s+maxState*e]);
                oldMsgs[s+maxState*e] = newMsgs[s+maxState*e];
                updatedMsg[s+maxState*e] = newMsgs[s+maxState*e];
            }
        }
        
        /* if sum(abs(newMsgs(:)-oldMsgs(:))) < 1e-4; break; */
        if(z < 1e-3)
        {
            break;
        }   
    }
        
    /* compute nodeBel */
    for(n = 0; n < nNodes; n++)
    {
        for(s = 0; s < nStates[n]; s++)
        {
            prodMsgs[s+maxState*n] = msgFromNode[n+nNodes*s];
            prodMsgForUptree[s+maxState*n] = msgFromNodeForUptree[n+nNodes*s];
            prodMsgForFactorMargs[s+maxState*n] = msgDownTree[n+nNodes*s];
        }
        
        for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
        {
            e = E[Vind]-1;
            n1 = edgeEnds[e]-1;
            n2 = edgeEnds[e+nEdges]-1;
            
            if (n == n2)
            {
                for(s = 0; s < nStates[n]; s++)
                {
                    prodMsgs[s+maxState*n] += newMsgs[s+maxState*e];
                    prodMsgForUptree[s+maxState*n] += newMsgs[s+maxState*e];
                    prodMsgForFactorMargs[s+maxState*n] += newMsgs[s+maxState*e];
                }
            }
            else
            {
                for(s = 0; s < nStates[n]; s++)
                {
                    prodMsgs[s+maxState*n] += newMsgs[s+maxState*(e+nEdges)];
                    prodMsgForUptree[s+maxState*n] += newMsgs[s+maxState*(e+nEdges)];
                    prodMsgForFactorMargs[s+maxState*n] += newMsgs[s+maxState*(e+nEdges)];
                }
            }
        }
        
        /* normalize */
        M1 = prodMsgs[0+maxState*n];
        for (s = 0; s < nStates[n]; s++)
        {
            nodeBel[n+nNodes*s] = prodMsgs[s+maxState*n];
            if (prodMsgs[s+maxState*n] > M1)
                M1 = prodMsgs[s+maxState*n];
        }
        z = M1;
        temp = 0.0;
        for (s = 0; s < nStates[n]; s++)
            temp += exp(nodeBel[n+nNodes*s] - M1);
        z += myLog(temp);
        for (s = 0; s < nStates[n]; s++)
            nodeBel[n+nNodes*s] -= z;
            
        M1 = prodMsgForUptree[0+maxState*n];
        for (s = 0; s < nStates[n]; s++)
        {
            nodeBelUptree[n+nNodes*s] = prodMsgForUptree[s+maxState*n];
            if (prodMsgForUptree[s+maxState*n] > M1)
                M1 = prodMsgForUptree[s+maxState*n];
        }
        z = M1;
        temp = 0.0;
        for (s = 0; s < nStates[n]; s++)
            temp += exp(nodeBelUptree[n+nNodes*s] - M1);
        z += myLog(temp);
        for (s = 0; s < nStates[n]; s++)
            nodeBelUptree[n+nNodes*s] -= z;
                
        M1 = prodMsgForFactorMargs[0+maxState*n];
        for (s = 0; s < nStates[n]; s++)
        {
            if (prodMsgForFactorMargs[s+maxState*n] > M1)
                M1 = prodMsgForFactorMargs[s+maxState*n];
        }
        z = M1;
        temp = 0.0;
        for (s = 0; s < nStates[n]; s++)
            temp += exp(prodMsgForFactorMargs[s+maxState*n] - M1);
        z += myLog(temp);
        for (s = 0; s < nStates[n]; s++)
        {
            prodMsgForFactorMargs[s+maxState*n] -= z;
        }
        
    }        
    

    /* Compute edgeBel */
    for(e = 0; e < nEdges; e++)
    {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        for(s1 = 0; s1 < nStates[n1]; s1++)
        {
            for(s2 = 0; s2 < nStates[n2]; s2++)
            {
                edgeBel[s1+maxState*(s2+maxState*e)] = nodeBel[n1+nNodes*s1] - newMsgs[s1+maxState*(e+nEdges)];
                edgeBel[s1+maxState*(s2+maxState*e)] += nodeBel[n2+nNodes*s2] - newMsgs[s2+maxState*e];
                edgeBel[s1+maxState*(s2+maxState*e)] += edgePot[s1+maxState*(s2+maxState*e)];
            }
        }
       
        M1 = edgeBel[maxState*(maxState*e)];
        for(s1 = 0; s1 < nStates[n1]; s1++)
            for(s2 = 0; s2 < nStates[n2]; s2++)
                if (edgeBel[s1+maxState*(s2+maxState*e)] > M1)
                    M1 = edgeBel[s1+maxState*(s2+maxState*e)];

        z = M1;
        temp = 0.0;
        for(s1 = 0; s1 < nStates[n1]; s1++)
            for(s2 = 0; s2 < nStates[n2]; s2++)
                    temp += exp(edgeBel[s1+maxState*(s2+maxState*e)] - M1);
        z += myLog(temp);
        
        for(s1 = 0; s1 < nStates[n1]; s1++)
            for(s2 = 0; s2 < nStates[n2]; s2++)
                edgeBel[s1+maxState*(s2+maxState*e)] -= z;
    }



    /* Compute Bethe Free Energy */  
    nodeTerm = 0.0;
    factorTerm = 0.0;
    
    for(n = 0; n < nNodes; n++)
    {
        nNbrs = V[n+1]-V[n];
        
       M1 = nodePot[n] + prodMsgForFactorMargs[maxState*n];
        for(s = 0; s < nStates[n]; s++)
        {
            nodeTerm += (nNbrs+nNeighbors[n]+1)*exp(nodeBel[n+nNodes*s])*nodeBel[n+nNodes*s];
            factorMarginals[s+maxState*n] = nodePot[n+nNodes*s] + prodMsgForFactorMargs[s+maxState*n];
            
            if (factorMarginals[s+maxState*n] > M1)
                M1 = factorMarginals[s+maxState*n];

        }
        zf = M1;
        temp = 0.0;
        for (s = 0; s < nStates[n]; s++)
            temp += exp(factorMarginals[s+maxState*n] - M1);
        zf += myLog(temp);
        

        for(s = 0; s < nStates[n]; s++)

        {
            factorMarginals[s+maxState*n] -= zf;
            temp = exp(factorMarginals[s+maxState*n]) * (-nodePot[n+nNodes*s] + factorMarginals[s+maxState*n]);

            factorTerm += temp;
        }
    }
    
    for(e = 0; e < nEdges; e++)
    {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        
        for(s1 = 0; s1 < nStates[n1];s1++)
        {
            for(s2 = 0; s2 < nStates[n2]; s2++)
            {
                factorTerm += exp(edgeBel[s1+maxState*(s2+maxState*e)]) * (-edgePot[s1+maxState*(s2+maxState*e)] + edgeBel[s1+maxState*(s2+maxState*e)]);
            }
        }
    }
    logZ[0] = nodeTerm - factorTerm;
    
   /* Free memory */
    mxFree(prodMsgs);
    mxFree(oldMsgs);
    mxFree(newMsgs);
    mxFree(tmp);
    mxFree(prodMsgForUptree);
    mxFree(prodMsgForFactorMargs);
    mxFree(factorMarginals);    
}

double myLog(double x)
{
    if (x == 0)
        return MINUSINF;
    else
        return log(x);
}
