#include <math.h>
#include "mex.h"
#include "matrix.h"

#define EPSILON log(1e-16)
#define MINUSINF -10000000

double myLog(double);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int D, i, j, d, ctr, lvl, lengthOfTopLevel, lengthOfNewTopLevel, *topLevel, *newTopLevel,
            *cards, *startIdxs, totalCumSum, start, end, m, dd, ch1, ch2, startCh1, startCh2,
            endCh1, endCh2, startP, endP, cardSum, calcMargs, hammingRadius, root;

    double *T, temp, *upMessages, *downMessages, *logNodePotsZeros, *logNodePotsOnes, *ch1Msg, *ch2Msg, 
            Z, *cardinalityFactor, *cardinalityFactorBelief, *downMsgFromZ, topFactorTerm,
            *down1Msg, *down2Msg, *downPMsg, nodeTerm, factorTerm, *factorBeliefs, 
            *nodeBeliefsForCh1, *nodeBeliefsForCh2, *rootNodeBeliefs, *myLogZ, *finalMsgs, *maxVal, maxVal1;
    
    /* Input: */
    D = ((int*)mxGetPr(prhs[0]))[0];
    logNodePotsZeros = (double*)mxGetPr(prhs[1]);
    logNodePotsOnes = (double*)mxGetPr(prhs[2]);
    calcMargs = ((int*)mxGetPr(prhs[3]))[0];
    hammingRadius = ((int*)mxGetPr(prhs[4]))[0];
    
    /* Output: */
    plhs[0] = mxCreateDoubleMatrix(D-1, 4, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(D, 2, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    /* make tree */
    T = mxGetPr(plhs[0]);
    finalMsgs = mxGetPr(plhs[1]);
    myLogZ = mxGetPr(plhs[2]);

    topLevel = mxCalloc(D, sizeof(int));
    lengthOfTopLevel = D;

    newTopLevel = mxCalloc(ceil((double)D/2), sizeof(int));
    lengthOfNewTopLevel = ceil((double)D/2);
    
    for (i = 0; i < D-1; i++)
        for(j = 0; j < 4; j++)
        {
            T[i+(D-1)*j] = 0;
        }
    
    ctr = 0;
    for (i = 0; i < lengthOfTopLevel; i++)
        topLevel[i] = i;
    for (i = 0; i < lengthOfNewTopLevel; i++)
        newTopLevel[i] = 0;
    lvl = 1;
    
    while (lengthOfNewTopLevel > 1)
    {
        for (d = 0; d < (int)lengthOfTopLevel/2; d++)
        {
            newTopLevel[d] = ctr+D;
            T[ctr + (D-1)*0] = topLevel[2*d];
            T[ctr + (D-1)*1] = topLevel[2*d+1];
            T[ctr + (D-1)*2] = lvl;
            T[ctr + (D-1)*3] = 0;
            
            if (T[ctr + (D-1)*0] < D)
                T[ctr + (D-1)*3] += 1;
            else
                T[ctr + (D-1)*3] += T[(int)T[ctr + (D-1)*0] - D + (D-1)*3];
            if (T[ctr + (D-1)*1] < D)
                T[ctr + (D-1)*3] += 1;
            else
                T[ctr + (D-1)*3] += T[(int)T[ctr + (D-1)*1] - D + (D-1)*3];
                
            ctr += 1;            
        }
        
        if ((lengthOfTopLevel % 2) == 1)
        {
            newTopLevel[lengthOfNewTopLevel-1] = topLevel[lengthOfTopLevel-1];
        }
        
        mxFree(topLevel);
        topLevel = mxCalloc(lengthOfNewTopLevel, sizeof(int));
        lengthOfTopLevel = lengthOfNewTopLevel;
        for (i = 0; i < lengthOfTopLevel; i++)
            topLevel[i] = newTopLevel[i];
        
        mxFree(newTopLevel);
        newTopLevel = mxCalloc(ceil((double)lengthOfTopLevel/2), sizeof(int));
        lengthOfNewTopLevel = ceil((double)lengthOfTopLevel/2);
        for (i = 0; i < ceil((double)lengthOfTopLevel/2); i++)
            newTopLevel[i] = 0;
        
        lvl += 1;
    }
    
    T[ctr + (D-1)*0] = topLevel[0];
    T[ctr + (D-1)*1] = topLevel[1];
    T[ctr + (D-1)*2] = lvl;
    
    if (T[ctr + (D-1)*0] < D)
        T[ctr + (D-1)*3] += 1;
    else
        T[ctr + (D-1)*3] += T[(int)T[ctr + (D-1)*0] - D + (D-1)*3];
    if (T[ctr + (D-1)*1] < D)
        T[ctr + (D-1)*3] += 1;
    else
        T[ctr + (D-1)*3] += T[(int)T[ctr + (D-1)*1] - D + (D-1)*3];

    
    /* Pass messages upwards: */
    cards = mxCalloc(2*D-1, sizeof(int));
    for (i = 0; i < 2*D-1; i++)
    {
        if (i < D)
            cards[i] = 1;
        else
            cards[i] = cards[(int)T[i-D + (D-1)*0]] + cards[(int)T[i-D + (D-1)*1]];;
    }
    
    cardSum = 0;
    for(i =  0; i < 2*D-1; i++)
        cardSum += cards[i]+1;
    
    startIdxs = mxCalloc(2*D, sizeof(int));
    totalCumSum = 0;
    for (i = 0; i < 2*D; i++)
    {
        temp = 0;
        for (j = 0; j < i; j++)
            temp += cards[j]+1;
        startIdxs[i] = temp;
        totalCumSum += temp;
    }

    upMessages = mxCalloc(cardSum, sizeof(double));
    downMessages = mxCalloc(cardSum, sizeof(double));
    
    for (d = 0; d < D; d=d+1)
    {
        upMessages[2*d] = logNodePotsZeros[d];
        upMessages[2*d+1] = logNodePotsOnes[d];
    }
    for (d = 2*D; d < cardSum; d++)
        upMessages[d] = 0;
    
    for (m = 0; m < D-1; m++)
    {
        dd = D + m;
        ch1 = T[m + (D-1)*0];
        ch2 = T[m + (D-1)*1];
        
        startCh1 = startIdxs[ch1];
        endCh1   = startCh1 + cards[ch1] + 1;
        startCh2 = startIdxs[ch2];
        endCh2   = startCh2 + cards[ch2] + 1;
        startP   = startIdxs[dd]; 
        endP     = startP + cards[dd] + 1;
        
        ch1Msg = mxCalloc(cards[ch1]+1, sizeof(double));
        ch2Msg = mxCalloc(cards[ch2]+1, sizeof(double));
        maxVal = mxCalloc(cards[ch1]+cards[ch2]+1, sizeof(double));

        for (i = 0; i <= cards[ch1]; i++)
            ch1Msg[i] = upMessages[startCh1+i];

        for (i = 0; i <= cards[ch2]; i++)
            ch2Msg[i] = upMessages[startCh2+i];

        for (i = 0; i <= cards[dd]; i++)
            maxVal[i] = ch1Msg[0] + ch2Msg[0];

        for (i = 0; i <= cards[ch1]; i++)
        {
            for (j = 0; j <= cards[ch2]; j++)
            {
                temp = ch1Msg[i] + ch2Msg[j];
                if (maxVal[i+j] < temp)
                    maxVal[i+j] = temp;
            }
        }
                
        for (i = 0; i <= cards[ch1]; i++)
        {
            for (j = 0; j <= cards[ch2]; j++)
            {
                temp = ch1Msg[i] + ch2Msg[j];
                upMessages[startP+i+j] += exp(temp - maxVal[i+j]);
            }
        }

        for (i = startP; i <= startP+cards[dd]; i++)
            upMessages[i] = maxVal[i-startP] + myLog(upMessages[i]);
        
        mxFree(maxVal);

        /* normalizing template: use this for all norm. */
        maxVal1 = upMessages[startP];
        for (i = startP; i <= startP+cards[dd]; i++)
            if (upMessages[i] > maxVal1)
                maxVal1 = upMessages[i];
        
        Z = maxVal1;
        temp = 0;
        for (i = startP; i <= startP+cards[dd]; i++)
            temp += exp(upMessages[i] - maxVal1);

        Z += myLog(temp);
        
        for(i = startP; i <= startP+cards[dd]; i++)
            upMessages[i] -= Z;
        /* normalizing template */
        
        mxFree(ch1Msg);
        mxFree(ch2Msg);
    }    

    /* Pass messages down from the Cardinality Potential: */
    cardinalityFactor = mxCalloc(D+1, sizeof(double));
    for (i = 0; i < D+1; i++)
        cardinalityFactor[i] = MINUSINF;
    cardinalityFactor[1]= 0;

    downMsgFromZ = mxCalloc(D+1, sizeof(double));
    for (i = 0; i < D+1; i++)
        downMsgFromZ[i] = cardinalityFactor[i];
    
    topFactorTerm = 0;
    
    if (calcMargs)
    {
        cardinalityFactorBelief = mxCalloc(D+1, sizeof(double));

        for (i = 0; i < D+1; i++)
            cardinalityFactorBelief[i] = upMessages[startP+i] + cardinalityFactor[i];
        /* */
        maxVal1 = cardinalityFactorBelief[0];
        for (i = 0; i < D+1; i++)
            if (cardinalityFactorBelief[i] > maxVal1)
                maxVal1 = cardinalityFactorBelief[i];
        
        Z = maxVal1;
        temp = 0;
        for (i = 0; i < D+1; i++)
        {
            temp += exp(cardinalityFactorBelief[i] - maxVal1);
        }
        Z += myLog(temp);
        
        for(i = 0; i < D+1; i++)
        {
            cardinalityFactorBelief[i] -= Z;
        }
        /* */        

        for ( i = 0; i < D+1; i++)
        {
            topFactorTerm += exp(cardinalityFactorBelief[i])*(-cardinalityFactor[i] + cardinalityFactorBelief[i]);
        }
    }
            
    /* Now pass the message down the tree: */
    nodeTerm = 0;
    factorTerm = 0;
    
    for (i = 0; i < cardSum; i++)
        downMessages[i] = 0;
    
    for (i = 0; i <= cards[dd]; i++)
        downMessages[startP+i] = downMsgFromZ[i];
    
    root = 1;
    
    for (m = D-2; m >= 0; m--)
    {
        dd = D + m;
        ch1 = (int)T[m + (D-1)*0];
        ch2 = (int)T[m + (D-1)*1];
        
        startCh1 = startIdxs[ch1];
        endCh1   = startCh1 + cards[ch1] + 1;
        startCh2 = startIdxs[ch2];
        endCh2   = startCh2 + cards[ch2] + 1;
        startP   = startIdxs[dd];
        endP     = startP + cards[dd] + 1;

        downPMsg = mxCalloc(cards[dd]+1, sizeof(double));
        for (i = 0; i <= cards[dd]; i++)
            downPMsg[i] = downMessages[startP+i];
        
        ch1Msg = mxCalloc(cards[ch1]+1, sizeof(double));
        ch2Msg = mxCalloc(cards[ch2]+1, sizeof(double));
        
        for (i = 0; i <= cards[ch1]; i++)
            ch1Msg[i] = upMessages[startCh1+i];

        for (i = 0; i <= cards[ch2]; i++)
            ch2Msg[i] = upMessages[startCh2+i];
        
        down1Msg = mxCalloc(cards[ch1]+1, sizeof(double));
        down2Msg = mxCalloc(cards[ch2]+1, sizeof(double));

        maxVal = mxCalloc(cards[ch1]+1, sizeof(double));
        for (i = 0; i <= cards[ch1]; i++)
            maxVal[i] = ch2Msg[0] + downPMsg[0];
            
        for (i = 0; i <= cards[ch1]; i++)
        {
            for (j = 0; j <= cards[ch2]; j++)
            {
                temp = ch2Msg[j] + downPMsg[i+j];
                if (temp > maxVal[i])
                    maxVal[i] = temp;
            }
        }            
            
        for (i = 0; i <= cards[ch1]; i++)
            for (j = 0; j <= cards[ch2]; j++)
                down1Msg[i] += exp(ch2Msg[j] + downPMsg[i+j] - maxVal[i]);

        for (i = 0; i <= cards[ch1]; i++)
            down1Msg[i] = maxVal[i] + myLog(down1Msg[i]);
        
        mxFree(maxVal);
        
        /* normalize down1Msg */
        maxVal1 = down1Msg[0];
        for (i = 0; i <= cards[ch1]; i++)
            if (down1Msg[i] > maxVal1)
                maxVal1 = down1Msg[i];

        Z = maxVal1;
        temp = 0;
        for (i = 0; i <= cards[ch1]; i++)
            temp += exp(down1Msg[i] - maxVal1);
        Z += myLog(temp);
        
        /*mexPrintf("ch1 = %d: ", ch1);*/
        for(i = 0; i <= cards[ch1]; i++)
        {
            down1Msg[i] -= Z;
            /*mexPrintf("%f ", down1Msg[i]);*/
        }
        /*mexPrintf("\n");*/
        /* */
        
        maxVal = mxCalloc(cards[ch2], sizeof(double));
        for (i = 0; i <= cards[ch2]; i++)
            maxVal[i] = ch1Msg[0] + downPMsg[0];

        for (i = 0; i <= cards[ch2]; i++)
        {
            for (j = 0; j <= cards[ch1]; j++)
            {
                temp = (ch1Msg[j] + downPMsg[i+j]);
                if (maxVal[i] < temp)
                    maxVal[i] = temp;
            }
        }
            
        for (i = 0; i <= cards[ch2]; i++)
            for (j = 0; j <= cards[ch1]; j++)
                down2Msg[i] += exp(ch1Msg[j] + downPMsg[i+j] - maxVal[i]);
                
        for (i = 0; i <= cards[ch2]; i++)
            down2Msg[i] = maxVal[i] + myLog(down2Msg[i]);
        mxFree(maxVal);
        /* normalize down2Msg */
        maxVal1 = down2Msg[0];
        for (i = 0; i <= cards[ch2]; i++)
        {
            if (down2Msg[i] > maxVal1)
                maxVal1 = down2Msg[i];
        }
        
        Z = maxVal1;
        temp = 0;
        for (i = 0; i <= cards[ch2]; i++)
        {
            temp += exp(down2Msg[i] - maxVal1);
        }
        Z += myLog(temp);
        
        for(i = 0; i <= cards[ch2]; i++)
        {
            down2Msg[i] -= Z;
        }
        /* */

        for (i = 0; i <= cards[ch1]; i++)
            downMessages[startCh1+i] = down1Msg[i];
        for (i = 0; i <= cards[ch2]; i++)
            downMessages[startCh2+i] = down2Msg[i];
                
        if (calcMargs)
        {
            factorBeliefs = mxCalloc((cards[ch1]+1)*(cards[ch2]+1), sizeof(double));
            temp = 0;
            for ( i = 0; i <= cards[ch1]; i++)
            {
                for (j = 0; j <= cards[ch2]; j++)
                {
                    factorBeliefs[(int)temp] = downPMsg[i+j] + ch1Msg[i] + ch2Msg[j];
                    temp += 1;
                }
            }

            /* */
            maxVal1 = factorBeliefs[0];
            for (i = 0; i < (cards[ch1]+1)*(cards[ch2]+1); i++)
                if (factorBeliefs[i] > maxVal1)
                    maxVal1 = factorBeliefs[i];

            Z = maxVal1;
            temp = 0;
            for (i = 0; i < (cards[ch1]+1)*(cards[ch2]+1); i++)
                temp += exp(factorBeliefs[i] - maxVal1);
            Z += myLog(temp);

            for(i = 0; i <= (cards[ch1]+1)*(cards[ch2]+1); i++)
                factorBeliefs[i] -= Z;
            /* */
            
            for (i = 0; i < (cards[ch1]+1)*(cards[ch2]+1); i++)
                factorTerm += exp(factorBeliefs[i])*factorBeliefs[i];
            
            if (ch1 >= D)
            {                
                nodeBeliefsForCh1 = mxCalloc(cards[ch1]+1, sizeof(double));
            
                for (i = 0; i <= cards[ch1]; i++)
                    nodeBeliefsForCh1[i] = ch1Msg[i] + down1Msg[i];
                /* */
                maxVal1 = nodeBeliefsForCh1[0];
                for (i = 0; i <= cards[ch1]; i++)
                    if (nodeBeliefsForCh1[i] > maxVal1)
                        maxVal1 = nodeBeliefsForCh1[i];

                Z = maxVal1;
                temp = 0;
                for (i = 0; i <= cards[ch1]; i++)
                    temp += exp(nodeBeliefsForCh1[i] - maxVal1);
                Z += myLog(temp);

                for(i = 0; i <= cards[ch1]; i++)
                {
                    nodeBeliefsForCh1[i] -= Z;
                }
                /* */
                for (i = 0; i <= cards[ch1]; i++)
                    nodeTerm += exp(nodeBeliefsForCh1[i])*nodeBeliefsForCh1[i];
                mxFree(nodeBeliefsForCh1);
            }
            
            if (ch2 >= D)
            {
                nodeBeliefsForCh2 = mxCalloc(cards[ch2]+1, sizeof(double));
                for (i = 0; i <= cards[ch2]; i++)
                {
                    nodeBeliefsForCh2[i] = ch2Msg[i] + downMessages[startCh2+i];
                }
                /* */
                maxVal1 = nodeBeliefsForCh2[0];
                for (i = 0; i <= cards[ch2]; i++)
                    if (nodeBeliefsForCh2[i] > maxVal1)
                        maxVal1 = nodeBeliefsForCh2[i];

                Z = maxVal1;
                temp = 0;
                for (i = 0; i <= cards[ch2]; i++)
                    temp += exp(nodeBeliefsForCh2[i] - maxVal1);
                Z += myLog(temp);

                for(i = 0; i <= cards[ch2]; i++)
                {
                    nodeBeliefsForCh2[i] -= Z;
                }
                /* */
                for (i = 0; i <= cards[ch2]; i++)
                    nodeTerm += exp(nodeBeliefsForCh2[i])*nodeBeliefsForCh2[i];
                mxFree(nodeBeliefsForCh2);
            }
            
            if (root)
            {
                rootNodeBeliefs = mxCalloc(cards[dd]+1, sizeof(double));
                for (i = 0; i <= cards[dd]; i++)
                {
                    rootNodeBeliefs[i] = upMessages[startP+i] + downPMsg[i];
                }
                /* */
                maxVal1 = rootNodeBeliefs[0];
                for (i = 0; i <= cards[dd]; i++)
                    if (rootNodeBeliefs[i] > maxVal1)
                        maxVal1 = rootNodeBeliefs[i];

                Z = maxVal1;
                temp = 0;
                for (i = 0; i <= cards[dd]; i++)
                    temp += exp(rootNodeBeliefs[i] - maxVal1);
                Z += myLog(temp);

                for(i = 0; i <= cards[dd]; i++)
                {
                    rootNodeBeliefs[i] -= Z;
                }
                /* */
                for (i = 0; i <= cards[dd]; i++)
                    nodeTerm += exp(rootNodeBeliefs[i])*rootNodeBeliefs[i];
                root = 0;
                mxFree(rootNodeBeliefs);
            }
        }

        mxFree(ch1Msg);
        mxFree(ch2Msg);
        mxFree(downPMsg);
        
    }

    for (i = 0; i < D; i++)
    {
        finalMsgs[i + D*0] = downMessages[2*i];
        finalMsgs[i + D*1] = downMessages[2*i+1];
    }
    myLogZ[0] = nodeTerm - factorTerm - topFactorTerm;
    
    mxFree(upMessages);
    mxFree(downMessages);
    mxFree(cardinalityFactor);
    mxFree(downMsgFromZ);
    mxFree(topLevel);
    mxFree(newTopLevel);
}

double myLog(double x)
{
    if (x == 0)
        return MINUSINF;
    else
        return log(x);
}
