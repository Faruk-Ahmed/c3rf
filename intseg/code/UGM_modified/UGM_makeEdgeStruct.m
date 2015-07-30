function [edgeStruct] = UGM_makeEdgeStruct(el,nStates,nNodes,useMex,maxIter)
% Modified version from Mark Schmidt's UGM library here:
% http://www.cs.ubc.ca/~schmidtm/Software/UGM.html 

if nargin < 3
    useMex = 1;
end
if nargin < 4
    maxIter = 100;
end

[V,E] = UGM_makeEdgeVE(el,nNodes,useMex);


edgeStruct.edgeEnds = el;
edgeStruct.V = V;
edgeStruct.E = E;
edgeStruct.nNodes = nNodes;
edgeStruct.nEdges = size(el,1);

% Handle other arguments
if isscalar(nStates)
   nStates = repmat(nStates,[nNodes 1]);
end
edgeStruct.nStates = int32(nStates(:));
edgeStruct.useMex = useMex;
edgeStruct.maxIter = int32(maxIter);

