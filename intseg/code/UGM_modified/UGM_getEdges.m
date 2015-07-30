function [edges] = UGM_getEdges(n,edgeStruct)
% From Mark Schmidt's UGM library here:
% http://www.cs.ubc.ca/~schmidtm/Software/UGM.html
edges = edgeStruct.E(edgeStruct.V(n):edgeStruct.V(n+1)-1);
