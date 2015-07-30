function [MASS, NODE_MARGINALS, EDGE_MARGINALS] = get_mass(nodePot, edgePot, edgeStruct, SOLUTION, HAMMING_RADIUS, noHOP, pf, lambda, solind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = '.';

% ne
nnodes = size(nodePot, 1);

ct_iter = 0;

ZERO_NODE_INDS = find(SOLUTION == 0);
ONE_NODE_INDS = find(SOLUTION == 1);

msgProd = zeros(size(nodePot));
for i = 1:size(nodePot, 1)
    M1 = max(nodePot(i, 1), nodePot(i, 2));
    Z = M1 + my_log(exp(nodePot(i,1)-M1) + exp(nodePot(i,2)-M1));
    msgProd(i, :) = [nodePot(i,1) - Z, nodePot(i,2) - Z];
end

NORM_NODE_POT = msgProd;
messages_DOWNTREE = ones(nnodes, 2);
msgProd_old = zeros(nnodes, 2);
messages_SP_old = -1;
oldMsgSumChange = 0;

useDamping = 0;
ALPHA = 0.3;
calc_margs = 1;
breakNow = 0;

%% Iterate till convergence of messages:
if ~noHOP
while 1

    ct_iter = ct_iter + 1;
    
    %% Perform 1 round of message passing (sum-product) on pairwise graph:
    [nodeBel, edgeBel, logZ_SP, new_msg, nodeBel_SP_UPTREE] = LBP_log(msgProd, ...
                                                                   NORM_NODE_POT, ...
                                                                   nodePot, ...
                                                                   my_log(messages_DOWNTREE), ...
                                                                   0*int32(ones(nnodes, 1)), ...
                                                                   edgePot, ... 
                                                                   int32(edgeStruct.edgeEnds),int32(edgeStruct.nStates),int32(edgeStruct.V),int32(edgeStruct.E),edgeStruct.maxIter, ...
                                                                   int32(calc_margs), int32(useDamping), messages_SP_old, ALPHA);

    nodeBel = exp(nodeBel);
    nodeBel = nodeBel ./ repmat(sum(nodeBel,2), 1, 2);
    nodeBel = nodeBel ./ repmat(sum(nodeBel,2), 1, 2);
    
    edgeBel = exp(edgeBel);
    nodeBel_SP_UPTREE = exp(nodeBel_SP_UPTREE);

    %% Impose cardinality potential:
    nodeBel_SP_UPTREE = nodeBel_SP_UPTREE ./ repmat(sum(nodeBel_SP_UPTREE, 2), 1, 2);
    [~, messages_downtree, logZ_card] = sum_cardinality(int32(nnodes), ...
                                                       [my_log(nodeBel_SP_UPTREE(ZERO_NODE_INDS, 1)); my_log(nodeBel_SP_UPTREE(ONE_NODE_INDS, 2))],...
                                                       [my_log(nodeBel_SP_UPTREE(ZERO_NODE_INDS, 2)); my_log(nodeBel_SP_UPTREE(ONE_NODE_INDS, 1))],...
                                                       int32(calc_margs), ...
                                                       int32(HAMMING_RADIUS));

    messages_DOWNTREE(ZERO_NODE_INDS, :) = exp(messages_downtree(1:length(ZERO_NODE_INDS), [1, 2]));
    messages_DOWNTREE(ONE_NODE_INDS, :) = exp(messages_downtree(length(ZERO_NODE_INDS)+1:end, [2, 1]));
    messages_DOWNTREE = messages_DOWNTREE ./ repmat(sum(messages_DOWNTREE, 2), 1, 2);                 

    %%
    msgProd = messages_DOWNTREE .* exp(NORM_NODE_POT);
    msgProd = my_log(msgProd ./ repmat(sum(msgProd, 2), 1, 2));

    %% check for convergence:
    msgSumChange = sum(abs(msgProd(:) - msgProd_old(:)));
    NODEBEL_NOTCONVERGED = msgSumChange > 1e-3;

    if ct_iter == 30 %abs(msgSumChange - oldMsgSumChange) < 1e-4
        useDamping = 1;
    end
    
    msgProd_old = msgProd;
    messages_SP_old = new_msg;
    oldMsgSumChange = msgSumChange;

    %% Compute masses if converged:
    if ~NODEBEL_NOTCONVERGED
        if breakNow
            break;
        end
        breakNow = 1;
    end    
    
    if ct_iter > 5000
        fprintf('not converged: pf = %d,  lambda = %.2f,  solind = %d\n', pf, lambda, solind);
        break;
    end
end
else
    logZ_card = 0;
    [nodeBel, edgeBel, logZ_SP, new_msg, nodeBel_SP_UPTREE] = LBP_log(msgProd, ...
                                                               NORM_NODE_POT, ...
                                                               nodePot, ...
                                                               my_log(messages_DOWNTREE), ...
                                                               0*int32(ones(nnodes, 1)), ...
                                                               edgePot, ... 
                                                               int32(edgeStruct.edgeEnds),int32(edgeStruct.nStates),int32(edgeStruct.V),int32(edgeStruct.E),edgeStruct.maxIter, ...
                                                               int32(calc_margs), int32(useDamping), messages_SP_old, ALPHA);

end

MASS = logZ_card + logZ_SP;
NODE_MARGINALS = nodeBel;
EDGE_MARGINALS = edgeBel;

end
