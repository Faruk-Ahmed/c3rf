function [msgs_from_xor, logZ_XOR] = get_msg_from_XOR_factor(msg_to_xor_factor, nnodes, nlabels, calc_margs)
%%
msgs_from_xor = zeros(nlabels*nnodes, 2);
logZ_XOR = 0;

for i = 1:nnodes
    Z = my_log(msg_to_xor_factor(nlabels*(i-1) + (1:nlabels)', 1)');
    O = my_log(msg_to_xor_factor(nlabels*(i-1) + (1:nlabels)', 2)');
    [~, msgs_from_XOR, logZ_XOR_contribution] = sum_cardinality_one_of_K(int32(nlabels), Z, O, int32(1), int32(1));
    msgs_from_xor(nlabels*(i-1) + (1:nlabels)', :) = exp(msgs_from_XOR) ./ repmat(sum(exp(msgs_from_XOR),2),1,2);
    logZ_XOR = logZ_XOR + logZ_XOR_contribution;
end

