function V = my_norm(A)

V = zeros(size(A));
for i = 1:size(A, 1)
    M1 = max(A(i, 1), A(i, 2));
    Z = M1 + my_log(exp(A(i,1) - M1) + exp(A(i,2) - M1));
    V(i, :) = A(i, :) - Z;
end
