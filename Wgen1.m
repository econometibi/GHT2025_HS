function W = Wgen1(n,m)
% circular network
% n is the dimension of W.
% 2*m is the number of connections
W1 = zeros(n,3*n);
for i = 1:n
    W1(i,n+i-m:n+i-1) = 1;
    W1(i,n+i+1:n+i+m) = 1;
end
W2 = W1(:,1:n)+W1(:,n+1:2*n)+W1(:,2*n+1:3*n);
W = zeros(n);
W(W2>0) = 1;
% W = sparse(W);