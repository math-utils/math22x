function [L,U,P] = RLU(A)
% Written by Jonathan Ong (jongleberry@gmail.com) for UC Berkeley - 
% Math 221: Matrix Computations/Numerical Linear Algebra
%
% This program implements the recursive LU decomposition algorithm.
% The code is written for educational purposes and is not optimized.
%
% This code assumes that A is nxm with n > m and m = 2^k for some integer k
%
% Note that this code is not optimized for memory transfer
% Passing through indexes or pointers to the actual matrix instead of
% passing the actual matrix and not creating new variables to divide A
% would enhance the programs efficiency. 
%
% Also note that A = Perm\(L*U)
% Where Perm is defined by the following two lines:
%   Perm = eye(n);
%   Perm = Perm(P,:);

[n,m] = size(A);
if m == 1
    P = 1:n;
    [NULL,index] = max(A(:,1));
    if length(index) ~= 1
        index = index(1);
    end
    A([1,index],:) = A([index,1],:);
    P([1,index]) = P([index,1]);
    L = A/A(1);
    U = A(1);
else
    A11 = A(1:m/2,1:m/2);
    A21 = A(m/2+1:n,1:m/2);
    [L1,U1,P1] = RLU([A11;A21]);
    A = A(P1,:);
    A12 = A(1:m/2,m/2+1:m);
    A22 = A(m/2+1:n,m/2+1:m);
    L11 = L1(1:m/2,:);
    L12 = L1(m/2+1:n,:);
    A12 = L11\A12;
    A22 = A22-L12*A12;
    [L2,U2,P2] = RLU(A22);
    P = P1([1:m/2,m/2+P2]);
    L1(m/2+1:end,:) = L1(m/2+P2,:);
    L = [L1,[zeros(size(A12));L2]];
    U = [U1,A12;zeros(size(U1)),U2];
end