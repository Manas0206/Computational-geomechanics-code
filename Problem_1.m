
clear; clc;
% Specify the size of n of a square matrix
n = 10;
% Generate a matrix A, and a vector b of random numbers
A = rand(n,n);
b = rand(n,1);

% LU factorization with partial pivoting
[L, U, P] = lu_factorization_partial_pivoting(A);

% forward substitution
x = fwd_sub(L, P, b);

% back substitution
x = back_sub(U, x);

% verification
% for checking the code, an nxn matrix Pbar is formed from the vector P
% if correct factorization is obtained, Pbar*A = L*U
Pbar = zeros(n,n);
for iRow=1:n
    Pbar(iRow,P(iRow)) = 1;
end
fprintf('Norm of A_err = %3.2e\n', norm((Pbar*A-L*U),2));
fprintf('Norm of x_err = %3.2e\n', norm((x-(A\b)),2));
