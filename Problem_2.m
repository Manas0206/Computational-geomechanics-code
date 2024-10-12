
n = 3;
A = rand(n,n);
b = rand(n,1);

% Rank-one update of A
u = rand(n,1); v = rand(n,1);
B = A - u*v';

% Solve Bx = b
% You may use the functions you wrote for
% LU factorization, forward/back substitution in Problem 1 or Homeworks 

%%LU of A
n=length(A);          %length of a matrix
L=eye(n);             %identity matrix
U=A;
for i=1:n-1  %loop over column
    for j=i+1:n %loop over row
        L(j,i)=U(j,i)/U(i,i);
        U(j,1:n)=U(j,1:n)-L(j,i)*U(i,1:n); %small u uper triangular matrix
        
    end
end


%Az=u
% forward substitution
y = fwd_sub(L, u);

% back substitution
x = back_sub(U, y);
z=x;

%Am=b
% forward substitution
y = fwd_sub(L, b);

% back substitution
x = back_sub(U, y);
m=x;

alpha=(v'*m)/(1-v'*z);
x=m+alpha*z;

fprintf('Norm of x_err = %3.2e\n', norm((x-(B\b)),2));
