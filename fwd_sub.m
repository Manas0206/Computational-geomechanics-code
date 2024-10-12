function [y] = fwd_sub(L, P, b);
% Forward substitution: Ly = Pb
% input arguments: Lower triangular factor L
%                : Column vector of row-exchange indices P
%                : Column vector of right-hand side b
% output arguments: Column vector of intermediate solution y = Ux
%new b from p
n=length(b);
M=b;
for i=1:1:n
    k=P(i);
    b(i)=M(k);
end

y=zeros(n,1);
for i=1:n
    y(i,1)=b(i,1)/L(i,i);
    for j=(i+1):n
        b(j,1)=b(j,1)-L(j,i)*y(i,1);
    end
end
end

