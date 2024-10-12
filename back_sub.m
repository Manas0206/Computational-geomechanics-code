function [y] = back_sub(U, y)
% Back substitution to solve Ux = y
% input arguments: Upper triangular factor U
%                : Column vector of intermediate solutions y
% output arguments: Column vector of final solution stored in y
% To minimize storage, y and x are stored in the same vectors
n=length(y);
%x=zeros(n,1);
for a=n:-1:1
    y(a,1)=y(a,1)/U(a,a);
    for b=(a-1):-1:1
        y(b,1)=y(b,1)-U(b,a)*y(a,1);
    end
end
end

