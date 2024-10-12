function [L,U,P] = lu_factorization_partial_pivoting(A)
% LU factorization of an nxn square matrix A with partial pivoting
% input arguments: nxn square matrix A
% output arguments: Lower triangular factor L
%                   Upper triangular factor U
%                   Column vector of row-exchange indices P
n=length(A);          %length of a matrix

%CALCULATION FOR p
for i=1:1:n
 P(i,1)=i;
end
k=A;
for i=1:n-1
     
    for j=i+1:n
        
        if k(j,i)>k(i,i)
            
            c=k(i,:);
            k(i,:)=k(j,:);
            k(j,:)=c;
            C=P(i);
            M=P(j);
            P(i)=M;
            P(j)=C;
             
        end
    end
end


%LU
A1=k;
L=eye(n); 
U=A1;
for i=1:n-1  %loop over column
    for j=i+1:n %loop over row
        L(j,i)=U(j,i)/U(i,i);
        U(j,1:n)=U(j,1:n)-L(j,i)*U(i,1:n); %small u uper triangular matrix
        
    end
end

end




