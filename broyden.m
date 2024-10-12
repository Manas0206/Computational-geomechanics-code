function NLSolve = broyden(NLSolve)
% Solve for eigenvalues and eigenvectors of symmetric/positive-definite 
% nxn square matrix A using Broyden algorithm
% Store your converged solution in: NLSolve.lambda_num = lambda_converged
%                                   NLSolve.x_num = x_converged

% Write code for your function here
fprintf('Solving using Broyden algorithm:\n')

u=NLSolve.x0; 
a=NLSolve.A;
n=length(NLSolve.A);
k=n+1;
max_iter=30;
lambda=NLSolve.lambda0;

%resudual calculation (4*1) if A is (3*3) as x'x=1 equation used
R0=a*u-lambda*u;
R0(k,1)=-1;
for i=1:n
    R0(k,1)=R0(k,1)+u(i,1)*u(i,1);
end

    
%calculation of jacobian (after derivative)
a=a-lambda*eye(n);
for j=1:n
    a(j,k)=-1*u(j,1);
end

for j=1:n
    a(k,j)=2*u(j,1);
end
a(k,k)=0; 
J0=a;



for z=0:max_iter
   

    if(z==0)
        %%LU 
    n1=length(J0);          %length of a matrix
    L=eye(n1);             %identity matrix
    U=J0;
    for i=1:n1-1  %loop over column
        for j=i+1:n1 %loop over row
            L(j,i)=U(j,i)/U(i,i);
            U(j,1:n1)=U(j,1:n1)-L(j,i)*U(i,1:n1); %small u uper triangular matrix
        
        end
    end

    % forward substitution
    y = fwd_sub(L, R0);

    % back substitution
    delta_u = -back_sub(U, y);
    
    
    else
        
    %resudual calculation (4*1) if A is (3*3) as x'x=1 equation used
    R=a*u-lambda*u; % u is n*1 size
    R(k,1)=-1;
    for i=1:n
         R(k,1)=R(k,1)+u(i,1)*u(i,1);
    end
    
    
  
      mm=-1*(R-R0-J0*delta_u);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      ss=delta_u'/(delta_u'*delta_u);
      J=J0-mm*ss; % j(n)=j(n-1)-u*v'  
       
      
         %%LU 
    n1=length(J);          %length of a matrix
    L=eye(n1);             %identity matrix
    U=J;
    for i=1:n1-1  %loop over column
        for j=i+1:n1 %loop over row
            L(j,i)=U(j,i)/U(i,i);
            U(j,1:n1)=U(j,1:n1)-L(j,i)*U(i,1:n1); %small u uper triangular matrix
        
        end
    end

    % forward substitution
    y = fwd_sub(L, R);

    % back substitution
    delta_u = -back_sub(U, y);
      
        
        
        
     R0=R; %storing R to R0
     J0=J; %storing J to J0
        
    end
    %updating u
    u(k,1)=lambda;   %passing eigen value through u at n+1 row col-1
   
    u=u+delta_u ; %u is (n+1)*1 size   %%%%%%%%%%%%%%%%55%%%%%%%%%%%%%
    lambda=u(k,1);
    u=u(1:n,1);
    a=NLSolve.A;
    
    if abs(norm((delta_u),2))<NLSolve.tol
        break;
    end
end
    NLSolve.lambda_num=lambda;
    NLSolve.x_num=u;
end

