
clear; clc;

% Specify size of square matrix A
n = 3;
B =[0.9133 0.5383 0.4427;0.1524 0.9961 0.1067;0.8258 0.0782 0.9619];%rand(3,3);% create a matrix B of random numbers

% Create a data structure to group together several variables
% NLSolve.A = Symmetric/Positive Definite square matrix of random numbers and size nxn
% NLSolve.NR = true --> Use N-R method;
% NLSolve.Broyden = true --> Use Broyden method;

% NLSolve.x(:,1) = Eigenvector (nx1) at (k)th iteration
% NLSolve.x(:,2) = Eigenvector (nx1) at (k+1)th iteration
% NLSolve.lambda(:,1) = Eigenvalue at (k)th iteration
% NLSolve.lambda(:,2) = Eigenvalue at (k+1)th iteration

% NLSolve.max_iter = Maximum number of iterations beyond which the algorithm quits
% NLSolve.tol = Convergence tolerance
% NLSolve.x0 = Initial guess for eigenvector
% NLSolve.lambda0 = Initial guess for eigenvalue
% NLSolve.lambda_num = Converged numerical solution for eigenvalue lambda
% NLSolve.x_num = Converged numerical solution for eigenvalue x
NLSolve = struct('A',B'*B,'NR',false,'Broyden',true,'x',zeros(n,2),'lambda',zeros(1,2), ...
'max_iter',30,'tol',1e-12,'x0',zeros(n,1),'lambda0',0,'lambda_num',0,'x_num',zeros(n,1));

% Initial guess
NLSolve.x0 = rand(n,1); % Initial guess for x
NLSolve.x0 = (NLSolve.x0)/(norm(NLSolve.x0,2)); % Normalize the vector
NLSolve.lambda0 = NLSolve.x0'*(NLSolve.A)*NLSolve.x0; % Initial guess for lambda = x0'*A*x (Rayleigh quotient)

if(NLSolve.NR)
    % Solve using N-R method
    NLSolve = newton_raphson(NLSolve);
else
    % Solve using Broyden method
    % NLSolve.x(:,1) = Eigenvector (nx1) at (k-1)th iteration
    % NLSolve.x(:,2) = Eigenvector (nx1) at (k)th iteration
    % NLSolve.x(:,3) = Eigenvector (nx1) at (k+1)th iteration
    % NLSolve.lambda(:,1) = Eigenvalue at (k-1)th iteration
    % NLSolve.lambda(:,2) = Eigenvalue at (k)th iteration
    % NLSolve.lambda(:,3) = Eigenvalue at (k+1)th iteration
    
    NLSolve.lambda = zeros(1,3);
    NLSolve.x = zeros(n,3);
    NLSolve = broyden(NLSolve);
end

% Verification
% Verify solution against MATLAB's inbuilt-function
[V, D] = eig(NLSolve.A);
[r,c] = find(abs(D-NLSolve.lambda_num)<1e-6);
fprintf('Error in computed eigenvalue: lambda_m-lambda = %3.2f\n', (D(r,c)-NLSolve.lambda_num));
fprintf('Error in computed eigenvector: x_m-x: %3.2f\n', norm((abs(V(:,c))-abs(NLSolve.x_num)),2));