function z  = create2dVAR(T, p, sigmaA, isCausalInteraction)
% generate T data points from a two-dimensional stationary VAR z = [x;y] 
% of lag-order p
% IN    T                   - number of data points 
%       p                   - lag order 
%       sigmaA              - standard deviation of AR coefficients
%       isCausalInteraction - if true: x Granger-causes y
%                             if false: x and y are independent 0
% OUT   z                   - [x; y] signal (2 x T matrix)  

%% for lag order p = 0, simply return Gaussian distributed random variables
if p == 0
    Sigma = diag(rand(1,2)); 
    z= transpose(mvnrnd([0; 0], Sigma,  T));
    return
end

%% create AR coefficients
%create 2px2p coefficient matrix boldA of the VAR(1) representation of the VAR(p) process
%the process is stationary if all eigenvalues of boldA are smaller than 1
boldA = 2*ones(2*p, 2*p); 
A = cell(1, p); 
while not(isempty(find(abs(eig(boldA)) > 1))) %require stationarity
    for j = 1:p
        A{j} = sigmaA*randn(2,2); %draw AR coefficients from gaussian distribution
        A{j}(1, 2) = 0;  %zero interaction from y to x
        if not(isCausalInteraction)
            A{j}(2, 1) = 0; %zero interaction from x to y 
        end
    end
    boldA(1:2, :) = cell2mat(A); 
    boldA(3:end,1:end-2) = eye(2*(p-1)); 
    boldA(3:end,end-1:end) = zeros(2*(p-1), 2); 
end

%% generate T + 100 samples 
Sigma = diag(rand(1,2)); 
z= transpose(mvnrnd([0; 0], Sigma,  T + 100)); 
for t = p+1:T+100   
    for j = 1:p 
        z(:,t) = z(:,t) + A{j} * z(:,t-j); 
    end 
end

%% cut of first 100 samples
z = z(:,101:end); 