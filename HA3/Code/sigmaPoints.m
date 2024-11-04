function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%

    n = length(x);
    P_sqrt = sqrtm(P);
    switch type        
        case 'UKF'
            
        % Compute sigma points using Unscented Transform
        SP = zeros(n,2*n+1);
        W = zeros(1,2*n+1);
        
        % initial
        W(1) = 1-n/3;
        SP(:,1)=x;
        
        for i = 2:(n+1)
            SP(:,i) = x + sqrt(n/(1-W(1)))*P_sqrt(:,i-1);
            SP(:,i+n) = x - sqrt(n/(1-W(1)))*P_sqrt(:,i-1);
            W(i)=(1-W(1))/(2*n);
            W(i+n)=(1-W(1))/(2*n);
        end
        %%%%%
        case 'CKF'
            
        % Compute sigma points using Cubature Rule
        SP = zeros(n,2*n);
        W = zeros(1,2*n);  
        
        
        for i=1:n
            SP(:, i) = x + sqrt(n)*P_sqrt(:,i);
            SP(:, i+n) = x - sqrt(n)*P_sqrt(:,i);
            W(i) = 1/(2*n);
            W(i+n) = 1/(2*n);         
        end
            
        otherwise
            error('Incorrect type of sigma point')
    end

end