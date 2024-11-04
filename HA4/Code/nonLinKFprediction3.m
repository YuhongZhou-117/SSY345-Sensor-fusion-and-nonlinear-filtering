function [x, P] = nonLinKFprediction3(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%
    n = length(x);
    switch type
        case 'EKF'
            
            [fx,Fx] = f(x);
            
            x = fx;
            P = Fx*P*Fx' + Q;
            
        case 'UKF'
            
            % Compute sigma points
            [SP,W] = sigmaPoints(x, P, 'UKF');
            
            x_pre = zeros(n,1); % predicted state mean
            P_pre = Q; % predicted state covariance
            
            for i = 1:length(W)
                x_pre = x_pre + f(SP(:,i))*W(i);
            end
            
            for k = 1:length(W)
                P_pre = P_pre + (f(SP(:,k))-x_pre)*(f(SP(:,k))-x_pre)'*W(k);
            end
            
            x = x_pre;
            P = P_pre;
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
                
        case 'CKF'
            
            % Compute sigma points
            [SP,W] = sigmaPoints(x, P, 'CKF');
            
            x_pre = zeros(n,1); % predicted state mean
            P_pre = Q; % predicted state covariance
            
            for i = 1:length(W)
                x_pre = x_pre + f(SP(:,i))*W(i);
            end
            
            for k = 1:length(W)
                P_pre = P_pre + (f(SP(:,k))-x_pre)*(f(SP(:,k))-x_pre)'*W(k);
            end
            
            x = x_pre;
            P = P_pre;
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end