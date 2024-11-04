function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, ...
                                     Ps_kplus1, ...
                                     xf_k, ... 
                                     Pf_k, ...
                                     xp_kplus1, ...
                                     Pp_kplus1, ...
                                     f, ...
                                     T, ...
                                     sigmaPoints, ...
                                     type)
%NONLINRTSSUPDATE Calculates mean and covariance of smoothed state
% density, using a non-linear Gaussian model.
%
%Input:
%   xs_kplus1   Smooting estimate for state at time k+1
%   Ps_kplus1   Smoothing error covariance for state at time k+1
%   xf_k        Filter estimate for state at time k
%   Pf_k        Filter error covariance for state at time k
%   xp_kplus1   Prediction estimate for state at time k+1
%   Pp_kplus1   Prediction error covariance for state at time k+1
%   f           Motion model function handle
%   T           Sampling time
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xs          Smoothed estimate of state at time k
%   Ps          Smoothed error convariance for state at time k

    n = size(Pf_k,1);
    switch type
        case 'EKF'
            [~,Ak] = f(xf_k,T);
            
            % smoothing gain
            Gk = Pf_k*Ak'*inv(Pp_kplus1);
            
            xs = xf_k + Gk*(xs_kplus1 - xp_kplus1);
            
            Ps = Pf_k - Gk*(Pp_kplus1 - Ps_kplus1)*Gk';
            
        case 'UKF'
            % calculate sigma points
            [SP,W] = sigmaPoints(xf_k, Pf_k, type);
            
            Pkk = zeros(n,n);
            for i = 1:length(W)
                Pkk = Pkk + (SP(:,i)-xf_k)*(f(SP(:,i),T)-xp_kplus1).' * W(i);
            end
            
            % smoothing gain
            Gk = Pkk*inv(Pp_kplus1);    
            
            xs = xf_k + Gk*(xs_kplus1 - xp_kplus1);
            
            Ps = Pf_k - Gk*(Pp_kplus1 - Ps_kplus1)*Gk';          
            
        case 'CKF'
            % calculate sigma points
            [SP,W] = sigmaPoints(xf_k, Pf_k, type);
            
            Pkk = zeros(n,n);
            for i = 1:length(W)
                Pkk = Pkk + (SP(:,i)-xf_k)*(f(SP(:,i),T)-xp_kplus1).' * W(i);
            end
            
            % smoothing gain
            Gk = Pkk*inv(Pp_kplus1);    
            
            xs = xf_k + Gk*(xs_kplus1 - xp_kplus1);
            
            Ps = Pf_k - Gk*(Pp_kplus1 - Ps_kplus1)*Gk';               
        otherwise
            error('Incorrect type of non-linear Kalman filter')  
    end
end