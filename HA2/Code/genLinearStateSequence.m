function X = genLinearStateSequence(x_0, P_0, A, Q, N)
%GENLINEARSTATESEQUENCE generates an N-long sequence of states using a 
%    Gaussian prior and a linear Gaussian process model
%
%Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   N           [1 x 1] Number of states to generate
%
%Output:
%   X           [n x N+1] State vector sequence
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [];
% x0
X(:,1) = mvnrnd(x_0,P_0)';
for i = 2 : N+1
   X(:,i) = A*X(:,i-1)+mvnrnd(zeros(1,length(Q)),Q)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end