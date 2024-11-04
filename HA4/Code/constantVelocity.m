function fx = constantVelocity(x, T)
[n,N] = size(x);
fx = zeros(n,N);
for i = 1:N
fx(:,i) = [x(1,i)+T*x(3,i);x(2,i)+T*x(4,i);x(3,i);x(4,i)];
end
end