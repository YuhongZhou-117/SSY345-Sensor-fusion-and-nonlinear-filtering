function hx = cvMeasurement1(x)
N = size(x,2);
hx = zeros(2,N);
for i = 1:N
hx(:,i) = [x(3,i);x(4,i)]; 
end
end