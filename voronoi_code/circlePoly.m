function [polyout] = circlePoly(point,radius)
cancel = false;
if (abs(point(1))==Inf)||(abs(point(2))==Inf)
    cancel = true;
end

    
if ~cancel
numpoints = 12;
theta = 2*3.14159/numpoints;
pointsmatrix = zeros(numpoints,2);
for i=1:numpoints
    
    y = radius*sin(theta*i);
    x = radius*cos(theta*i);
    pointsmatrix(i,1) = x + point(1);
    pointsmatrix(i,2) = y + point(2);
end

polyout = polyshape(pointsmatrix(:,1),pointsmatrix(:,2));
else
    polyout = polyshape([0,0],[0,0]);
end
end