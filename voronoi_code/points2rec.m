function [rectangleout] = points2rec(pA,pB,thickness)
%takes in two points and a thickness and returns a rectangle polyshape. The
%two points draw a line that cuts a line right down the middle of the
%rectangle and the thickness defines the width of the rectangle
%pA = [1,2];
%pB = [3,4];
RA = pB - pA; %vector from p1 to p2
RA_mag = (RA(1)^2+RA(2)^2)^0.5;
RB = thickness*0.5*[-RA(2)/RA_mag,RA(1)/RA_mag];

p1 = pA + RB;
p2 = pA - RB;
p3 = pB + RB;
p4 = pB - RB;

rectanglematrix = [p1;p3;p4;p2];

rectangleout = polyshape(rectanglematrix(:,1),rectanglematrix(:,2));

end