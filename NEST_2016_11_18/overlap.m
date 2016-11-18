function [overlapping,area]=overlap(x1,x2,y1,y2,a1,a2,b1,b2,phi_shape1,phi_shape2)

% Overlap determines whether two ellipses centred on (x1,y1) and (x2,y2)
% and of major and minor half-axes (a,b), major axis rotated by an angle
% phi_shape to the x-axis overlap, and if they do, the area of overlap.

% units of area are determined from the units of the input dimensions, the
% default assumption is that all input sizes, positions are in m, so area
% in in m^2.

% Default assumption is that they don't overlap 
overlapping=false;
area=0.0;

% Separation
dx=x2-x1;
dy=y2-y1;

% Create a polygon for ellipse 1:
angle=linspace(0,2*pi,180);
x=a1*cos(angle);
y=b1*sin(angle);
crot=cos(phi_shape1);
srot=sin(phi_shape1);
x1test=crot*x-srot*y;
y1test=srot*x+crot*y;

% Create a polygon for ellipse 2:
x=a2*cos(angle);
y=b2*sin(angle);
crot=cos(phi_shape2);
srot=sin(phi_shape2);
x2test=crot*x-srot*y+dx;
y2test=srot*x+crot*y+dy;

% Do they overlap?
[in]=inpolygon(x1test,y1test,x2test,y2test);
% If there is an overlap then what's the area of overlap?
if numel(x1test(in)) > 1
    overlapping=true;
    area=polyarea(x1test(in),y1test(in));
end
end
