function [xe,ye,ze,dve] = ellipse_elements(a,b,t,v,phi_a,npoints)

% ellipse_elements subdivides an elliptical cross-section column into
% roughly equal volume elements distributed roughly evenly over the volume.

% a: major axis half-length (equivalent to radius for a circular cross section)
% b: minor axis half-length 
% t: half-height (t is defined equivalently to a and b)
% v is the volume, assumed to be 2*t*pi*a*b
% phi_a: (radians) the angle of the major axis anti-clockwise from the x-axis
% npoints: the number of elements along the radius. The routine subdivides
% the volume in approximately equal intervals radially, vertically and
% around each radial strip. Therefore for a=b=t there are approximately
% 2*pi*npoints^3, ie for npoints = 10 there are about 6000 elements.  

% Set up the radius values:
na=npoints;
da=a/na;
ai=((1:na)-0.5)*da; % the first point will be in the centre but add that later to avoid numerical problems
bi=ai*b/a;

%then the phi values:
% calculate perimeter by method 3 of http://www.mathsisfun.com/geometry/ellipse-perimeter.html
hi=((ai-bi).*(ai-bi))./((ai+bi).*(ai+bi));
perim_i=pi*(ai+bi).*(1+3*hi./(10+sqrt(4-3*hi)));
% the points should be spaced around the elipse at about the same spacing as they are radially 
nphi_i=max(8,round(npoints*perim_i/a));
dphi_i=2*pi./nphi_i;

% aim for the same point spacing through the thickness as radially:
nz=max(3,round(t/da));
dz=2.0*t/nz;
zearray=(((1:nz)-0.5)*dz) - t;

% Create the elements

iaphi=0; % the number of (a,phi) values for which we have already created nz elements

for ia=1:na
    clear phi x y xphi yphi dxp dyp xr yr dxr dyr darea

% generate evenly spaced points around the circumference 
    phi=(1:2*nphi_i(ia))*dphi_i(ia)/2;
% repeat the first and last points after and before to enable differencing
    phi=[phi(2*nphi_i(ia)) phi phi(1)];
    
% co-ordinates of the points on the perimter:
    x=ai(ia)*cos(phi);
    y=bi(ia)*sin(phi);
% rotate the ellipse by phi_a:
    crot=cos(phi_a);
    srot=sin(phi_a);
    xphi=crot*x-srot*y;
    yphi=srot*x+crot*y;

% compute the x and y displacements along the perimeter for this sub-element
% using the locations of the sub-elements to either side
    dxp=xphi(3:2:2*nphi_i(ia)+2);
    dxm=xphi(1:2:2*nphi_i(ia));
    dx=dxp-dxm;
    dyp=yphi(3:2:2*nphi_i(ia)+2);
    dym=yphi(1:2:2*nphi_i(ia));
    dy=dyp-dym;

% compute the x and y displacements of the radial vector. The location of
% the sub-element is:
    xr=xphi(2:2:2*nphi_i(ia)+1);
    yr=yphi(2:2:2*nphi_i(ia)+1);
% the radial vector points along this direction, so (xr,yr) needs scaling
% down. At each angle the step in radius is 1/(ia-0.5)*radius
	rscale=1/(ia-0.5);
    dxr=xr*rscale;
    dyr=yr*rscale;
    
% sub-element area is the cross-product of the radial and circumferential vectors:
    darea=dxr.*dy-dyr.*dx;

% the sub-element locations and areas are calculated for this
% circumferential strip, now generate the points for all the sub-elements
% up through the thickness for each angle at this radius:
    for iphi=1:nphi_i(ia)
% Calculate the start and end points in the result array of all sub-elements
% of this column of nz elements at (ai(ia),phi(iphi))  
        istart=iaphi*nz+1;
        iaphi=iaphi+1;
        iend=iaphi*nz;
% Then fill this block of the result array
        xe(istart:iend)=xr(iphi);
        ye(istart:iend)=yr(iphi);
        ze(istart:iend)=zearray;
        dve(istart:iend)=darea(iphi)*dz;
    end % end iphi loop
end % end ia loop

% The total volume of the sub-elements will not quite exactly sum to the
% correct volume, correct this so that the total volume of the elements is
% the value specified by the calling program. 
% This should be a very small correction, probably never really necessary.
scale=v/sum(dve);
dve=dve*scale;

% all done.
end
