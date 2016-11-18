function [dxx, dyy, dzz, dxy, dxz, dyz] = magnetostatic_interaction_tensor_dipole(a1, b1, t1, v1, phi_a1, a2, b2, t2, v2, phi_a2, dx, dy, dz, npoints)

% Computes the 3x3 magnetostatic interaction tensor D where H_2=D*M_1 where:
% H2 is the interaction field experienced by the element number 2.
% M is the magnetisation of the source element number 1
% H and M are in in 3-D cartesian co-ordinates and SI units
% D is assumed symmetric so that Dxy=Dyx, Dxz=Dzx, Dyz=Dzy

% a,b,t are the semi-major and semi-minor axes and half-thickness for elements 1 and 2 respectively.
% v is the volume, assumed to be 2*t*pi*a*b
% phi_a is the angle between the semi-major axis and the x-direction, in the x-y plane
% The thickness is 2t in the z-direction, the cross-section and location of
% the element are invariant in the z-direction. 

% dx, dy, dz are the separation of the centres of the two elements (ie r2-r1) in cartesian co-ordinates and SI units 
% The two elements may touch but may not overlap, there is no overlap checking in this function.

% All quantities are assumed to be in SI units. It should be possible to
% scale all dimensions equally and to scale H and M equally, anyone wishing
% to do that should test it first. Note that using CGS units is not an
% equal scaling of M and H and some additional scaling wold be required. 

% The function performs a simple division of both the elements into
% approximately equal volume sub-elements and then performs a double
% integral of the dipole field over the volumes to obtain the interaction
% factors

% npoints is the number of sub-elements along the major axis of the element
% in the integral. Since the perimeter is ~2*pi*a then for an element with
% a~b~t there will be ~6*npoints^3 sub-elements. The factor calculation is
% a product of all pairs of sub-elements in both ellipsoids and so there will be 
% around 36*npoints^6 calculations required for each tensor component, about
% 200*npoints^6 in all. npoints~10 seems adequate, npoints<10 is likely to
% work well but accuracy shoud be checked (agains a larger value of
% npoints).

% Sub-divide element 1:
[xe1,ye1,ze1,dve1] = ellipse_elements(a1,b1,t1,v1,phi_a1,npoints);

% Sub-divide element 2:
[xe2,ye2,ze2,dve2] = ellipse_elements(a2,b2,t2,v2,phi_a2,npoints);

% Element 2 is at a position (dx,dy,dz) wrt element 1, which is at (0,0,0)
% Move the locations of the sub-elements of ellipse 2 to their positions
% relative to the element 1 origin.
xe2=xe2+dx;
ye2=ye2+dy;
ze2=ze2+dz;

% Compute the demag factors as a double integral over the two ellipses:
dxx=0.0;
dxy=0.0;
dxz=0.0;
dyy=0.0;
dyz=0.0;
dzz=0.0;

% Assume that each sub-element of ellipse 1 acts as a dipole source creating a
% field at each sub-element of ellipse 2. Weight each contribution by the 
% volume of each sub-element of ellipse 2 and ellipse 1.

% The inner loop (over all sub-elements of the field element (2)) is vectorised
for i1=1:size(xe1,2)
    dxe=xe2-xe1(i1);
    dxe2=dxe.*dxe;
    dye=ye2-ye1(i1);
    dye2=dye.*dye;
    dze=ze2-ze1(i1);
    dze2=dze.*dze;
    dre2=dxe2+dye2+dze2;
    dre3=dre2.*sqrt(dre2);
    mult3=dve1(i1)*dve2./(4.0*pi*dre3);
    mult5=3.0*mult3./dre2;
    dxx=dxx+sum(mult5.*dxe.*dxe-mult3);
    dxy=dxy+sum(mult5.*dxe.*dye);
    dxz=dxz+sum(mult5.*dxe.*dze);
    dyy=dyy+sum(mult5.*dye.*dye-mult3);
    dyz=dyz+sum(mult5.*dye.*dze);
    dzz=dzz+sum(mult5.*dze.*dze-mult3);
end

% The tensors are computed as a double sum over sub-elements of the two ellipses.
% Divide out the volume of ellipse 2 to get the field per unit volume. This
% means that the volume of ellipse 1 is included in D 
% Consequently H=DM where H and M are in A/m.

dxx=dxx/v2;
dxy=dxy/v2;
dxz=dxz/v2;
dyy=dyy/v2;
dyz=dyz/v2;
dzz=dzz/v2;

% all done
end
