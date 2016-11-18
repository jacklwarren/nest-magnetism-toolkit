% Simple LLG integrator for macrospins.

%--------------------------------------------------------------------------
% S0: Introduction
%--------------------------------------------------------------------------

% To use this code read the introduction (this section), set the settable
% parameters in sections S1 - S6 and then run the main program (this).

% Note that the core code is quite flexible. In this illustrative example 
% there is limited flexibility by changing numerical parameters that are 
% separated out at the top of the code, greater flexibility is available
% by editing the code, which I have tried to self-document as well as I can.

% The core LLG integrator integrates M in spherical polar co-ordinates. For
% some applied field sequences (rotating fields etc) it's also useful to
% use polar co-ordinates. The code therefore requires some conversions.
% Co-ordinates are defined as in Matlab functions cart2sph and sph2cart:
% phi = azimuth = angle in the x-y plane from the +x axis towards the +y axis, in radians.
% theta = elevation = angle in degrees from the x-y plane towards +z axis, in radians.
% r = magnitude of the vector.
% so:
% (x,y,z) = (r.cos(theta).cos(phi), r.cos(theta).sin(phi), r.sin(theta))
% (phi,theta,r) = (atan2(y,x), asin(z/r), sqrt(x^2+y^2+z^2))
% (cart2sph and sph2cart should be used for co-ordinate transforms).

% Array indexing: (time,element,co-ordinate). This means that when
% interpolating time, time is the fastest moving index. At a fixed time
% element number is the fastest moving so that A(i,j) where i is element and j
% is co-ordinate number, if accessed linearly, is a list of co-ordinate i=1
% for each of the elements, then co-ordinate j for each element etc. 
% This should optimise array multiplications in the dM/dt routine.

clear all

global LLG_nelements
global LLG_alpha LLG_gamma
global LLG_Ms
global LLG_cryst LLG_K1_cart 
global LLG_shape LLG_nxx LLG_nyy LLG_nxy LLG_nzz
global LLG_magnx LLG_dxx LLG_dxy LLG_dxz LLG_dyy LLG_dyz LLG_dzz
global LLG_ex_up 
% add any global variables or constants that need to be shared between the 
% main code (this) and ANY functions here:

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Start of user-settable parameters, in 6 sections S1 - S6
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

disp('Initialising constants, computing pre-factors etc') 
start_time_cpu=cputime;
start_time_clock = clock;

%--------------------------------------------------------------------------
% S1: fundamental constants
%--------------------------------------------------------------------------

% Fundamental constants
mu_0=4*pi*1.0e-7;

% LLG constants:
gamma_g = -2.21e5; %gyromagnetic constant (negative)
alpha_g = 1.0; %damping term (positive)(should be ~0.2 for MRAM-like materials)
% calculate the constants for the LLG equation when reformulated into the LL form
% Note the signs of the constants used here, these are reflected in the code that computes dM/dt
% In the normal vector cross-product form of the equation, if the constants are positive then both of the terms MxH
% and MxMxH are negative. Here the constants are negative and so both terms
% of the equation are positive. 
LLG_gamma=gamma_g/(1+alpha_g^2); %precession constant in LL formulation of the LLG. If gamma_g is negative, so is gamma
LLG_alpha=alpha_g*LLG_gamma; %damping constant in LL formulation LLG. If gamma_g is negative, so is alpha
% Note that the 1/Ms term does not appear in alpha in the spherical polar form of the LLG

%--------------------------------------------------------------------------
% S2 structural information
%--------------------------------------------------------------------------

% Compute for a 10x10x2 array:
nelements_x=10;
nelements_y=10;
nelements_z=2;
LLG_nelements=nelements_x*nelements_y*nelements_z; %total number of magnetic elements

% Actual locations are generated later once the element dimesions are defined

% Shape of islands - this code is built for elliptical cross-section in the
% x-y plane, with vertical sides in the z direction, ie elliptic cylinders.
% Circular cross-section works as a special case of an ellipse (b=a). 

% The self-demag factors are computed using a function self_demagfactors.m, 
% created by Josephat Kalezhi, that computes self-demag factors for different
% shapes of elements. Specify that geometry:
% *DO NOT CHANGE*
islandgeo = 2; % Element geometry is elliptic cylinders 
alpha = 0.999; % ratio of top to bottom axes, only used for elliptic cones. 
% *For other shapes this code would need to be modified consistently throughout*
% *and in particular a function to compute magnetostatic and exchange* 
% *interactions between other shapes would be required.*

% Elliptical cylinder properties:
a_mean = 10.0e-9; % (m) mean semi-major axis length
a_sigma = 0.0e-9; % (m) sigma semi-major axis length
b_mean = 10.0e-9; % (m) mean semi-minor axis length
b_sigma = 0.0e-9; % (m) sigma semi-minor axis length
t_mean = 5.0e-9;  % (m) mean half-height (thickness)
t_sigma = 0.0e-9; % (m) sigma height (thickness)
% Allow the major axis to be at an angle phi to the x-axis
phi_shape_mean = 0.0; % (radians), mean angle between the major axis and the x-axis, from x anti-clockwise to the major axis
phi_shape_sigma = 0.0; % (radians), sigma of the angle between the major axis and the x-axis

LLG_shape=true;% is shape anisotropy (self-demag) to be included? (true/false)
demag_tol = 1e-6;  % tolerance used in the integrals of self-demag factors.

% Separations of elements - needed if there are interactions, otherwise unused.
% Note that this is the spacing of the element centres and that the closest
% is 2a in the x-direction, 2b in the y-direction and 2t in  the z-direction,
% ignoring any possible rotation of the element major axis.
x_spacing=a_mean*2.1;
y_spacing=b_mean*2.1;
z_spacing=t_mean*2.02;

% Specify whether periodic boundary conditions are required along each axis
wrap_x=true;
wrap_y=true;
wrap_z=false;

%--------------------------------------------------------------------------
% S3 Magnetostatic parameters
%--------------------------------------------------------------------------

Ms_mean=1.0e6; % Mean value of Ms in A/m
Ms_sigma=0.0e6; % Standard deviation of Ms in A/m

LLG_magnx=true;% are magnetostatic interactions included? (true/false)
npoints=10; % accuracy of tensor calculation if there are magnetostatic interactions 

%--------------------------------------------------------------------------
% S4 Inter-layer exchange parameters
%--------------------------------------------------------------------------
% Hexg=M*Cexg*M;
% For surface to surface separation d then:
% Cexg={C0+C1*exp(-d/D1)+C2*exp(-2*d/D2)}*cos{(2*pi*d/lambda)+phi}
% C0, C1 and C2 are all per unit area.
LLG_ex_up=true;
exg_terms=[true,false,false,false]; % include the following terms:
exg_exp_0=0.1/(pi*a_mean*b_mean); % (/m^2) C0 in Hexg=C0*M
exg_exp_1=[0.0,0.1e-9]; % [/m^2,m]: [C1,D1] in Hexg=C1*exp(-d/D1)*M
exg_exp_2=[0.0,0.1e-9]; % [/m^2,m]: [C2,D2] in Hexg=C2*exp(-2*d/D2)*M
exg_cos=[1.0e-9,0.0]; % [m,radians]: [lambda,phi] in cos{(2*pi*d/lambda)+phi}
% Note that to use the cosine term at least one of C0, C1 or C2 must be used.
% The direction of the cosine can be inverted by the sign of lambda and 
% the cosine term can be made into a sin by use of the phase angle phi.

% Setting exg_terms=[true,false,false,false] means that Cexg = C0 which can
% be used for exchange constant calculated by any external method.

% As yet there is no variability of exchange coupling between elements,
% this could be included by following a similar process to other variable
% constants

%--------------------------------------------------------------------------
% S5 Crystalline anisotropy
%--------------------------------------------------------------------------
% K1 (vector in spherical polar co-ordinates):
LLG_cryst=true;% is crystalline anisotropy to be included? (true/false)
K1_mean_phi=0.0; % radians
K1_sigma_phi=0.0; % radians
%K1_sigma_phi=0.0873; % radians
K1_mean_theta=0.0; % radians
K1_sigma_theta=0.0; % radians
%K1_sigma_theta=0.0873; % radians
K1_mean_r=1.0e6; % Warning: Setting this to zero would cause a serious problem in scaling the anisotropy field, use a realistic very low value. 
K1_sigma_r=0.0e6; % J/m^3

%--------------------------------------------------------------------------
% S6 Applied field
%--------------------------------------------------------------------------
% Applied field sequence
ntimes=5000; % number of times at which H is defined and at which a value of M is returned by the ODE solver 
loop_half_time=1e-9; % the time taken for half a loop
% Vector t_H defines those times 
t_H=linspace(0,5,ntimes)*loop_half_time;
H_app_max=2.0*K1_mean_r/(mu_0*Ms_mean);% maximum magnitude of the applied field
% In this example multiple loops are computed for fields at different
% angles, so the field initialisation has to be done in the main code. For
% simpler cases there could be more user settable parameters.

%--------------------------------------------------------------------------

% If you create additional effective field terms then add settable parameters here.
% Where scalar properties differ between elements create a single column matrix (LLG_nelements,1)
% for each property.

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% End of user-settable parameters
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Start initialisation of parameters, factors etc
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% Create the locations of the elements

% To locate on a 3-D cubic grid:
[x,y,z] = meshgrid((0:nelements_x-1)*x_spacing,(0:nelements_y-1)*y_spacing,(0:nelements_z-1)*z_spacing);

x=reshape(x,[LLG_nelements,1]);
y=reshape(y,[LLG_nelements,1]);
z=reshape(z,[LLG_nelements,1]);
% Identify which layer they are on, counting from 0:
layer=round(z/z_spacing); %this can be used to set different parameters for different layers

% Set up the values of Ms
LLG_Ms(1:LLG_nelements,1)=Ms_mean+Ms_sigma*randn(size(1:LLG_nelements)); %A/m

%--------------------------------------------------------------------------
% Crystalline anisotropy:
%--------------------------------------------------------------------------
if LLG_cryst
display('Calculating crystalline anistotropy factors');
% Initialise crystaline anisotropy
% K1 in polar co-ordinates:
K1_sph(1:LLG_nelements,1)=K1_mean_phi+K1_sigma_phi*randn(size(1:LLG_nelements));
K1_sph(1:LLG_nelements,2)=K1_mean_theta+K1_sigma_theta*randn(size(1:LLG_nelements));
K1_sph(1:LLG_nelements,3)=K1_mean_r+K1_sigma_r*randn(size(1:LLG_nelements));
Hk=(2.0*K1_sph(:,3))./(mu_0*LLG_Ms);
Hk_cryst=mean(Hk);
% Scale the magnitude of K1 so that H=(K.m)K where K is the scaled
% anisotropy vector and m is the unit magnetisation vector.
Hanis_r=2.0./(mu_0.*K1_sph(:,3).*LLG_Ms.*LLG_Ms);
% Multiply K_cart by sqrt(Hanis_r)because in the LLG routine we calculate
% (K.M)K and we need to embed the multiplicative constant H_anis_r. Do that
% while converting K1 to cartesian co-ordiantes:
[LLG_K1_cart(:,1),LLG_K1_cart(:,2),LLG_K1_cart(:,3)]=sph2cart(K1_sph(:,1),K1_sph(:,2),K1_sph(:,3).*sqrt(Hanis_r));
end

%--------------------------------------------------------------------------
% Shape anisotropy:
%--------------------------------------------------------------------------
if LLG_shape
display('Calculating shape anistotropy factors');
% Initialise shape anisotropy demag factors (N):
a(1:LLG_nelements,1)=a_mean+a_sigma*randn(size(1:LLG_nelements));
b(1:LLG_nelements,1)=b_mean+b_sigma*randn(size(1:LLG_nelements));
t(1:LLG_nelements,1)=t_mean+t_sigma*randn(size(1:LLG_nelements));
v(1:LLG_nelements,1)=2*t.*pi.*a.*b;
[Nx Ny Nz] = self_demagfactors(a, b, t, alpha, islandgeo, demag_tol);
Nx=Nx';
Ny=Ny';
Nz=Nz';

Ndiff=Nx-Ny; % (should be Ny-Nx) but negate here so that the field calculated in the dM/dt routine is in the correct direction

% self_demagfactors.m computes nxx, nyy assuming that the major axis is in the x-direction. If not then rotate the axes:
%...but in the x-y plane the major axis of the element is at an angle to the x-axis, so need to transform the x-y self-demag factors:
phi_shape(1:LLG_nelements,1)=phi_shape_mean+phi_shape_sigma*randn(size(1:LLG_nelements));
sin_phi_shape=sin(phi_shape);
cos_phi_shape=cos(phi_shape);
LLG_nxx = Ndiff.*sin_phi_shape.*sin_phi_shape - Nx;
LLG_nxy = -Ndiff.*sin_phi_shape.*cos_phi_shape;
LLG_nyy = Ndiff.*cos_phi_shape.*cos_phi_shape - Nx;
LLG_nzz = -Nz;

K1_shape=0.5*mu_0*LLG_Ms.*LLG_Ms.*abs(min(Ny,Nz)-Nx); %Factor half as it's a self-energy
% in 2-D K1_shape would be ~(Ny-Nx) but remember that M could rotate through the z-axis instead of the y-axis, so choose the easiest one (smallest of Ny, Nz).
Hk_shape=mean((2.0*K1_shape)./(mu_0*LLG_Ms));

Hk=Hk_cryst+Hk_shape;
end

%--------------------------------------------------------------------------
% If the moments interact then we need to specify element positions
%--------------------------------------------------------------------------
if LLG_magnx || LLG_ex_up
    
try % try to load pre-calculated interactions...
    
pre_warn=num2str(now); %set the last warning message to something unlikely
lastwarn(pre_warn); 
disp('Attempting to load interaction factors from file interactions.mat')
load interactions.mat LLG_dxx LLG_dxy LLG_dxz LLG_dyy LLG_dyz LLG_dzz
% See if any warnings were issued, if so then assert an error to force recompute
assert(strcmp(lastwarn,pre_warn),'Failed to load interaction factors, computing...')
% Check that the array size is correct, if not then assert an error to force recompute
% If x failed to load that will have already caused an error, so only check size if variable x exists
if exist('LLG_dxx','var')
assert(size(LLG_dxx,1)==LLG_nelements,'Size of factor arrays loaded does not match problem size')
end

catch % An error ocurred when loading so recompute the interaction factors.
    
% Say what's know about why they failed to load:
disp(lasterr)

% Implement periodic boundary conditions - assume that the world wraps
% around in specified dimensions.
% If periodic boundary conditions are required for axis k (k=x.y. or z) 
% then set k_wrap to the positive size of the system in the k-direction. 
if wrap_x
    x_wrap=nelements_x*x_spacing;
end
if wrap_y
    y_wrap=nelements_y*y_spacing;
end
if wrap_z
    z_wrap=nelements_z*z_spacing;
end

%--------------------------------------------------------------------------
% Magnetostatic interactions
%--------------------------------------------------------------------------
if LLG_magnx
display('Calculating magnetostatic interaction factors');
% Field H is computed as Hx=Dxx*Mx etc.
% M is a column matrix of size (LLG_nelements,1)
% D is a square matrix of size (LLG_nelements,LLG_nelements)
% H is a column matrix of size (LLG_nelements,1)

% Hence H(i)=D(i,j).M(j) computes the contribution to the field experienced
% by element i arising from the magnetisation of element j

% Symmetry: D is defined such that H and M are in A/m ie per volume. 
% The field should be proportional to the volume of the source.
% If the two volumes are the same then Dji=Dij.
% If the volumes are different then Dji=Dij*vi/vj

   for j=1:LLG_nelements
       LLG_dxx(j,j)=0.0; % no self-demag, only element-element interactions
       for i=j+1:LLG_nelements % use symmetry, only compute for j>i
           % find the separation from j to i: dr_ij=(dx,dy,dz)
           % wrap the separation if required, signalled by a positive wrap
           % size, separately selected for each co-ordinate
           dx=x(i,1)-x(j,1);% Distance from source(j) to field (i)
           if wrap_x && (abs(dx) > x_wrap/2.0)
               dx=-1*sign(dx)*(x_wrap-abs(dx));
           end
           
           dy=y(i,1)-y(j,1);% Distance from source(j) to field (i)
           if wrap_y && (abs(dy) > y_wrap/2.0)
               dy=-1*sign(dy)*(y_wrap-abs(dy));
           end
           
           dz=z(i,1)-z(j,1);% Distance from source(j) to field (i)
           if wrap_z && (abs(dz) > z_wrap/2.0)
               dz=-1*sign(dz)*(z_wrap-abs(dz));
           end
           
           % Compute the magnetostatic interaction factors
           [LLG_dxx(i,j),LLG_dyy(i,j),LLG_dzz(i,j),LLG_dxy(i,j),LLG_dxz(i,j),LLG_dyz(i,j)] = magnetostatic_interaction_tensor_dipole(a(j,1), b(j,1), t(j,1), v(j,1), phi_shape(j,1), a(i,1), b(i,1), t(i,1), v(i,1), phi_shape(i,1), dx, dy, dz, npoints);
           
           % Fill the full matrix using symmetry:
           LLG_dxx(j,i)=LLG_dxx(i,j)*v(i)/v(j);
           LLG_dxy(j,i)=LLG_dxy(i,j)*v(i)/v(j);
           LLG_dxz(j,i)=LLG_dxz(i,j)*v(i)/v(j);
           LLG_dyy(j,i)=LLG_dyy(i,j)*v(i)/v(j);
           LLG_dyz(j,i)=LLG_dyz(i,j)*v(i)/v(j);
           LLG_dzz(j,i)=LLG_dzz(i,j)*v(i)/v(j);
       end
   end

end

%--------------------------------------------------------------------------
% Inter-layer exchange 
%--------------------------------------------------------------------------

if LLG_ex_up
display('Calculating interlayer exchange interaction factors');
   if ~LLG_magnx % if the interaction tensors haven't been set up yet, 
                 % create them as zero arrays of the right size:
       LLG_dxx=zeros(LLG_nelements);
       LLG_dxy=zeros(LLG_nelements);
       LLG_dxz=zeros(LLG_nelements);
       LLG_dyy=zeros(LLG_nelements);
       LLG_dyz=zeros(LLG_nelements);
       LLG_dzz=zeros(LLG_nelements);
   end
   % Then add the exchange constant C to the interaction tensor D:
   noverlap=0;
   for i=1:LLG_nelements
       for j=i+1:LLG_nelements % use symmetry, only compute for j>i
           [overlapping,area]=overlap(x(i,1),x(j,1),y(i,1),y(j,1),a(i,1),a(j,1),b(i,1),b(j,1),phi_shape(i,1),phi_shape(j,1));
           if overlapping
               noverlap=noverlap+1;
               dz=abs(z(j,1)-z(i,1))-(t(i)+t(j));% separation of the top surface of element(j) and the bottom of element (i)
               Cexg=interlayer_exchange(dz,exg_terms,exg_exp_0,exg_exp_1,exg_exp_2,exg_cos);
               % interlayer_exchange returns the exchange constant per unit area 
               Cexg_a=Cexg*area;
               % The calculated Cexg_a is used in H=Cexg_a*M where H and M are
               % vectors and Cexg_a is a scalar, so the implementation is
               % just to add Cexg_a to the leading diagonal of the
               % interaction matrix Dij (and since the overlap area and 
               % separation are the same for ji as for ij then do both 
               % at the same time) 
               LLG_dxx(i,j)=LLG_dxx(i,j)+Cexg_a;
               LLG_dxx(j,i)=LLG_dxx(j,i)+Cexg_a;
               LLG_dyy(i,j)=LLG_dyy(i,j)+Cexg_a;
               LLG_dyy(j,i)=LLG_dyy(j,i)+Cexg_a;
               LLG_dzz(i,j)=LLG_dzz(i,j)+Cexg_a;
               LLG_dzz(j,i)=LLG_dzz(j,i)+Cexg_a;
           end
       end
   end       
end

disp('Interactions computed, saving for future re-use')
save interactions.mat LLG_dxx LLG_dxy LLG_dxz LLG_dyy LLG_dyz LLG_dzz

end

end
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% End initialisation of parameters, factors etc
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

end_time_cpu=cputime;
cpu_time_setup=end_time_cpu-start_time_cpu;
disp(['Setup completed in ',num2str(cpu_time_setup),'seconds'])

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Compute loops
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%--------------------------------------------------------------------------
% Open and set up figures to plot results
%--------------------------------------------------------------------------
figure(3)
clf
hold on
grid on
g1=hggroup;
g2=hggroup;
g3=hggroup;
title('Plot of M as a function of field');
xlabel('Happlied(t)/Hk'); 
ylabel('M(t)/Ms');

figure(4)
clf
hold on
grid on
title('Coercivity as a function of angle');
xlabel('Angle of field to anisotropy axis (degrees)'); 
ylabel('Coercivity');

%--------------------------------------------------------------------------
% Compute a loop for each angle
%--------------------------------------------------------------------------
for iangle=1:19

start_time_cpu=cputime;
start_time_clock = clock;

angle(iangle)=(iangle-1)*5*pi/180.0; % angle in radians

% Avoid having the applied field exactly along the easy axis, move it slightly:
if iangle == 1
    angle(iangle)=0.5*pi/180.0;
elseif iangle == 19
    angle(iangle)=89.5*pi/180.0;
end
display(['Computing for angle ',num2str(angle(iangle)*180/pi)]);

%--------------------------------------------------------------------------
% Setup the initial magnetisation and H(t) for each loop
%--------------------------------------------------------------------------

% Initialise the magnetisation:
% %M_init is the initial value of M, M_init(:,1:3) phi,theta and |M| components of each of the LLG_nelements.
M_init(1:LLG_nelements,3)=LLG_Ms;
M_init(1:LLG_nelements,1)=angle(iangle)*ones(size(LLG_Ms));
M_init(1:LLG_nelements,2)=zeros(size(LLG_Ms));
% Having set up M in spherical polar, transform to cartesian to integrate:
% M_init_cart is the initial value of M, M_init_cart(:,1:3) are the x,y and z components of each of the LLG_nelements.
[M_init_cart(:,1),M_init_cart(:,2),M_init_cart(:,3)]=sph2cart(M_init(:,1),M_init(:,2),M_init(:,3));

% Define applied field as a function of time:
% Vector t_H defines the times at which the value of the applied field H is defined
% Specify the corresponding H values in spherical polar form. 
% Split the time into 4 equal sections:
% a) Field in direction angle(iangle), magnitude H_app_max (constant)
% b) Field in direction angle(iangle), ramping linearly from H_app_max to 0
% c) Field in direction -angle(iangle), ramping linearly from 0 to H_app_max 
% d) Field in direction -angle(iangle), magnitude H_app_max (constant)
Happlied_t=zeros(ntimes,LLG_nelements,3); % create an array of the right size.
Happlied_t(1:ntimes/2,:,1)=angle(iangle); % Sets H_applied to be in the +angle direction for the first part of the loop
Happlied_t(ntimes/2+1:ntimes,:,1)=pi+angle(iangle); % Sets H_applied to be in the -angle direction for the first part of the loop:
Happlied_t(:,:,2)=0; % H_applied is in the x-y plane at all times
Happlied_t(1:ntimes/4,:,3)=H_app_max*ones(ntimes/4,LLG_nelements); % section a) |H|=H_app_max
Happlied_t(3*ntimes/4+1:ntimes,:,3)=H_app_max*ones(ntimes/4,LLG_nelements);% section d) |H|=H_app_max
H_rampdown=linspace(H_app_max,0.0,ntimes/4);% create linear ramp up
H_rampup=linspace(0.0,H_app_max,ntimes/4);% create linear ramp down
for ielement=1:LLG_nelements
Happlied_t(ntimes/4+1:ntimes/2,ielement,3)=H_rampdown; % for each element, section b) ramps magnitude down 
Happlied_t(ntimes/2+1:3*ntimes/4,ielement,3)=H_rampup; % for each element, section c) ramps magnitude up
end

%--------------------------------------------------------------------------
% Everything is set up, now compute the loop
%--------------------------------------------------------------------------
% Solve the LLG to calculate M(t)
[t_solved, M_solved] = ode113(@(t,M) dM_dt(t,M,t_H,Happlied_t),t_H,M_init(:,1:2)); % Solve ODE

end_time_cpu=cputime;
end_time_clock=clock;
cpu_time_elapsed(iangle)=end_time_cpu-start_time_cpu;
clock_time_elapsed(iangle)=etime(end_time_clock,start_time_clock);
disp(['Solved in ',num2str(cpu_time_elapsed(iangle)),'seconds'])

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Loop computed
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Post process the loop: compute coercivity, draw graphs etc
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% The ODE functions use single row vecors, reshape the result array into 
% the array structure described in the header:
M_solved=reshape(M_solved,ntimes,LLG_nelements,2);
% recreate the magnitude of M for the solution (ODE assumes |M| = Ms = constant).
M_solved(:,:,3)=repmat(LLG_Ms',ntimes,1);

% Convert the resulting M(t) from polar to cartesian form:
[M_x,M_y,M_z]=sph2cart(M_solved(:,:,1),M_solved(:,:,2),M_solved(:,:,3));

% Compute an overall average of M/Ms for the assembly:
if LLG_nelements < 2 % if there is only one element then it's easy:
    M_x_tot=M_x(:,1)'/LLG_Ms(1);
    M_y_tot=M_y(:,1)'/LLG_Ms(1);
    M_z_tot=M_z(:,1)'/LLG_Ms(1);
else % otherwise do a weighted sum according to volume:
    M_x_tot=(M_x*v(:,1))/sum(LLG_Ms.*v(:,1));
    M_y_tot=(M_y*v(:,1))/sum(LLG_Ms.*v(:,1));
    M_z_tot=(M_z*v(:,1))/sum(LLG_Ms.*v(:,1));
end

% Convert the applied field at the sample times from polar to cartesian in case required
[Happ_x,Happ_y,Happ_z]=sph2cart(Happlied_t(:,:,1),Happlied_t(:,:,2),Happlied_t(:,:,3));

% In case the results are needed for other post-processing, save them:
save M_vs_t M_x M_y M_z M_x_tot M_x_tot M_x_tot Happ_x Happ_y Happ_z t_solved x y z

% The applied field at the applied field angle is just the magnitude of the applied field, the sign is derived from the in-plane angle:
H_angle=Happlied_t(:,1,3).*sign(cos(Happlied_t(:,1,1)));

% If you want to find M  in the direction of H:
% Create a unit vector in the applied field direction:
% [H_angle_unit_x,H_angle_unit_y,H_angle_unit_z]=sph2cart(Happlied_t(1,1,1),Happlied_t(1,1,2),1.0);
% Then find the projection of the magnetisation onto the field direction (dot-product):
% M_angle=(M_x_tot*H_angle_unit_x+M_y_tot*H_angle_unit_y+M_z_tot*H_angle_unit_z);

% To find the coercivity map H and M onto the easy axis direction
% Create a unit vector in the easy axis direction:
[K1_angle_unit_x,K1_angle_unit_y,K1_angle_unit_z]=sph2cart(K1_mean_phi,K1_mean_theta,1.0);
% Then find the projection of the magnetisation onto the easy axis direction (dot-product):
M_angle=(M_x_tot*K1_angle_unit_x+M_y_tot*K1_angle_unit_y+M_z_tot*K1_angle_unit_z);
% M_angle=M_x_tot; (simpler way to do it if the easy axis is in the x-direction)

% and (if you wanted it) similarly for H:
% H_angle=(Happ_x*K1_angle_unit_x+Happ_y*K1_angle_unit_y+Happ_z*K1_angle_unit_z);

% find the point where the magnetisation along the easy-axis switches direction
iz=1;
while (~((M_angle(iz+1) < 0) && (M_angle(iz) >= 0))) && (iz < ntimes-5)
iz=iz+1;
end
M_sample=M_angle(iz-5:iz+5);
H_sample=H_angle(iz-5:iz+5);
Hc(iangle)=abs(interp1(M_sample,H_sample,0.0));
% and that's the coercivity

% % Plot the solution M(t) as a function of time:
% figure(1)
% clf
% hold on
% grid on
% g1=hggroup;
% g2=hggroup;
% g3=hggroup;
% plot(t_solved, M_x,'r-o','Parent',g1);
% plot(t_solved, M_y,'g-o','Parent',g2);
% plot(t_solved, M_z,'b-o','Parent',g3);
% title('Plot of M as a function of time');
% xlabel('Time'); 
% ylabel('M(t) (A/m)');
% legend([g1,g2,g3],'Mx','My','Mz','Location','southeast');
% 
% figure(2)
% clf
% hold on
% grid on
% g1=hggroup;
% g2=hggroup;
% g3=hggroup;
% plot(t_H,Happ_x(:,1),'r-o','Parent',g1);
% plot(t_H,Happ_y(:,1),'g-o','Parent',g2);
% plot(t_H,Happ_z(:,1),'b-o','Parent',g3);
% title('Plot of Happlied as a function of time');
% xlabel('Time'); 
% ylabel('Happlied(t) (A/m)');
% legend([g1,g2,g3],'Hx','Hy','Hz','Location','southeast');

% Plot the solution M(t) as a function of field, both along the easy axis:
figure(3)
%plot(H_angle(ntimes/2:ntimes,1)/Hk,M_angle(ntimes/2:ntimes),'r-o','Parent',g1);
plot(H_angle,M_angle,'r-o','Parent',g1);
plot(0.0-Hc(iangle),0.0,'bx')
% plot(H_angle(:,1)/Hk,M_y_tot,'g-o','Parent',g2);
% plot(H_angle(:,1)/Hk,M_z_tot,'b-o','Parent',g3);

end

figure(3)
legend([g1,g2,g3],'Mx/Ms','My/Ms','Mz/Ms','Location','southeast');

figure(4)
plot(angle*180/pi,Hc,'r-+')



