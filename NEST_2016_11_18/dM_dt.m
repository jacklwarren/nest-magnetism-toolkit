function dMdt = dM_dt(t,M,t_H,Happlied_t)

% Computes dM/dt in polar co-ordinates using the LL form of the LLG.
% M, t are inputs from the calling routine. The routine returns dM/dt(t,M).
% Happlied_t is the applied field at the times in the array t_H and is
% passed in from the main routine. 
% alpha, gamma are the constants in the LL formulation of the LLG equation.
% All physical vectors are passed in in polar co-ordinates (see below).

% Note that H_applied is specified in polar co-ordinates, so it is
% interpolated in polar form, then converted to cartesian to add to the
% other firled terms, then converted back to polar form. This is
% inefficient computationally but allows rotating fields to be specified accurately in
% a small number of inputs. That may not be the most efficient compromise
% in some cases.

% co-ordinates: defined as in Matlab functions cart2sph and sph2cart.
% phi = azimuth = angle in the x-y plane from the +x axis towards the +y axis, in radians.
% theta = elevation = angle in degrees from the x-y plane towards +z axis, in radians.
% r = magnitude of the vector.
% so:
% (x,y,z) = (r.cos(theta).cos(phi), r.cos(theta).sin(phi), r.sin(theta))
% (phi,theta,r) = (atan2(y,x), asin(z/r), sqrt(x^2+y^2+z^2))
% but cart2sph and sph2cart should be used for co-ordinate transforms.

% Array indexing: (time,element,co-ordinate). This means that when
% interpolating time, time is the fastest moving index. At a fixed time the
% element number is the fastest moving index so that A(i,j) where i is element and j
% is co-ordinate number, if accessed linearly, is a list of co-ordinate i=1
% for each of the elements, then co-ordinate j etc. This should optimise
% array multiplications in the dM/dt routine.

global LLG_nelements
global LLG_alpha LLG_gamma
global LLG_Ms
global LLG_cryst LLG_K1_cart 
global LLG_shape LLG_nxx LLG_nyy LLG_nxy LLG_nzz
global LLG_magnx LLG_dxx LLG_dxy LLG_dxz LLG_dyy LLG_dyz LLG_dzz
global LLG_ex_up 

% add any further global variables or constants that are needed in the dM/dt routine here:

% ODE routines provide the time dependent variable as a row vector.
% Separate out the two co-ordinates phi and theta:
M_polar=reshape(M,LLG_nelements,2);
% Calculate the cartesian coordinates of the UNIT magnetisation vector m:
[M_cart(:,1),M_cart(:,2),M_cart(:,3)]=sph2cart(M_polar(:,1),M_polar(:,2),LLG_Ms);

% Interpolate the supplied field values (t_H,Happlied_t) to find Happ at time t
Happ_phi = interp1(t_H,Happlied_t(:,:,1),t);
Happ_theta = interp1(t_H,Happlied_t(:,:,2),t);
Happ_r = interp1(t_H,Happlied_t(:,:,3),t);
% Convert the applied field from spherical polar to cartesian co-ordinates
[Happ_x,Happ_y,Happ_z]=sph2cart(Happ_phi,Happ_theta,Happ_r);

% add any effective field terms, calculating [Hx,Hy,Hz] for all elements here:

if LLG_cryst
% Anisotropy field (NOTE THAT LLG_K1_cart IS PRE_SCALED to include the multiplicative 
% factor that is needed to compute field in A/m):
    Hanis_dot=dot(LLG_K1_cart,M_cart,2);
    Hanis_x=Hanis_dot.*LLG_K1_cart(:,1);
    Hanis_y=Hanis_dot.*LLG_K1_cart(:,2);
    Hanis_z=Hanis_dot.*LLG_K1_cart(:,3);
else
    Hanis_x=zeros(size(Happ_x));
    Hanis_y=zeros(size(Happ_x));
    Hanis_z=zeros(size(Happ_x));
end
    
if LLG_shape
    Hshape_x = LLG_nxx.*M_cart(:,1)+LLG_nxy.*M_cart(:,2);
    Hshape_y = LLG_nxy.*M_cart(:,1)+LLG_nyy.*M_cart(:,2);
    Hshape_z = LLG_nzz.*M_cart(:,3);
else
    Hshape_x=zeros(size(Happ_x));
    Hshape_y=zeros(size(Happ_x));
    Hshape_z=zeros(size(Happ_x));
end

if LLG_magnx || LLG_ex_up

% Interaction field computed by a full-matrix N*N method. The code stores
% the full N*N interaction matrix. The storage and thus the
% memory-processor data transfer time could be reduced by only storing the
% upper or lower triangle of the matrix D at the expense of time lost in
% translating indices. The full N*N matrix is stored assuming that Matlab
% code will only be used for small memory problems where the size of D will
% be small and D will probably be cached. For larger problems different
% optimisation is sensible.

% Magnetostatic interaction field H is computed as Hx=Dxx*Mx etc.
% M is a column matrix of size (LLG_nelements,1)
% D is a square matrix of size (LLG_nelements,LLG_nelements)
% H is a column matrix of size (LLG_nelements,1)

% Exchange interaction field H is computed as H=C*M where H and M are 
% vectors and C is a scalar. This is implemented by adding C to the leading
% diagonal of the magnetostatic interaction matrix D during the setup phase. 

% Hence H(i)=D(i,j).M(j) computes the contribution to the field experienced
% by element i arising from the magnetisation of element j

    HMint_x = LLG_dxx*M_cart(:,1)+LLG_dxy*M_cart(:,2)+LLG_dxz*M_cart(:,3);
    HMint_y = LLG_dxy*M_cart(:,1)+LLG_dyy*M_cart(:,2)+LLG_dyz*M_cart(:,3);
    HMint_z = LLG_dxz*M_cart(:,1)+LLG_dyz*M_cart(:,2)+LLG_dzz*M_cart(:,3);
else
    HMint_x=zeros(size(Happ_x));
    HMint_y=zeros(size(Happ_x));
    HMint_z=zeros(size(Happ_x));
end

%sum all effective field terms
H_x = Happ_x'+Hanis_x+Hshape_x+HMint_x; %+...any other effective field terms
H_y = Happ_y'+Hanis_y+Hshape_y+HMint_y; %+...any other effective field terms
H_z = Happ_z'+Hanis_z+Hshape_z+HMint_z; %+...any other effective field terms

% Calculate dM/dt:
% Calculate sin and cos of the magnetisation angles:
cph_M=cos(M_polar(:,1)); %cos(phi_M)
sph_M=sin(M_polar(:,1)); %sin(phi_M)
cth_M=cos(M_polar(:,2)); %cos(theta_M)
sth_M=sin(M_polar(:,2)); %sin(theta_M)

% find vector H expressed in the polar co-ordinates of vector M: [H_M_phi, H_M_theta, H_M_r] 
H_M_theta = cth_M.*H_z - sth_M.*cph_M.*H_x - sth_M.*sph_M.*H_y;
H_M_phi = cph_M.*H_y - sph_M.*H_x;

% Trap very small values of cos(theta_M) which can cause problems in the d(phi)/dt term of the LLG
% For multiple spins use array operations:
theta_near_z= abs(cth_M) < 0.0001; %logical variable, true if cos(theta) is small
% Aim is to stop cos(theta) being too small (1/cos(theta) in d(phi)/dt causes an infinity)
% sign must be preserved to keep the correct direction of rotation
dir=sign(cth_M); %dir is -1,0,1 for cos(theta): <0, 0, >0 
dir=dir+(abs(dir) < 0.5); %need to deal with dir=0. If dir=0 then add 1, otherwise add 0. 
% The line above converts dir=0 to dir=1 and leaves everything else unchanged.
cth_M = ~theta_near_z.*cth_M + theta_near_z.*dir.*0.0001; %Two complementary conditions added:
% if theta isn't near to the z-axis leave cos(theta) unchanged, 
% if theta is near to the z-axis set it to the lower limit value multiplied by the correct sign 

% Calculate dM/dt. Note that gamma and alpha are both negative.
dMdt_array(:,1) = (-LLG_alpha*H_M_phi - LLG_gamma*H_M_theta)./cth_M; 
dMdt_array(:,2) = LLG_gamma*H_M_phi - LLG_alpha*H_M_theta;

% ODE routine uses a row vector, convert the 2-D array back to the form
% assumed by the calling routine (may not really be necessary).
dMdt=reshape(dMdt_array,LLG_nelements*2,1);

