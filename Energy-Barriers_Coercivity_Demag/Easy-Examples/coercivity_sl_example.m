%% Basic example of calculating switching field and coercivity for single layer
% This script demonstrates the basic function of hc_search.m function to function coercivity in a single layer case. % % Switching field is given by hc_search.m when temperature = 0. The method given here accounts % of a demagnetising % % field of the object and requires the demagfactor.m function; this can be left out and the method is shown below.

%% Choose Island Geometry
% The demagfactors.m function supports three different kinds of geometry.
% Here elliptic cylinder has been used, alternatively a truncated elliptic
% cone or prism could be used.
% ELLIPTIC CYLINDER:
% Specify island parameters
% Comment/uncomment the following lines to activate/deactivate this
% geometry

islandgeo = 2;
a = 8/2; % semi major axis in nm 
b = a; % semi minor axis in nm
alpha = 1; % any number will do, used for truncated elliptic cone only, needed to complete island_prop
% For single layer case only one thickness must be specified, for ecc two
% thicknesses must be specified, see demagfactors_ecc.m example file.
t_tot = 6; % total thickness in nm
vol= pi*a.*b.*t_tot;
area = pi*a*b;

% Example for truncated elliptic cone:

% islandgeo = 1; % 1 represents truncated elliptic cone
% a = 12.5/2; % semi major axis in nm
% b = a; % semi minor axis in nm
% alpha = 0.66; % ratio of top to bottom semi major axis

% Example for prism:
% islandgeo = 3;
% a = (18/sqrt(2))/2;   % half length in nm: 50% fill factor
% b = a;  % half width in nm
% alpha = 1; % any number will do, used for truncated elliptic cone only, needed to complete island_prop

%% Specify magnetic properties

demag_tol = 1e-6; % tolerance used in some integrals of demagnetising factors
ms = 700e3; % saturation magnetisation in A/m:  700emu/cc 
k1 =7.80e5; % in J/m^3 muo*ms*hk/2;  %  crystalline anisotropy constant in J/m^3: since hk = 2*k1./(muo*ms) => k1 = muo*ms*hk/2

%calculate demagnetising factors  
[nxx nyy nzz] = demagfactors(a, b, t_tot, alpha, islandgeo, demag_tol); 

%% Specify uniform applied field magnitude and angles 
% Specify an applied field, this can be expressed as an array of values.
% Calculating for energy barrier will be made using an applied field
% reduced by the effective field, a reduced head field can be used from the
% offset if this is accounted for, i.e. h_reduced is commented out. The
% field must be negatively applied in an array, e.g. -1 to 1.
% 
% The energy barrier is calculated in reduced units, therefore the applied
% field must also be in reduced units. Since the effective anisotropy field is
% calculated from:
%
% $$ H_{K,eff} = \frac{2 K_{1,eff}}{\mu_o M_s}
%              = H_k + M_s(N_{xx} \cos(\phi_H)^2 + N_{yy} \sin(\phi_H)^2 -
%              N_{zz}) $$
%
% An effective anisotropy can be calculated from our given value, which is
% given by:
%
% $$ K_{1,eff}= K_1 + \frac{1}{2}\mu_0 M_s^2(N_{xx} \cos(\phi_H)^2 + N_{yy}
% \sin(\phi_H)^2 - N_{zz}) $$
%


%Specify angle of applied field, can be expressed as a single value or as
%an array the same size as the applied field,h. 
theta_H = -pi; % field polar angle in radians (the usual field angle). can be an array of values, same size as h

phi_H = pi/2; % field azimuthal angle in radians (can be an array of values, same size as h)

% The permeability of free space is given by
 muo = 4*pi*1e-7; %in SI units

%Accounting for demagnetising factors, i.e. including shape anisotropy, the
%anisotropy is given as:
k = k1 %+ (1/2)*muo*ms.^2.*(nxx.*cos(phi_H).^2 + nyy.*sin(phi_H).^2 - nzz);
 
%The anisotropy field is calculated as: 
hk = 2*k/(muo*ms);

%% Calculate coercivity and switching field
kB = 1.3806503*1e-23; % J/T; % Boltzmann constant in SI units
temp = 0:10:300; % absolute temperature in Kelvin

duration = 1*1e-9; %duration in seconds
attfreq = 1e11; % attempt frequency in Hz;




% Calculating switching field. The coercivity will equal the coercivity at
% zero temperature. 


% Calculating coercivity 
for i = 1:length(temp)

[hc(i), tprob, xo, x]= hc_search(theta_H, hk, muo, ms, vol, kB, temp(i), duration, attfreq);

end

%% Plot coercivity against temperature
figure(1)

fontsiz= 10; 
fontsiz_text= 10; 
fontsiz_leg=10; 
linewid= 5; 
markersizb=10;



clf
grid on
hold on

xlabel('Temperature(K)','FontSize',fontsiz);
ylabel('Coercivity (A/m)','FontSize',fontsiz);
plot(temp, hc,'-o','LineWidth',linewid,'Color','blue','MarkerSize',markersizb,'MarkerFaceColor','blue');
