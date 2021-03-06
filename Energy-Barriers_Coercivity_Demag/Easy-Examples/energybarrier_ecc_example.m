%% Example of calculating energy barrier for ecc 
% This script demonstrates the basic function  energy_barrier_hess.m
% function. The function finds the energy barrier of a two-layer shape. The
% method given here accounts of a demagnetising field of the object and
% requires the demagfactor.m function; this can be left out and the method
% is shown below.

clear all;
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
alpha = Inf; % any number will do, used for truncated elliptic cone only, needed to complete island_prop
% For single layer case only one thickness must be specified, for ecc two
% thicknesses must be specified, see demagfactors_ecc.m example file.
t_h = 4;%6;%4; % hard layer thickness in nm , t_h + t_s = t % if working with single layer, this should be equal to t_tot!
t_s = 2;%0;%; % soft layer thickness in nm % if working with singlelayer, this should be zero!
t_tot = t_h + t_s; % total thickness in nm, used for island thickness for single layer, or can use t_h
% A spacing between layers can be specified in nm, here set to 0 nm
spacing = 0; 
vol= pi*a.*b.*t_tot;
area = pi*a*b;
v_h = pi*a.*b.*t_h; % volume in nm^3
v_s = pi*a.*b.*t_s;

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
m_h = 1000*1e3;  % hard layer saturation magnetisation in A/m: % 400 emu/cm^3
m_s = 1400*1e3; % soft layer saturation magnetisation in A/m: % 1000 emu/cm^3


k_h_c = 1.16E+06;   %  hard layer crystalline anisotropy in J/m^3   % = 15*1e6 erg/cm^3  %1 erg/(cm^3)=0.1joules/(m^3)
k_s_c = 0.1*1*1e6;    % soft layer crystalline anisotropy in J/m^3   % = 1*1e6 erg/cm^3

j_xc = 5*1e-3;        % exchange coupling in J/m^2   % = 5erg/cm^2   % 1erg/(cm^2)=0.001joules/(m^2)

demag_tol = 1e-6; % tolerance used in some integrals of demagnetising factors


%Calculate demagnetising factors. Calculation of energy barrier can be
%performed without demagnetising factors. However, here they have been
%included to demonstrate how they can be incorporated. To remove them do
%not use an effective anisotropy, or start with an applied field that is
%already reduced.
[nxx_h  nyy_h nzz_h] = demagfactors(a, b, t_h, alpha, islandgeo, demag_tol);
[nxx_s  nyy_s nzz_s] = demagfactors(a, b, t_s, alpha, islandgeo, demag_tol);

%% Specify uniform applied field magnitude and angles 
% Comment out if using single layer Specify an applied field, this can be
% expressed as an array of values. Calculating for energy barrier will be
% made using an applied field reduced by the effective field, a reduced
% head field can be used from the offset if this is accounted for, i.e.
% h_h_reduced and h_s_reduced are commented out. A field must be specified
% for both hard and soft layers, but it can be the same value. The field
% must be negatively applied in an array, e.g. -1 to 1.
% The function energybarrier.m uses reduced units. Therefore the applied
% field must by calculated in reduced units. Since the effective field is
% calculated from:
%
% $$ H_{K,eff} = \frac{2 K_{1,eff}}{\mu_o*M_s}
%              = H_k + M_s(N_{xx} \cos(\phi_H)^2 + N_{yy} \sin(\phi_H)^2 -
%              N_{zz}) $$
%
% An effective anisotropy can be calculated from our given values, which is
% given by:
%
% $$ K_{1,eff}= K_1 + 0.5\mu_0*M_s^2(N_{xx} \cos(\phi_H)^2 + N_{yy}
% \sin(\phi_H)^2 - N_{zz}) $$
%
% The effective anisotropy can be used to calculate the effective field and
% the reduced field is then:
% 
% $$ H_{reduced} = \frac{H}{H_{K,eff}} $$
%
% For an ECC, or bilayer, system this must be done for both hard and soft
% layers. 


h_h = (-8e5:0.8e5:8e5); % negative applied field in A/m (can be an array of values)
h_s = (-5e5:0.5e5:5e5); % negative applied field in A/m (can be an array of values)

theta_H_h = -pi:pi/10:pi;% field polar angle in radians (the usual field angle). can be an array of values, same size as h
theta_H_s = -pi/2:pi/20:pi/2; % field polar angle in radians (the usual field angle). Can be an array of values, same size as h

phi_H_h = pi/2; % field azimuthal angle in radians (can be an array of values, same size as h)
phi_H_s = pi/2; 

% The permeability of free space
muo = 4*pi*1e-7; %in SI units

%For both hard and soft layer the anisotropy, accounting for the
%demagnetising field is then given by:
k_h = k_h_c  + (1/2)*muo*m_h.^2.*(nxx_h.*cos(phi_H_h).^2 + nyy_h.*sin(phi_H_h).^2 - nzz_h)
k_s = k_s_c  + (1/2)*muo*m_s.^2.*(nxx_s.*cos(phi_H_s).^2 + nyy_s.*sin(phi_H_s).^2 - nzz_s)


%% Calculate energy barriers
% Note there are two energy barriers, barrier1_reduced is the usual value.
% A complex value of energy barrier indicates the magnetisation has
% reversed, i.e. no energy barrier

% The energy barrier for ECC must specify more details than in the single
% layer case (where only the applied field and angle must be known). This
% is due to the effect of the two layers interacting with each other.

[energy_barrier1, energy_barrier2]= energybarrier_ecc(h_h, h_s, theta_H_h, theta_H_s, j_xc, k_h_c, k_s_c, m_h, m_s,v_h,v_s, area,t_h, t_s, a,spacing)


% Calculate energy barrier in k_{B}T.  
kB = 1.3806503*1e-23; % Boltzmann constant in J/T;
temp = 300; % temperature in Kelvin

energy_barrier1_kBT = energy_barrier1/(kB*temp); % energy barrier in units of kBT
energy_barrier2_kBT = energy_barrier2/(kB*temp); % energy barrier in units of kBT

%% Plot graphs
fontsiz= 20; % for intermag
linewid= 5; 
markersizb=10;



figure(1)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz)

xlabel('Applied Field (A/m)','FontSize',fontsiz);
ylabel('Energy Barrier (J)','FontSize',fontsiz);

plot(h_h, energy_barrier1,'-o','LineWidth',linewid,'Color','blue','MarkerSize',markersizb,'MarkerFaceColor','blue');
plot(h_s, energy_barrier1,'-o','LineWidth',linewid,'Color','red','MarkerSize',markersizb,'MarkerFaceColor','red');

plot(h_h, energy_barrier2,'-o','LineWidth',linewid,'Color','green','MarkerSize',markersizb,'MarkerFaceColor','green');
plot(h_s, energy_barrier2,'-o','LineWidth',linewid,'Color','magenta','MarkerSize',markersizb,'MarkerFaceColor','magenta');
plot_leg = legend('Hard Layer Energy barrier 1','Soft Layer Energy barrier 1', 'Hard Layer Energy barrier 2','Soft Layer Energy barrier 2');   
set(plot_leg,'FontSize',fontsiz);

figure(2)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz)

xlabel('Field Angle (rads)','FontSize',fontsiz);
ylabel('Energy Barrier (k_B T)','FontSize',fontsiz);

plot(h_h, energy_barrier1_kBT,'-o','LineWidth',linewid,'Color','blue','MarkerSize',markersizb,'MarkerFaceColor','blue');
plot(h_s, energy_barrier1_kBT,'-o','LineWidth',linewid,'Color','red','MarkerSize',markersizb,'MarkerFaceColor','red');

plot(h_h, energy_barrier2_kBT,'-o','LineWidth',linewid,'Color','green','MarkerSize',markersizb,'MarkerFaceColor','green');
plot(h_s, energy_barrier2_kBT,'-o','LineWidth',linewid,'Color','magenta','MarkerSize',markersizb,'MarkerFaceColor','magenta');
plot_leg = legend('Hard Layer Energy barrier 1','Soft Layer Energy barrier 1', 'Hard Layer Energy barrier 2','Soft Layer Energy barrier 2');   
set(plot_leg,'FontSize',fontsiz);
