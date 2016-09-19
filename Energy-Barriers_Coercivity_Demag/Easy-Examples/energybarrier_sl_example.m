%% Example of calculating energy barrier for singe layer
% This script demonstrates the basic functions of the energybarrier_sl.m
% function. The function finds the energy barrier single layer an arbitrary
% ferromagnetic shape subject to an applied field. However, the method given here
% accounts of a demagnetising field of the object and requires the
% demagfactor.m function; this can be left out and the method is shown
% below.

% Clear all previous values
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
k =7.80e5; % in J/m^3 muo*ms*hk/2;  %  crystalline anisotropy constant in J/m^3: since hk = 2*k1./(muo*ms) => k1 = muo*ms*hk/2

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
% The effective anisotropy can be used to calculate the effective field and
% the reduced field is then:
% 
% $$ H_{reduced} = \frac{H}{H_{K,eff}} $$
%

% An array of applied fields
h = (-5e5:0.5e5:5e5);% in A/m

%Specify angle of applied field, can be expressed as a single value or as
%an array the same size as the applied field,h. 
theta_H = -pi:(pi/10):pi; % field polar angle in radians (the usual field angle). can be an array of values, same size as h

phi_H = pi/2; % field azimuthal angle in radians (can be an array of values, same size as h)



%Accounting for demagnetising factors, i.e. including shape anisotropy, the
%anisotropy is given as:
%k = k1 + (1/2)*muo*ms.^2.*(nxx.*cos(phi_H).^2 + nyy.*sin(phi_H).^2 - nzz);
 

%% Calculate energy barriers
% Note there are two energy barriers and barrier1_reduced is the usual
% value. A complex value of energy barrier indicates the magnetisation has
% reversed, i.e. no energy barrier


[energy_barrier1, energy_barrier2] = energybarrier_sl(h, theta_H, k, ms, vol); % reduced energy barrier


% Calculate energy barriers in k_{B}T. 
kB = 1.3806503*1e-23; % Boltzmann constant in J/T;
temp = 300; % temperature in Kelvin
% 
%In units of k_{B}T, the energy barrier is given as:
energy_barrier1_kBT = energy_barrier1/(kB*temp) % energy barrier in units of kBT
energy_barrier2_kBT = energy_barrier2/(kB*temp) % energy barrier in units of kBT

%% Plot energy barrier against field 



fontsiz= 30; 
fontsiz_text= 10; 
fontsiz_leg=10; 
linewid= 5; 
markersizb=10;


figure(1)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz)

xlabel('Applied Field (A/m)','FontSize',fontsiz);
ylabel('Energy Barrier (k_B T)','FontSize',fontsiz);
plot(h, energy_barrier1,'-o','LineWidth',linewid,'Color','red','MarkerSize',markersizb,'MarkerFaceColor','red');
plot(h, energy_barrier2,'-o','LineWidth',linewid,'Color','blue','MarkerSize',markersizb,'MarkerFaceColor','blue');
plot_leg = legend('Energy barrier 1','Energy barrier 2');   
set(plot_leg,'FontSize',fontsiz);

figure(2)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz)

xlabel('Field Angle (rads)','FontSize',fontsiz);
ylabel('Energy Barrier (k_B T)','FontSize',fontsiz);
plot(theta_H, energy_barrier1_kBT,'-o','LineWidth',linewid,'Color','red','MarkerSize',markersizb,'MarkerFaceColor','red');
plot(theta_H, energy_barrier2_kBT,'-o','LineWidth',linewid,'Color','blue','MarkerSize',markersizb,'MarkerFaceColor','blue');
plot_leg = legend('Energy barrier 1','Energy barrier 2');   
set(plot_leg,'FontSize',fontsiz);
%% References
% [1] Kalezhi, Josephat, et al. "A statistical model of write-errors in bit
% patterned media." Journal of Applied Physics 111.5 (2012): 053926.
