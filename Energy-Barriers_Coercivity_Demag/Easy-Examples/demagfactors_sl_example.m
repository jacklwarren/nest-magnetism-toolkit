%% Basic Example of Calculating Demagnetising Factors for single layer
% This is an example of calculating demagnetising factors using the
% demagfactors function for truncated elliptic cones, elliptic cylinders
% and prisms. An elliptic cylinder has been specified below as a single
% layer configuration. All units are assumed to be in SI units, and
% dimensions are specified in nm.


% First, clear all previous values
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


% For single layer case only one thickness must be specified, for ecc two
% thicknesses must be specified, see demagfactors_ecc.m example file.
t_tot = 6; % total thickness in nm

%% Specify single layer magnetic properies

demag_tol = 1e-6; % tolerance used in some integrals of demagnetising factors
ms = 700e3; % saturation magnetisation in A/m:  700emu/cc 
k1 =7.80e5; % in J/m^3 muo*ms*hk/2;  %  crystalline anisotropy constant in J/m^3: since hk = 2*k1./(muo*ms) => k1 = muo*ms*hk/2

% calculate demagnetising factors
[nxx nyy nzz] = demagfactors(a, b, t_tot, alpha, islandgeo, demag_tol) 
 


%%

