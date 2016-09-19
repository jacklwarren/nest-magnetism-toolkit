%% Energy barrier calculation for single layer
% Calculates the energy barrier of a single layer island using the method
% of calculation given by Kalezhi et al[1]. Energy of coupled composite
% islands in an applied magnetic field are calculated here assumes that the island is
% uniformly magnetised and that the change in magnetisation angle
% occurs at the interface.
% For a single layer the energy barrier can be expressed as:
%
% $$E = KV\sin^2(\theta) - mu_0 M V \cos(\theta - \theta_H) $$
%
% Where $$K$ is the total anisotropy of the material (including shape and
% crystalline).
% $$V$ is the layer volume
% $$M$ is the saturation magnetisation
% $$\theta$ is the polar angles of the magnetisation
% $$\theta_H$ is the applied field polar angle
% $H$ is the magnitude of the applied field
% The energy can be expressed in reduced units as
%
% $$ E_{reduced} = 1/2 \sin^2\theta - h \cos(\theta_H - \theta) $$
%
% Where E_reduced = Energy/(mu_0 M_s V H_K^{eff}), h = H/H_K^{eff}
% is the reduced applied field, \theta is the magnetisation angle
% and \theta_H is the field angle. The energy barrier can be
% obtained from the critical points, which are obtained by
% differentiating E_{reduced} with respect to \theta and setting to
% zero. 
%
% $$ \frac{dE_{reduced}}{d\theta} = \sin\theta\cos\theta - h\sin(\theta_H - \theta) = 0 $$
%
% The energy barrier can be numerically approximated by 
%
% $$ E_{barrier,approx} = \frac{mu_0 V M_s H_a}{2} (1 - \frac{h}{h_K(\theta_H)} )^{0.86 + 1.14 h_k(\theta_H)} $$
%
% Where V is the volume, M_s is the saturation magnetisation, H_a
% is the anisotropy field, h is the reduced applied field, \theta_H
% is the applied field angle and 
% $$ h_K(\theta_H) = [\cos^{2/3}\theta_H + \sin^{2/3}\theta_H]^{-3/2} $$
%
% which is the normalised switching field. 
% In reduced units the energy barrier is given here as:
%
% $$E_{barrier,approx} = \frac{1}{2} (1 -\frac{h}{h_K(\theta_H)})^{0.86+1.14h_k(\theta_H)} $$
%
% The energy barrier is thoroughly derived by Kalezhi et al. [p. 106-113,
% 2] and is given here as:
%
% $$ E_{barrier,1} = \frac{s_1 + h\sin\theta_H}{2} \sqrt{2h^2\sin^2\theta_H - \frac{4(h^2-1}{3} - t + \frac{2h\sin\theta_H(h^2\cos^2\theta_H + 1)}{s_1}}$$ 
%
% $$ + h\cos\theta_H\sqrt{2h^2\cos^2\theta_H - \frac{4(h^2-1}{3} - t + \frac{2h\cos\theta_H(h^2\sin^2\theta_H + 1)}{sqrt{t + h^2\cos^2\theta_H - (2/3)(h^2-1)}}}$$
%
% Where
%
% $$ s_1 = \sqrt{h^2\sin\theta_H + t - \frac{2(h^2-1)}{3}} $$
% 
% $$ t = v^{1/3} + w^{1/3} $$
%
% $$ v = \left[h^2\sin\theta_H\cos\theta_H + \sqrt{[h^2\sin\theta_H\cos\theta_H]^2 + [(h^2 - 1)/3]^2}\right]^2$$
%
% $$ w = \left[h^2\sin\theta_H\cos\theta_H - \sqrt{[h^2\sin\theta_H\cos\theta_H]^2 + [(h^2 - 1)/3]^2}\right]^2$$
%
%

%% Energy barrier calculation
% The energy barrier is calculated here for given values of the applied
% field, appliedfield, with associated applied field angles, theta_H. 
% For the output in SI units (Joules), anisotropy, saturation magetisation,
% and volume of the island must also be stated. Where anisotropy = k, 
% saturation magnetisation = ms and island volume = vol. 


function [energy_barrier1, energy_barrier2] = energybarrier_sl(appliedfield, theta_H, k, ms, vol)
% This functions returns the barrier ideally for h positive, to find the other barrier replace h by
% -h

% The permeability of free space is given by
 muo = 4*pi*1e-7; %in SI units

%The anisotropy field is calculated as: 
hk = 2*k/(muo*ms);

%The reduced field is given as:
h_reduced = appliedfield./hk;

%Set array sizes
thetavec =ones(size(h_reduced)).*theta_H;
barrier1 = zeros(size(h_reduced));
barrier2 = zeros(size(h_reduced));

for j=1:length(h_reduced)
   
    h=h_reduced(j);
    theta=thetavec(j);
    
    % Check if indeterminate. If hx=hy=hz=0, polar and azimuthal angles
    % become indeterminate => heff becomes indeterminate too
    if isnan(h) || isnan(theta)
        disp('indeterminate')
        barrier1(j) = 0.5;
        barrier2(j) = 0.5;
    % Check if maxima and minima are degenerate     
    elseif sin(theta)==0 || sign(h)==0 
        disp('maxima and minima are degenerate ')
        if abs(h) <= sField(theta) || sign(h)==0
            barrier1(j) = ((1 + h).^2)/2;
            barrier2(j) = ((1 - h).^2)/2;
        else
            barrier1(j) = 2.*abs(h);
            barrier2(j) = 2.*abs(h);
        end
    % Else calculate energy barrier    
    else
        disp('calculating energy barrier ')    
        % Assemble energy barrier calculation by parts
        v = (h.^2.*sin(theta).*cos(theta) + sqrt((h.^2.*sin(theta).*cos(theta)).^2 + ((h.^2 -1)/3).^3 )).^2;
        w = (h.^2.*sin(theta).*cos(theta) - sqrt((h.^2.*sin(theta).*cos(theta)).^2 + ((h.^2 -1)/3).^3 )).^2;

        % use real parts 
        vcrt = nthroot(real(v), 3).*(real(v)==v) + v.^(1/3).*(real(v)~=v);
        wcrt = nthroot(real(w), 3).*(real(w)==w) + w.^(1/3).*(real(w)~=w);

        t = vcrt + wcrt;

        c1 = sqrt(t + h.^2.*sin(theta).^2 - 2*(h.^2 -1)/3);
        den1 = sqrt(t + h.^2.*sin(theta).^2 - 2*(h.^2 -1)/3);
        den2 = sqrt(t + h.^2.*cos(theta).^2 - 2*(h.^2 -1)/3);

        
        den1 = eps*(abs(den1)<=eps) + den1;
        den2 = eps*(abs(den2)<=eps) + den2;

        srt1 = sqrt(2*h.^2.*sin(theta).^2 - 4*(h.^2 -1)/3 - t + 2*h.*sin(theta).*(h.^2.*cos(theta).^2 + 1)./den1);
        srt2 = sqrt(2*h.^2.*cos(theta).^2 - 4*(h.^2 -1)/3 - t + 2*h.*cos(theta).*(h.^2.*sin(theta).^2 + 1)./den2);

        barrier1(j) = (c1 + h.*sin(theta)).*srt1/2 + h.*cos(theta).*srt2;

        h =-h;

        srt1 = sqrt(2*h.^2.*sin(theta).^2 - 4*(h.^2 -1)/3 - t + 2*h.*sin(theta).*(h.^2.*cos(theta).^2 + 1)./den1);
        srt2 = sqrt(2*h.^2.*cos(theta).^2 - 4*(h.^2 -1)/3 - t + 2*h.*cos(theta).*(h.^2.*sin(theta).^2 + 1)./den2);

        barrier2(j) = (c1 + h.*sin(theta)).*srt1/2 + h.*cos(theta).*srt2;
    end
end

%Convert the energy barrier from reduced units to SI. 
energy_barrier1 = muo*ms*vol*(1e-27).*hk.*barrier1; % energy barrier in Joules
energy_barrier2 = muo*ms*vol*(1e-27).*hk.*barrier2; % energy barrier in Joules

end
%% Switching Field Calculation
function h  = sField(theta)
t = nthroot(tan(theta), 3);
a = (1 - t.^2 + t.^4).^(1/2);
b = 1 + t.^2;
h = real(a./b); % real and positive, only dealing with magnitudes
end

%% References
% [1] Josephat Kalezhi, Simon J. Greaves, Yasushi Kanai, Manfred E. Schabes,
% Michael Grobis, and Jim J. Miles. A statistical model of write-errors in
% bit patterned media. Journal of Applied Physics, 111(5):053926, 2012.
%
% [2] Josephat Kalezhi. Modelling Data Storage in NanoIsland Magnetic
% Materials. PhD thesis, University of Manchester, 2011.
