%% Finds the coercivity of a single layer ferromagnetic object
% The coercivity is the field value when the probability of switching is
% 0.5. Therefore when the probability is given by:
%
% $$ p_{switch} = 1 - \exp(-vt) $$
%
% where
%
% $$ v = f_0 \exp(-\frac{E_{barrier,1}(h,\theta)}{k_B T}) + f_0
% \exp(-\frac{E_{barrier,2}(h,\theta)}{k_B T}) $$
%
% $$f_0$ is the attempt frequency, $$k_B$ is the Boltzmann constant,
% $$E_{barrier,1}$ and $$E_{barrier,2}$ are the energy barriers, $$t$ is
% the time elapsed and $$T$ is the temperature. The coercivity can then be
% obtained from:
%
% $$ v(h_c,\theta) = -\frac{\ln 2}{t} $$
%
% where $$h_c$ is the coercivity. 
% 
% The switching field is calculated when the temperature is given as 0 K. 

%% Search for the coercivity
function [hc, tprob, xo, x] = hc_search(theta_h, h_keff, mu_o, m_s, vol, kb, temp, duration, attfreq)


xo = log(2)/(duration*attfreq);

    if temp > 0
        tempcheck = 1;
        disp('Temperature > 0 K')
    elseif temp == 0
        tempcheck = 2;
        disp('Temperature == 0 K')
    else
        msg = 'Negative temperature. Temperature must equal to or greater than zero';
        error(msg)
    end
    
    
    switch tempcheck
        case 1 
            disp('Calculating coercivity')
            
            h = -1e-6;
            tol = 1e-6;
            h_delta = 1e-2;
            x_delta = 10;

            while (abs(x_delta)>= tol)
                tprob = switchrateIntegrand(h, theta_h, h_keff, mu_o, m_s, vol, kb, temp);
                x = tprob -xo;
                x_delta = x/xo;
                if abs(tprob)<= xo
                    %disp('abs(tprob)<= xo')
                    h = h - h_delta;       
                    h_temp = h;
                elseif  abs(tprob)== 2
                    %disp('abs(tprob)== 2')
                    h = h_temp;
                else
                    %disp('abs(tprob)>= xo')
                    h_delta = h_delta/2;
                    h = h + h_delta;
                end
            end
            case 2
                disp('Calculating coercivity')
                h = sField(theta_h);
                tprob = 0;
                x = 0;
    end
% Coercivity in SI units
    hc = abs(h).*h_keff;
end

%% Calculate switching field
function h  = sField(theta)
t = nthroot(tan(theta), 3);
a = (1 - t.^2 + t.^4).^(1/2);
b = 1 + t.^2;
h = real(a./b); % real and positive, only dealing with magnitudes
end


%% Switching rate integral
function tprob = switchrateIntegrand(h, theta_h, h_keff, mu_o, m_s, vol, kb, temp)
theta_h = theta_h.*ones(size(h));
[eb1 eb2] = scaledbarrier(h, theta_h, h_keff, mu_o, m_s, vol, kb, temp);

eb1_check = isreal(eb1);
eb1 = eb1.*eb1_check;

eb2_check = isreal(eb2);
eb2 = eb2.*eb2_check;

tprob = exp(-eb1)+ exp(-eb2);

end

%% Scales the energy barrier
% The energy barrier that is calculated is in reduced units, this function
% converts the energy barrier into Joules. 
function [eb1 eb2] = scaledbarrier(h, theta_h, h_keff, mu_o, m_s, vol, kb, temp)
h_norm = -abs(h);
%h_norm = -abs(h)./h_keff

[ebarrier1, ebarrier2] = energybarrier(h_norm, theta_h);
eb1 = mu_o*m_s*vol*(1e-27)*h_keff.*ebarrier1/(kb*temp);
eb2 = mu_o*m_s*vol*(1e-27)*h_keff.*ebarrier2/(kb*temp);

end

%% Calculates the energy barrier
% This function is the same as energybarrier.m.  
% The energy barrier is calculated here for given values of the applied
% field, hv. Applied field values for a given head position can be
% calculated from scaledbarrier_sl.m that will generate hv from an imported
% field of head field values and associated field angles; theta_H_vec.
function [barrier1, barrier2] = energybarrier(hv, phiv)
% This functions returns the barrier ideally for h positive, to find the other barrier replace h by
% -h

%Set array sizes
phivec =ones(size(hv)).*phiv;
barrier1 = zeros(size(hv));
barrier2 = zeros(size(hv));


for j=1:length(hv)
    h=hv(j);
    phi=phivec(j);
    % Check if indeterminate. If hx=hy=hz=0, polar and azimuthal angles
    % become indeterminate => heff becomes indeterminate too
    if isnan(h) || isnan(phi)
        %disp('indeterminate')
        barrier1(j) = 0.5;
        barrier2(j) = 0.5;
    % Check if maxima and minima are degenerate 
    elseif sin(phi)==0 || sign(h)==0 
        %disp('maxima and minima are degenerate ')
        if abs(h) <= sField(phi) || sign(h)==0
            barrier1(j) = ((1 + h).^2)/2;
            barrier2(j) = ((1 - h).^2)/2;
        else
            barrier1(j) = 2.*abs(h);
            barrier2(j) = 2.*abs(h);
        end
    % Else calculate energy barrier   
    else
        %disp('calculating energy barrier')
        % Assemble energy barrier calculation by parts
        v = (h.^2.*sin(phi).*cos(phi) + sqrt((h.^2.*sin(phi).*cos(phi)).^2 + ((h.^2 -1)/3).^3 )).^2;
        w = (h.^2.*sin(phi).*cos(phi) - sqrt((h.^2.*sin(phi).*cos(phi)).^2 + ((h.^2 -1)/3).^3 )).^2;

        % use real parts 
        vcrt = nthroot(real(v), 3).*(real(v)==v) + v.^(1/3).*(real(v)~=v);
        wcrt = nthroot(real(w), 3).*(real(w)==w) + w.^(1/3).*(real(w)~=w);

        t = vcrt + wcrt;

        c1 = sqrt(t + h.^2.*sin(phi).^2 - 2*(h.^2 -1)/3);
        den1 = sqrt(t + h.^2.*sin(phi).^2 - 2*(h.^2 -1)/3);
        den2 = sqrt(t + h.^2.*cos(phi).^2 - 2*(h.^2 -1)/3);

        den1 = eps*(abs(den1)<=eps) + den1;
        den2 = eps*(abs(den2)<=eps) + den2;

        srt1 = sqrt(2*h.^2.*sin(phi).^2 - 4*(h.^2 -1)/3 - t + 2*h.*sin(phi).*(h.^2.*cos(phi).^2 + 1)./den1);
        srt2 = sqrt(2*h.^2.*cos(phi).^2 - 4*(h.^2 -1)/3 - t + 2*h.*cos(phi).*(h.^2.*sin(phi).^2 + 1)./den2);

        barrier1(j) = (c1 + h.*sin(phi)).*srt1/2 + h.*cos(phi).*srt2;

        h =-h;

        srt1 = sqrt(2*h.^2.*sin(phi).^2 - 4*(h.^2 -1)/3 - t + 2*h.*sin(phi).*(h.^2.*cos(phi).^2 + 1)./den1);
        srt2 = sqrt(2*h.^2.*cos(phi).^2 - 4*(h.^2 -1)/3 - t + 2*h.*cos(phi).*(h.^2.*sin(phi).^2 + 1)./den2);

        barrier2(j) = (c1 + h.*sin(phi)).*srt1/2 + h.*cos(phi).*srt2;
    end
end

end

%% Calculate switching field
function h  = sField(theta)
t = nthroot(tan(theta), 3);
a = (1 - t.^2 + t.^4).^(1/2);
b = 1 + t.^2;
h = real(a./b); % real and positive, only dealing with magnitudes
end

%% See also
%# [barrier1, barrier2] = energybarrier(hv, theta_H_vec)
