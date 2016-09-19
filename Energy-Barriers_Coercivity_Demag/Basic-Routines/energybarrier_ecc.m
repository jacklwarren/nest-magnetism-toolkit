function [energy_barrier1, energy_barrier2]= energybarrier_ecc(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, j_xc, k_h, k_s, m_h, m_s, v_h, v_s, area, t_h, t_s, a, spacing)

%% Energy barrier calculation for exchange coupled composites  
% Calculates the energy barrier of a double layer island using the method
% of calculation given in (Kalezhi et al,2012). Energy of coupled composite
% islands in an applied magnetic field are calculated here using a two-spin
% approximation that assumes that the hard and soft layers of the islands
% are uniformly magnetised and that the change in magnetisation angle
% occurs at the interface.
% For a single layer the energy barrier can be expressed as:
%
% $$E = K V \sin^2(\theta_h) - \mu_0 M V\cos(\theta - \theta_H)$$
%
% For two layers the total energy assuming in-plane rotation of the field
% is expressed as:
%
% $$E = K_h V_h\sin^2(\theta_h) + K_s V_s\sin^2(\theta_s) -
% JA\cos(\theta_s - \theta_h) - \mu_0 M_h V_h\cos(\theta_h - \theta_H)$$
%
% $$-\mu_0 M_s V_s\cos(\theta_s - \theta_H)-
% \mu_0\frac{M_s V_s M_h V_h}{4\pi d^3}[2 S_2\cos(\theta_s)cos(\theta_h) -
% S_1\cos(\theta_s)cos(\theta_h)]$$
%
% Where $$K$ is the total anisotropy of the material (including shape and
% crystalline). $$V$ is the layer volume, $$t$ is the thickness of the
% material $$A$ is the cross-sectional area, $$d$ is the spacing between
% the centres of soft and hard layers, $$M$ is the saturation
% magnetisation. $$\theta_h$ and $$\theta_s$ are the polar angles of the
% magnetisation in the hard and soft layers respectively. $$\theta_H$ is
% the applied field polar angle, $$H$ is the magnitude of the applied
% field, $$J$ is the interlayer exchange coupling constant.The subscripts
% $$_s$ and $$_h$ denote the soft and hard layers. $$S_1$ and $$S_2$
% represent factors that depend on the relative position and shape of the
% layers and are calculated elsewhere.
%
% The energy barrier is calculated here for given values of the applied
% field in the hard and soft layers, hp_h_vec and hp_s_vec. Applied field
% values for a given head position can be calculated from scaledbarrier.m
% that will generate hp_h_vec and hp_s_vec from an imported field of
% head field values and associated field angles; theta_H_s_vec and
% theta_H_s_vec. Alternatively these values can be input directly as an
% array. 
%
% To simplify the equation, the terms used for magnetostatic interaction
% between layers.
%
% $$ H_{dip} = M_s V_s /(4\pi d^3 H_{k,h})$$
%
% Where $$M_s$ is the saturation magnetisation for soft layer, similarly
% $$M_h$ is the saturation magnetisation for the hard layer. $$V_s$ is the
% volume of soft layer and $$d$ is the spacing between the centres of the
% layers
% 
% $$H_{k,h} = 2K_h/(\mu_0 M_h) + M_h(N_{xx,h} \cos(\phi_H)^2 +
% N_{yy,h}\sin(\phi_H)^2 - N_{zz,h})$$
%
% Where $$N_{xx,h}$ ,  $$N_{yy,h}$ and $$N_{zz,h}$ are demagnetisation factors that can be
% calculated from the function demagfactors.m from the island properties
%
% $$ k = (K_s V_s)/(K_h V_h)$$
%
% $$ m = (M_s V_s)/(M_h V_h) $$
%
%
% The energy minima are calculated as:
%
% $$  E' = \frac{1}{2}\sin^2\theta_h + \frac{K_s V_s}{K_h V_h}\sin^2\theta_s - H_h\cos(\theta_h - \theta_H) - \frac{M_sV_s}{M_h V_h}H_s\cos(\theta_s - \theta_H) - J_{exch}\cos(\theta_s - \theta_h) - \frac{M_s V_s}{4\pi d^3 H_{k,h}} (2 S_2\cos(\theta_s)\cos(\theta_h) - S_1\sin(\theta_s)sin(\theta_h))$$
% 
%
% The saddle point is calculated as:
% 
% $$  E'(x) = \sin\theta_h\cos\theta_h +  H_h\sin(\theta_h - \theta_H) - J_{exch}\sin(\theta_s- \theta_h) + \frac{M_s V_s}{4\pi d ^3 H_{k,h}}(2 S_2\cos\theta_s\sin\theta_h+ S_1\sin\theta_s\cos\theta_h)$$
%
% $$ E'(y) =  \frac{K_s V_s}{K_h V_h}\sin\theta_s\cos\theta_s + \frac{M_s V_s}{M_h V_h}H_s\sin(\theta_s - \theta_H) + J_{exch}\sin(\theta_s-\theta_h) + H_{dip}(2 S_2\sin\theta_s\cos\theta_h+ S_1\cos\theta_s\sin\theta_h)$$ 
%
% 
%  
% 
% 

%% Initialise
% The permeability of free space
muo = 4*pi*1e-7; 

%Set array sizes
ebarrier1 = zeros(size(theta_H_h_vec));
ebarrier2 = zeros(size(theta_H_h_vec));

tol = 1e-6; % the tolerance used between two angles (e.g. if theta_1 - theta_2 <= tol then they are the same)
stop_count = 10; %if the energy barrier is stuck and the minimum cannot be found within the stop_count then it exits

r = 1e-20; %eps, rho value for Fourier calculation for magnetostatic interaction

%Set options for finding 
options = optimset('Display','off','TolFun',1e-15, 'TolX',1e-15);

% The effective field in the hard layer is given as:
h_k_h = 2*k_h/(muo*m_h); % in A/m
h_k_s = 2*k_s/(muo*m_s); % in A/m

% The reduced field is given by:
h_h_reduced = hp_h_vec./h_k_h;
h_s_reduced = hp_s_vec./h_k_s;

% Reduced interlayer exchange 
j_exch = area.*j_xc./(2*h_k_h.*v_h).*1e9;

% Overall anisotropy can be calculated as
% If given there is only on value of anisotropy and no demagnetising
% factors used then use this value only
k =k_s.*v_s./(k_h.*v_h);

% Saturation magnetisation overall
m = m_s.*v_s/(m_h.*v_h);

% calculate distance between soft and hard layer centres in nm, spacing has
% been set to zero here so it can be ignored
d = ((t_h + t_s)/2)+spacing;

% Effect of interaction between layers (constants)
h_dip = m_s.*v_s./(4*pi*d.^3.*h_k_h);


%% Magnetostatic interaction
% Magnetostatic interaction is included in the equation of energy as the
% term
%
% $$ -\mu_0 \frac{M_{s,s}V_{s}M_{s,h}V_{h}}{4\pi d^3}*
% [2S_2\cos\theta_s\cos\theta_h - S_1\sin\theta_s\sin\theta_h] $$
%
% Where M_{s,s} is the saturation magnetisation for the soft layer, M_{s,h}
% is the saturation magnetisation for the hard layer. V is the volume of
% each layer, where _{h} or _{s} denote hard and soft layers. mu_0 is
% permeability of free space, d is the spacing between the two layers and
% theta_h and theta_s are the field angle in the hard and soft layer. The
% two factors, S_1 and S_2, are calculated from a Fourier space integral
% that depend on the relative shape and position of the layers. More
% details are given in [1]. The factors are calculated below.
%
[s1_factor, s2_factor, ~, ~] = S_factors_fast(d, r, a, a, t_h, t_s);


%% Calculate energy barrier loop
    if outofboundscheck(h_h_reduced, h_s_reduced, theta_H_h_vec, theta_H_s_vec) == 1   
        for i=1:length(theta_H_h_vec)
            %% First Minimum
            % Begin searching for the first minimum of the function between
            % 0 and 0.
            theta_H_h = theta_H_h_vec(i);
            theta_H_s = theta_H_s_vec(i);
            hp_h = h_h_reduced(i);
            hp_s = h_s_reduced(i);


            [th_m1,~,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),[-1e-6,-1e-6],options);

            n_count = 0;
            while (exitflag~=1 && n_count<=stop_count)
                disp('energy barrier calculation: stuck at 1st minimum')
                [th_m1,~,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_m1,options);
                n_count = n_count +1;
            end
            %% Second Minimum
            % Begin searching for the first minimum between of the function
            % between -pi and -pi.
            [th_m2,~,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),[-pi,-pi],options);

            n_count = 0;
            while (exitflag~=1 && n_count<=stop_count)
                disp('energy barrier calculation: stuck at 2nd minimum')
                [th_m2,~,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_m2,options);
                n_count = n_count +1;
            end

            %% Direct root for finding saddlepoint
            % The saddle point of the function is found from the psoitions
            % of the first and second minimum

            % th_m1 is the first minimum th_m2 is the second minimum

            % Checks that the two angles are the same within the given
            % tolerance if they are then the saddle point is searched for
            % at half the value of the second minimum.
            if (abs(th_m2(1)-th_m1(1))<= tol && abs(th_m2(2)-th_m1(2))<= tol)
                th_s0 = th_m2/2;
            else
                th_s0 = (th_m1 + th_m2)/2;
            end
            % Searching for the saddle point using the starting location
            % found from above
            opts = optimset('Display','off','TolFun',1e-15, 'TolX',1e-15);
            [th_s,~,exitflag] = fsolve(@(x)direct_root(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_s0,opts);  % Call optimizer

            n_count = 0;
            set = [1 2 3 -3];

            while ((~ismember(exitflag, set)) && n_count<=stop_count)
                disp('energy barrier calculation: stuck at saddle point')
                th_s0 = th_s;
                [th_s,~,exitflag] = fsolve(@(x)direct_root(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_s0,opts);  % Call optimizer
                n_count = n_count +1;
            end

            %% Find energy barrier
            x(1) = th_s(1); % first root of the magnetisation angle theta
            x(2) = th_s(2); % second root of the magnetisation angle theta
            x(3) = hp_h; % applied field in the hard layer
            x(4) = hp_s; % applied field in the soft layer

            % If the determinant of the Hessian, delta_s, is greater than
            % zero then the first energy barrier is zero and the second is
            % calculated from the energy barrer calculation. Otherwise both
            % energy barriers are calculated from the energy barrier
            % calculation. Finally if the energy barrier is calculated as a
            % negative value it is assumed to be equal to zero.
            delta_s = hessian_det(x, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor);

            if delta_s >= 0 || (abs(th_m2(1)-th_m1(1))<= tol && abs(th_m2(2)-th_m1(2))<= tol)
                ebarrier1(i) = 0;
                ebarrier2(i) = energy(hp_h, hp_s, theta_H_h, theta_H_s, th_s(1), th_s(2), j_exch, h_dip, k, m, s1_factor, s2_factor) - energy(hp_h, hp_s, theta_H_h, theta_H_s, th_m2(1), th_m2(2), j_exch, h_dip, k, m, s1_factor, s2_factor);
            else
                ebarrier1(i) = energy(hp_h, hp_s, theta_H_h, theta_H_s, th_s(1), th_s(2), j_exch, h_dip, k, m, s1_factor, s2_factor) - energy(hp_h, hp_s, theta_H_h, theta_H_s, th_m1(1), th_m1(2), j_exch, h_dip, k, m, s1_factor, s2_factor);
                ebarrier2(i) = energy(hp_h, hp_s, theta_H_h, theta_H_s, th_s(1), th_s(2), j_exch, h_dip, k, m, s1_factor, s2_factor) - energy(hp_h, hp_s, theta_H_h, theta_H_s, th_m2(1), th_m2(2), j_exch, h_dip, k, m, s1_factor, s2_factor);
            end

            if ebarrier1(i)< 0
                ebarrier1(i)=0;
            end
        end
    else 
        msg = ('Input field array and field angle array are of different dimensions. Please use the same array size for both values.');
        error(msg);
        return;
    end
    
    %Convert the energy barrier from reduced units to SI.     
    energy_barrier1 = 2*k_h*v_h*(1e-27).*ebarrier1; % energy barrier in Joules
    energy_barrier2 = 2*k_h*v_h*(1e-27).*ebarrier2; % energy barrier in Joules

end 

%% Calculating first and second minimum
function out = energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor)
    theta_h = x(1);
    theta_s = x(2);
    out = (1/2)*(sin(theta_h).^2 + k.*sin(theta_s).^2) - hp_h.*cos(theta_h - theta_H_h) - m.*hp_s.*cos(theta_s - theta_H_s) - j_exch.*cos(theta_s -theta_h) - h_dip.*(2*s2_factor*cos(theta_s).*cos(theta_h) - s1_factor.*sin(theta_s).*sin(theta_h));
end

%% To find saddlepoint
% Use fsolve of the function to find the saddle point of the energy barrier
% theta_h = x(1), field angle for the hard layer
% theta_s = x(2), field angle for the soft layer
function out = direct_root(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor)
    out = [sin(x(1))*cos(x(1)) + hp_h*sin(x(1)-theta_H_h) - j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*sin(x(1)) + s1_factor*sin(x(2)).*cos(x(1)));
          k*sin(x(2))*cos(x(2)) + m*hp_s*sin(x(2)-theta_H_s) + j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*sin(x(2))*cos(x(1))+ s1_factor*cos(x(2)).*sin(x(1)))];
end


%% Calculate energy
function out = energy(hp_h, hp_s, theta_H_h, theta_H_s, theta_h, theta_s, j_exch, h_dip, k, m, s1_factor, s2_factor)
    out = (1/2)*(sin(theta_h).^2 + k.*sin(theta_s).^2) - hp_h.*cos(theta_h - theta_H_h) - m.*hp_s.*cos(theta_s - theta_H_s) - j_exch.*cos(theta_s -theta_h)- h_dip.*(2*s2_factor*cos(theta_s).*cos(theta_h) - s1_factor*sin(theta_s).*sin(theta_h));
end

%% Hessian determinant
% Returns the determinant of the energy barrier hessian matrix
function [delta, e_x, e_y, e_xx, e_yy, e_xy] = hessian_det(x, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor)

    e_x = sin(x(1))*cos(x(1)) + x(3)*sin(x(1)-theta_H_h) - j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*sin(x(1)) + s1_factor*sin(x(2)).*cos(x(1)));
    e_y = k*sin(x(2))*cos(x(2)) + m*x(4)*sin(x(2)-theta_H_s)+ j_exch*sin(x(2)-x(1))+ h_dip*(2*s2_factor*sin(x(2))*cos(x(1)) + s1_factor*cos(x(2)).*sin(x(1)));

    e_xx = cos(x(1))^2 - sin(x(1))^2 + x(3)*cos(x(1) - theta_H_h) + j_exch*cos(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*cos(x(1)) - s1_factor*sin(x(2)).*sin(x(1)));
    e_yy = k*(cos(x(2))^2 - sin(x(2))^2) + m*x(4)*cos(x(2)- theta_H_s) + j_exch*cos(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*cos(x(1))- s1_factor*sin(x(2)).*sin(x(1)));
    e_xy  = -j_exch*cos(x(2) - x(1)) + h_dip*(-2*s2_factor*sin(x(2))*sin(x(1)) + s1_factor*cos(x(2))*cos(x(1)));
    delta = e_xx*e_yy - e_xy*e_xy;
end

%% Check array out of bounds, 
% Could change this so that both arrays are made equal? 
function  check = outofboundscheck(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec)
    if (((length(hp_h_vec))==(length(hp_s_vec)))&&((length(hp_h_vec))==(length(theta_H_h_vec)))&&((length(theta_H_s_vec))==(length(theta_H_h_vec))))
        check = 1;
    % Index out of bounds
    else
        check = 2;
    end    
end


%% Calculate magnetostatic interaction factors between layers using a Fourier space integral
function [S_1_factor, S_2_factor, S_3_factor, S_4_factor] = S_factors_fast(Z, R_rho, a_1, a_2, t_1,t_2)

    k_rho_min = 0;  
    k_rho_max = Inf; 

    rel_tol = 1e-8;
    abs_tol = 1e-12;
    max_cnt = 1000;

 
    S_1_intgral = quadgk(@(x)S_1_integrand(x, a_1, a_2, t_1, t_2, R_rho, Z), k_rho_min, k_rho_max,'RelTol', rel_tol,'AbsTol', abs_tol,'MaxIntervalCount', max_cnt);

    S_2_intgral = quadgk(@(x)S_2_integrand(x, a_1, a_2, t_1, t_2, R_rho, Z), k_rho_min, k_rho_max,'RelTol', rel_tol,'AbsTol', abs_tol, 'MaxIntervalCount', max_cnt);

    S_3_intgral = quadgk(@(x)S_3_integrand(x, a_1, a_2, t_1, t_2, R_rho, Z), k_rho_min, k_rho_max,'RelTol', rel_tol,'AbsTol', abs_tol, 'MaxIntervalCount', max_cnt);

    S_4_intgral = quadgk(@(x)S_4_integrand(x, a_1, a_2, t_1, t_2, R_rho, Z), k_rho_min, k_rho_max,'RelTol', rel_tol,'AbsTol', abs_tol, 'MaxIntervalCount', max_cnt);

    V_1 = pi*a_1^2*t_1;
    V_2 = pi*a_2^2*t_2;

    R = sqrt(Z.^2 + R_rho.^2);
    S_1_factor = 8*pi*R^3*4*a_1*a_2*S_1_intgral/(R_rho*V_1*V_2);
    S_2_factor = 8*pi*R^5*4*a_1*a_2*S_2_intgral/((R_rho^2 - 2*Z^2)*V_1*V_2);
    S_3_factor = 8*pi*R^5*4*a_1*a_2*S_3_intgral/(3*R_rho^2*V_1*V_2);
    S_4_factor = 8*pi*R^5*4*a_1*a_2*S_4_intgral/(3*Z*R_rho*V_1*V_2);

end

function out = S_1_integrand(k_rho, a_1, a_2, t_1, t_2, R_rho, Z)
    z_val = k_rho.*R_rho;
    z_val_1 = k_rho.*a_1;
    z_val_2 = k_rho.*a_2;

    [J_1] = besselj(1,z_val);
    [J_1_1] = besselj(1,z_val_1);
    [J_1_2] = besselj(1,z_val_2);

    k_z_intgral_num_1 = S_1_k_z_intgral_num(abs(Z + (t_1 + t_2)/2), k_rho);
    k_z_intgral_num_2 = S_1_k_z_intgral_num(abs(Z + (t_1 - t_2)/2), k_rho);
    k_z_intgral_num_3 = S_1_k_z_intgral_num(abs(Z - (t_1 - t_2)/2), k_rho);
    k_z_intgral_num_4 = S_1_k_z_intgral_num(abs(Z - (t_1 + t_2)/2), k_rho);

    k_z_intgral_num = -k_z_intgral_num_1 + k_z_intgral_num_2 + k_z_intgral_num_3 - k_z_intgral_num_4;
    
    out = J_1.*J_1_1.*J_1_2.*k_z_intgral_num./(k_rho.^3);

end

function out = S_1_k_z_intgral_num(alpha_p, beta_p)
    out = pi*(-alpha_p.*beta_p - exp(-alpha_p.*beta_p))/8;
end

function out = S_2_integrand(k_rho, a_1, a_2, t_1, t_2, R_rho, Z)

    z_val = k_rho.*R_rho;
    z_val_1 = k_rho.*a_1;
    z_val_2 = k_rho.*a_2;
    
    [J_0] = besselj(0,z_val);
    [J_1_1] = besselj(1,z_val_1);
    [J_1_2] = besselj(1,z_val_2);
    
    k_z_intgral_num_1 = S_2_k_z_intgral_num(abs(Z + (t_1 + t_2)/2), k_rho);
    k_z_intgral_num_2 = S_2_k_z_intgral_num(abs(Z + (t_1 - t_2)/2), k_rho);
    k_z_intgral_num_3 = S_2_k_z_intgral_num(abs(Z - (t_1 - t_2)/2), k_rho);
    k_z_intgral_num_4 = S_2_k_z_intgral_num(abs(Z - (t_1 + t_2)/2), k_rho);
    
    k_z_intgral_num = -k_z_intgral_num_1 + k_z_intgral_num_2 + k_z_intgral_num_3 - k_z_intgral_num_4;
    
    out = J_0.*J_1_1.*J_1_2.*k_z_intgral_num./(k_rho.^2);

end

function out = S_2_k_z_intgral_num(alpha_p, beta_p)
    out = pi*exp(-alpha_p.*beta_p)/8;
end

function out = S_3_integrand(k_rho, a_1, a_2, t_1, t_2, R_rho, Z)

    z_val = k_rho.*R_rho;
    z_val_1 = k_rho.*a_1;
    z_val_2 = k_rho.*a_2;
    
    [J_2] = besselj(2,z_val);
    [J_1_1] = besselj(1,z_val_1);
    [J_1_2] = besselj(1,z_val_2);
    
    k_z_intgral_num_1 = S_1_k_z_intgral_num(abs(Z + (t_1 + t_2)/2), k_rho);
    k_z_intgral_num_2 = S_1_k_z_intgral_num(abs(Z + (t_1 - t_2)/2), k_rho);
    k_z_intgral_num_3 = S_1_k_z_intgral_num(abs(Z - (t_1 - t_2)/2), k_rho);
    k_z_intgral_num_4 = S_1_k_z_intgral_num(abs(Z - (t_1 + t_2)/2), k_rho);
    
    k_z_intgral_num = -k_z_intgral_num_1 + k_z_intgral_num_2 + k_z_intgral_num_3 - k_z_intgral_num_4;
    
    out = J_2.*J_1_1.*J_1_2.*k_z_intgral_num./(k_rho.^2);

end

function out = S_4_integrand(k_rho, a_1, a_2, t_1, t_2, R_rho, Z)

    z_val = k_rho.*R_rho;
    z_val_1 = k_rho.*a_1;
    z_val_2 = k_rho.*a_2;
    
    [J_1] = besselj(1,z_val);
    [J_1_1] = besselj(1,z_val_1);
    [J_1_2] = besselj(1,z_val_2);
    
    k_z_intgral_num_1 = sign(Z + (t_1 + t_2)/2).*S_4_k_z_intgral_num(abs(Z + (t_1 + t_2)/2), k_rho);
    k_z_intgral_num_2 = sign(Z + (t_1 - t_2)/2).*S_4_k_z_intgral_num(abs(Z + (t_1 - t_2)/2), k_rho);
    k_z_intgral_num_3 = sign(Z - (t_1 - t_2)/2).*S_4_k_z_intgral_num(abs(Z - (t_1 - t_2)/2), k_rho);
    k_z_intgral_num_4 = sign(Z - (t_1 + t_2)/2).*S_4_k_z_intgral_num(abs(Z - (t_1 + t_2)/2), k_rho);
    
    k_z_intgral_num = -k_z_intgral_num_1 + k_z_intgral_num_2 + k_z_intgral_num_3 - k_z_intgral_num_4;
    
    out = J_1.*J_1_1.*J_1_2.*k_z_intgral_num./(k_rho.^2);

end

function out = S_4_k_z_intgral_num(alpha_p, beta_p)
    out = pi*(1 - exp(-alpha_p.*beta_p))/8;
end

%% References
% [1] Kalezhi, J. & Miles, J.J., 2011. An Energy Barrier Model for Write
% Errors in Exchange-Spring Patterned Media. Magnetics, IEEE Transactions
% on, 47(10), pp.2540-2543.
%
% [2] Kalezhi, J., Belle, B.D. & Miles, J.J., 2010. Dependence of
% Write-Window on Write Error Rates in Bit Patterned Media. Magnetics, IEEE
% Transactions on, 46(10), pp.3752-3759.
%
% [3] Kalezhi, J., Miles, J.J. & Belle, B.D., 2009. Dependence of Switching
% Fields on Island Shape in Bit Patterned Media. Magnetics, IEEE
% Transactions on, 45(10), pp.3531-3534.
