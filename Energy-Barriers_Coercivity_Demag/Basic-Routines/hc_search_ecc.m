function [h tprob xo x] = hc_search_ecc(theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing, temp, kb, duration, attfreq)


    xo = log(2)/(1*duration*attfreq);
    if (temp > 0)
        tempcheck = 1;
        disp('Temperature > 0 K')
    else (temp == 0)
        tempcheck = 2;
        disp('Temperature == 0 K')
    end
    % Set options, to display iteration, change 'off' to 'iter'
    options = optimset('Display','off','TolFun',1e-50, 'TolX',1e-15);
        switch tempcheck
            case 1 
                disp('Calculating coercivity')
                ho = -1e-6;
                

                [h, f, exitflag] = fzero(@(x)switchrateIntegrand_function_ecc(x, theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing, temp, kb,xo), ho, options);

                x =f;
                hp_h_vec = h;
                hp_s_vec = h;
                tprob =switchrateIntegrand_ecc(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing, temp, kb);
            case 2
                disp('Calculating coercivity')
                h_min = -2; 
                h_max = 0;
                h_mid = -0.5;

                
                [h, f, exitflag] = fsolve(@(x)h_switch_search_function(x, theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing), h_mid, options);
                hp_h_vec = h;
                hp_s_vec = h;
                x =f;    
                [ebarrier1 ebarrier2]= energy_barrier_check(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, j_xc, k_h, k_s, m_h, m_s,v_h,v_s, area,t_h, t_s, a,spacing);
                tprob = 0;
            otherwise
                msg = 'Negative temperature. Temperature must equal to or greater than zero';
                error(msg)
        end
end

function out =switchrateIntegrand_function_ecc(h, theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing, temp, kb,xo)
hp_h_vec = h;
hp_s_vec = h;

tprob = switchrateIntegrand_ecc(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing, temp, kb);
out = tprob - xo;

end


function tprob =switchrateIntegrand_ecc(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing, temp, kb)
theta_H_h_vec = theta_H_h_vec.*ones(size(hp_h_vec));
theta_H_s_vec = theta_H_s_vec.*ones(size(hp_s_vec));

[eb1 eb2] = scaledbarrier_ecc(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing, temp, kb);

eb1_check = isreal(eb1);
eb1 = eb1.*eb1_check;

eb2_check = isreal(eb2);
eb2 = eb2.*eb2_check;

tprob = exp(-eb1)+ exp(-eb2);
end

function [eb1 eb2] = scaledbarrier_ecc(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing, temp, kb)
[ebarrier1, ebarrier2]= energy_barrier_check(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, j_xc, k_h, k_s, m_h, m_s,v_h,v_s, area,t_h, t_s, a,spacing);
eb1 = 2*k_h.*v_h*(1e-27).*ebarrier1/(kb*temp);
eb2 = 2*k_h.*v_h*(1e-27).*ebarrier2/(kb*temp);
end

function out = h_switch_search_function(h, theta_H_h_vec, theta_H_s_vec, k_h, k_s, m_h, m_s, j_xc, area, v_h, v_s, t_h, t_s, a, spacing)
hp_h_vec = h;
hp_s_vec = h;

[ebarrier1 ebarrier2]= energy_barrier_check(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, j_xc, k_h, k_s, m_h, m_s,v_h,v_s, area,t_h, t_s, a,spacing);
out = ebarrier1 -(sign(ebarrier1)-1).*h;

end
