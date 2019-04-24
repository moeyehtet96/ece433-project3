L_ac = 1e-3;
L_dc = 5e-3;
w_ac = 100;
f_ac = w_ac/(2*pi);
T_ac = 1/f_ac;
sigma = 0.01;
tau = 1e-5;

i_a(1) = 0;
i_b(1) = 0;
i_c(1) = 0;

i_d1(1) = 0;
i_d3(1) = 0;
i_d5(1) = 0;

i_dc(1) = 0;

v_dc_p(1) = 0;

del_t = T_ac/1000;
t_end = 10*T_ac;
t(1) = 0;
k = 1;

while t(k) < t_end
    % Logic for i_d1, i_d3, i_d5
    if i_a(k) > 0
        i_d1(k+1) = i_a(k);
    else
        i_d1(k+1) = 0;
    end
    
    if i_b(k) > 0
        i_d3(k+1) = i_b(k);
    else
        i_d3(k+1) = 0;
    end
    
    if i_c(k) > 0
        i_d5(k+1) = i_c(k);
    else
        i_d5(k+1) = 0;
    end
    
    % i_dc calculation
    i_dc(k+1) = i_d1(k+1) + i_d3(k+1) + i_d5(k+1);
    
    % Logic for i_a, i_b and i_c
    if i_a(k) > sigma
        v_ag(k) = v_dc_p(k);
    elseif i_a(k) < sigma
        v_ag(k) = 0;
    else 
        v_ag(k) = (v_dc_p(k)*i_a(k))/(2*sigma) + v_dc_p(k)/2; 
    end
    
    if i_b(k) > sigma
        v_bg(k) = v_dc_p(k);
    elseif i_a(k) < sigma
        v_bg(k) = 0;
    else 
        v_bg(k) = (v_dc_p(k)*i_b(k))/(2*sigma) + v_dc_p(k)/2;
    end
    
    if i_c(k) > sigma
        v_cg(k) = v_dc_p(k);
    elseif i_a(k) < sigma
        v_cg(k) = 0;
    else 
        v_cg(k) = (v_dc_p(k)*i_c(k))/(2*sigma) + v_dc_p(k)/2;
    end
    
    % e_a, e_b, e_c calculation
    e_a(k) = 100*cos(100*t(k));
    e_b(k) = 100*cos(100*t(k) - 2*pi/3);
    e_c(k) = 100*cos(100*t(k) + 2*pi/3);
    
    % v_a, v_b, v_c calculation
    v_a(k) = (2/3)*v_ag(k) - (1/3)*v_bg(k) - (1/3)*v_cg(k);
    v_b(k) = (2/3)*v_bg(k) - (1/3)*v_ag(k) - (1/3)*v_cg(k);
    v_c(k) = (2/3)*v_cg(k) - (1/3)*v_bg(k) - (1/3)*v_ag(k);
    
    % i_a, i_b, i_c calculation
    i_a(k+1) = i_a(k) + del_t*(e_a(k) - v_a(k))/L_ac;
    i_b(k+1) = i_b(k) + del_t*(e_b(k) - v_b(k))/L_ac;
    i_c(k+1) = i_c(k) + del_t*(e_c(k) - v_c(k))/L_ac;
    
    % v_dc_prime calculation
    v_dc_p(k+1) = (1/((del_t/tau) + 1)) * ((L_dc/tau)*(i_dc(k+1)-i_dc(k)) + (del_t/tau)*i_dc(k+1)*R_load + v_dc_p(k));
    
    % loop update
    t(k+1) = t(k) + del_t;
    k = k + 1;
end

plot(t,v_dc_p)