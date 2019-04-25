clear all
clc

L_ac = 1e-3;
L_dc = 5e-3;
w_ac = 100;
T_ac = (2*pi)/w_ac;
sig = 0.01;
tau = 1e-5;
R_load = 0.2865;

del_t = T_ac/10000;
t_end = 10*T_ac;
t(1) = 0;
k = 1;

e_a(1) = 100;
e_b(1) = -50;
e_c(1) = -50;

V_dc_p(1) = 0;
i_dc(1) = 0;

i_a(1) = 0;
i_b(1) = 0;
i_c(1) = 0;

i_d1(1) = 0;
i_d3(1) = 0;
i_d5(1) = 0;

V_ag(1) = 0;
V_bg(1) = 0;
V_cg(1) = 0;

while t(k) < t_end
    % diode currents logic
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
    
    % output current
    i_dc(k+1) = i_d1(k+1) + i_d3(k+1) + i_d5(k+1);
    
    % Logic for V_ag, V_bg and V_cg
    if i_a(k) > sig
        V_ag(k) = V_dc_p(k);
    elseif i_a(k) < -sig
        V_ag(k) = 0;
    else
        V_ag(k) = ((V_dc_p(k)*i_a(k))/(2*sig)) + (V_dc_p(k)/2);
    end
    
    if i_b(k) > sig
        V_bg(k) = V_dc_p(k);
    elseif i_b(k) < -sig
        V_bg(k) = 0;
    else
        V_bg(k) = ((V_dc_p(k)*i_b(k))/(2*sig)) + (V_dc_p(k)/2);
    end
    
    if i_c(k) > sig
        V_cg(k) = V_dc_p(k);
    elseif i_c(k) < -sig
        V_cg(k) = 0;
    else
        V_cg(k) = ((V_dc_p(k)*i_c(k))/(2*sig)) + (V_dc_p(k)/2);
    end
    
    % Va, Vb, Vc calculations
    V_a(k) = (2/3)*V_ag(k) - (1/3)*V_bg(k) - (1/3)*V_cg(k);
    V_b(k) = (2/3)*V_bg(k) - (1/3)*V_ag(k) - (1/3)*V_cg(k);
    V_c(k) = (2/3)*V_cg(k) - (1/3)*V_ag(k) - (1/3)*V_bg(k);
    
    e_a(k) = 100*cos(100*t(k));
    e_b(k) = 100*cos((100*t(k)) - (2*pi/3));
    e_c(k) = 100*cos((100*t(k)) + (2*pi/3));
    
    % Backward Euler
    i_a(k+1) = i_a(k) + del_t*((e_a(k) - V_a(k))/L_ac);
    i_b(k+1) = i_b(k) + del_t*((e_b(k) - V_b(k))/L_ac);
    i_c(k+1) = i_c(k) + del_t*((e_b(k) - V_c(k))/L_ac);
    
    %V_dc_p(k+1) = (1/((del_t/tau)+1)) * (V_dc_p(k) + ((L_dc*(i_dc(k+1)-i_dc(k)))/tau) + ((del_t*i_dc(k+1)*R_load)/tau));
    V_dc_p(k+1) = R_load*i_dc(k+1);
    
    t(k+1) = t(k) + del_t;
    k = k+1;
end

subplot(3,1,1)
plot(t,i_a)
title('i_a')

subplot(3,1,2)
plot(t,V_dc_p)
title("V_dc_prime")

subplot(3,1,3)
plot(t,i_dc)
title("i_dc")