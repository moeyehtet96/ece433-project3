clear all
clc

V_ac = 100;
w_ac = 100;
L_ac = 1e-3;

V_dc_12 = 9*sqrt(3)*V_ac/(4*pi);
V_dc_23 = 9*V_ac/(4*pi);

V_dc_mode2 = V_dc_23:0.00005:V_dc_12;
V_dc_mode1 = V_dc_12:0.00005:3*sqrt(3)*V_ac/pi;

i_dc_mode1 = (pi/(3*w_ac*L_ac)).*((3*sqrt(3)*V_ac/pi) - V_dc_mode1);
alpha = acos((2*pi.*V_dc_mode2)/(9*V_ac)) - pi/6;
i_dc_mode2 = ((sqrt(3)*V_ac)/(2*w_ac*L_ac)).*sin(alpha + pi/6);

i_dc = [i_dc_mode1, i_dc_mode2];
V_dc = [V_dc_mode1, V_dc_mode2];

p_out = V_dc_mode2.*i_dc_mode2;
r_load = V_dc_mode2./i_dc_mode2;

[p_max, p_max_index] = max(p_out)
r_p_max = r_load(p_max_index)

theta_ac = 0:5:360;
gamma = 45;
V_dc_b = ((3*sqrt(3)*V_ac)/(2*pi))*(cosd(gamma) + 1);
i_dc_b = (pi/(3*w_ac*L_ac))*(((3*sqrt(3)*V_ac)/pi) - V_dc_b);
i_ac_b = zeros(length(theta_ac));

pos = find(theta_ac <= 60 | theta_ac >= 300+gamma);
neg = find(theta_ac >= 120+gamma & theta_ac <= 240);
neg_slope_1 = find(theta_ac > 60 & theta_ac < 60+gamma);
neg_slope_2 = find(theta_ac > 120 & theta_ac < 120+gamma);
pos_slope_1 = find(theta_ac > 240 & theta_ac < 240+gamma);
pos_slope_2 = find(theta_ac > 300 & theta_ac < 300+gamma);
zero = find((theta_ac >= 60+gamma & theta_ac <= 120) | (theta_ac >= 240+gamma & theta_ac <= 300));

i_ac_b(pos) = i_dc_b;
i_ac_b(neg) = -i_dc_b;
i_ac_b(neg_slope_1) = -(i_dc_b/45).*(theta_ac(neg_slope_1)-60) + i_dc_b;
i_ac_b(neg_slope_2) = -(i_dc_b/45).*(theta_ac(neg_slope_2)-120);
i_ac_b(pos_slope_1) = (i_dc_b/45).*(theta_ac(pos_slope_1)-240) - i_dc_b;
i_ac_b(pos_slope_2) = (i_dc_b/45).*(theta_ac(pos_slope_2)-300);
i_ac_b(zero) = 0;

alpha = 10;
V_dc_c = ((9*V_ac)/(2*pi))*cosd(alpha+30);
i_dc_c = ((sqrt(3)*V_ac)/(2*w_ac*L_ac))*sind(alpha+30);
i_ac_c = zeros(length(theta_ac));

pos_c = find(theta_ac >= alpha & theta_ac <= 60+alpha);
neg_c = find(theta_ac >= 180+alpha & theta_ac <= 240+alpha);
neg_slope_c = find(theta_ac > 60+alpha & theta_ac < 180+alpha);
pos_slope_c_1 = find(theta_ac < alpha);
pos_slope_c_2 = find(theta_ac > 240+alpha);

i_ac_c(pos_c) = i_dc_c;
i_ac_c(neg_c) = -i_dc_c;
i_ac_c(neg_slope_c) = -(i_dc_c/60)*(theta_ac(neg_slope_c)-60-alpha) + i_dc_c;
i_ac_c(pos_slope_c_1) = (i_dc_c/60)*(theta_ac(pos_slope_c_1)+60-alpha);
i_ac_c(pos_slope_c_2) = (i_dc_c/60)*(theta_ac(pos_slope_c_2)-240-alpha) - i_dc_c;

n = find(V_dc == V_dc_12);
m = find(V_dc == V_dc_23);
[V_dc_max,o] = max(V_dc);
p = find(theta_ac == 60);
q = find(theta_ac == 60+alpha);

figure;
plot(i_dc_mode1,V_dc_mode1)
hold on
plot(i_dc_mode2,V_dc_mode2)
legend('Mode 1', 'Mode 2')
xlabel('i_d_c (A)')
ylabel('V_d_c (V)')
title('Output DC voltage vs Output DC current')
text(i_dc(o),V_dc(o),[num2str(V_dc_max),' V'])
text(i_dc(n),V_dc(n),[num2str(V_dc_12), ' V'])
text(i_dc(m),V_dc(m),[num2str(V_dc_23), ' V'])
saveas(gcf,'v-vs-i','tiffn')

figure;
plot(r_load,p_out)
xlabel('r_l_o_a_d (\Omega)')
ylabel('p_o_u_t (W)')
title('Power Output vs Load Resistance')
text(r_p_max,p_max,['R_l_o_a_d at maximum power = ', num2str(r_p_max), ' \Omega'])
saveas(gcf,'p-vs-r','tiffn')

figure;
plot(theta_ac, i_ac_b)
xlabel('\theta_a_c (deg)')
ylabel('i_a (A)')
title('Phase a current vs \theta_a_c for Step (b)')
text(theta_ac(p),i_ac_b(p),[num2str(i_ac_b(p)),' A'])
saveas(gcf,'ia-vs-theta-b','tiffn')

figure;
plot(theta_ac, i_ac_c)
xlabel('\theta_a_c (deg)')
ylabel('i_a (A)')
title('Phase a current vs \theta_a_c for Step (c)')
text(theta_ac(q),i_ac_c(p),[num2str(i_ac_c(q)),' A'])
saveas(gcf,'ia-vs-theta-c','tiffn')
        