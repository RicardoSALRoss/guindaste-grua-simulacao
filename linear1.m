clear ; close all ; clc
syms l1 m1 l2 m2 g rho Cd A real;
syms x dx ddx gamma dgamma ddgamma theta1 dtheta1 ddtheta1 phi1 dphi1 ddphi1 theta2 dtheta2 ddtheta2 phi2 dphi2 ddphi2;

%rho   = sym(1.225);   
%Cd    = sym(1.2);     
%A     = sym(15.6); 
vwind = sym(10);  

rT = [x*cos(gamma); x*sin(gamma); 0];
rm1 = rT  + [l1*sin(theta1)*cos(phi1); l1*sin(theta1)*sin(phi1); -l1*cos(theta1)];
rm2 = rm1 + [l2*sin(theta2)*cos(phi2); l2*sin(theta2)*sin(phi2); -l2*cos(theta2)];

vr1 = diff(rm1, x)*dx + diff(rm1, gamma)*dgamma + diff(rm1, theta1)*dtheta1 + diff(rm1, phi1)*dphi1;
vrf1 = simplify(formula(vr1));

vr2 = diff(rm2, x)*dx + diff(rm2, gamma)*dgamma + diff(rm2, theta1)*dtheta1 + diff(rm2, phi1)*dphi1 + diff(rm2,theta2)*dtheta2 + diff(rm2,phi2)*dphi2;
vrf2 = simplify(formula(vr2));

T1 = simplify( 0.5 * m1 * (vrf1(1)^2 + vrf1(2)^2 + vrf1(3)^2 ));   
T2 = simplify( 0.5 * m2 * (vrf2(1)^2 + vrf2(2)^2 + vrf2(3)^2 ));
T = simplify(T1 + T2);

V1 = m1 * g * (-l1*cos(theta1));
V2 = m2 * g * (-l1*cos(theta1)-l2*cos(theta2));
V = simplify(V1 + V2);

vrel2 = [ vwind - vrf2(1) ; -vrf2(2) ; -vrf2(3) ];

Fd = 0.5 * rho * Cd * A * vrel2 * sqrt(vrel2(1)^2 + vrel2(2)^2 + vrel2(3)^2);

L = simplify(T - V);
s = [theta1; theta2; phi1; phi2];
ds= [dtheta1; dtheta2; dphi1; dphi2];

d_v2d_dtheta1 = diff(vrf2, dtheta1);
d_v2d_dtheta2 = diff(vrf2, dtheta2);
d_v2d_dphi1 = diff(vrf2, dphi1);
d_v2d_dphi2 = diff(vrf2, dphi2);

Q = sym(zeros(length(s),1));

Q(1) = simplify(Fd.' * d_v2d_dtheta1);
Q(2) = simplify(Fd.' * d_v2d_dtheta2);
Q(3) = simplify(Fd.' * d_v2d_dphi1);
Q(4) = simplify(Fd.' * d_v2d_dphi2);

function f = EL(s,ds,L,Q)
dds = str2sym( "d" + string(ds) );
f = sym(zeros(length(s),1));
for i = 1:length(s)
    d_s  = diff(L, s(i));
    d_ds = diff(L,ds(i));
    d_ds_dt = jacobian(d_ds, [s;ds]) * [ds;dds];
    f(i) = simplify(expand( solve(Q(i) == d_ds_dt - d_s, dds(i))));
end
end

VF = EL(s,ds,L,Q);

sol = solve([VF(1) == ddtheta1, VF(2) == ddtheta2, VF(3) == ddphi1, VF(4) == ddphi2], [ddtheta1, ddtheta2, ddphi1, ddphi2]);

% Soluções
ddtheta1_sol = simplify(sol.ddtheta1);
ddtheta2_sol = simplify(sol.ddtheta2);
ddphi1_sol = simplify(sol.ddphi1);
ddphi2_sol = simplify(sol.ddphi2);

syms_zero = [dx, ddx, dgamma, ddgamma, dtheta1, ddtheta1, dtheta2, ddtheta2, dphi1, ddphi1, dphi2, ddphi2];
zero_vals = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

ddtheta1_zero = subs(ddtheta1_sol, syms_zero, zero_vals);
ddtheta2_zero = subs(ddtheta2_sol, syms_zero, zero_vals);
ddphi1_zero = subs(ddphi1_sol, syms_zero, zero_vals);
ddphi2_zero = subs(ddphi2_sol, syms_zero, zero_vals);

% Simplificação
ddtheta1_zero = simplify(ddtheta1_zero);
ddtheta2_zero = simplify(ddtheta2_zero);
ddphi1_zero = simplify(ddphi1_zero);
ddphi2_zero = simplify(ddphi2_zero);

% Parâmetros físicos fornecidos
l1_val = 50;      % m
l2_val = 25;       % m
m1_val = 100;     % kg
m2_val = 2000;   % kg
g_val = 9.8;      % m/s²
A_val = 15.6;       % m²
Cd_val = 1.2;     % adimensional
rho_val = 1.225;  % kg/m³

% Substituir todos os parâmetros físicos
ddtheta1_param = subs(ddtheta1_zero, [gamma, l1, l2, m1, m2, g, A, Cd, rho], ...
                     [pi/4, l1_val, l2_val, m1_val, m2_val, g_val, A_val, Cd_val, rho_val])

ddtheta2_param = subs(ddtheta2_zero, [gamma, l1, l2, m1, m2, g, A, Cd, rho], ...
                     [pi/4, l1_val, l2_val, m1_val, m2_val, g_val, A_val, Cd_val, rho_val])

ddphi1_param = subs(ddphi1_zero, [gamma, l1, l2, m1, m2, g, A, Cd, rho], ...
                   [pi/4, l1_val, l2_val, m1_val, m2_val, g_val, A_val, Cd_val, rho_val])

ddphi2_param = subs(ddphi2_zero, [gamma, l1, l2, m1, m2, g, A, Cd, rho], ...
                   [pi/4, l1_val, l2_val, m1_val, m2_val, g_val, A_val, Cd_val, rho_val])

%ddtheta1_param = subs(ddtheta1_zero, [x, phi1, phi2, gamma, l1, l2, m1, m2, g, A, Cd, rho], ...
%                     [40, 0, 0, 0, l1_val, l2_val, m1_val, m2_val, g_val, A_val, Cd_val, rho_val])

%ddtheta2_param = subs(ddtheta2_zero, [x, phi1, phi2, gamma, l1, l2, m1, m2, g, A, Cd, rho], ...
%                     [40, 0, 0, 0, l1_val, l2_val, m1_val, m2_val, g_val, A_val, Cd_val, rho_val])

%ddphi1_param = subs(ddphi1_zero, [x, phi1, phi2, gamma, l1, l2, m1, m2, g, A, Cd, rho], ...
%                   [40, 0, 0, 0, l1_val, l2_val, m1_val, m2_val, g_val, A_val, Cd_val, rho_val])

%ddphi2_param = subs(ddphi2_zero, [x, phi1, phi2, gamma, l1, l2, m1, m2, g, A, Cd, rho], ...
%                   [40, 0, 0, 0, l1_val, l2_val, m1_val, m2_val, g_val, A_val, Cd_val, rho_val])

symvar(ddtheta1_param)
symvar(ddtheta2_param)
symvar(ddphi1_param)
symvar(ddphi2_param)
