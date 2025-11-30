clear ; close all ; clc
syms l1 m1 l2 m2 g real;
syms x dx ddx gamma dgamma ddgamma theta1 dtheta1 ddtheta1 phi1 dphi1 ddphi1 theta2 dtheta2 ddtheta2 phi2 dphi2 ddphi2;

rho   = sym(1.225);   
Cd    = sym(1.2);     
A     = sym(15.6); 
vwind = sym(0);  

%rT = [x*cos(gamma); x*sin(gamma); 0];
%rm1 = rT  + [cos(gamma); sin(gamma); 0](l1*sin(theta1)*cos(phi1)) + [-sin(gamma); cos(gamma); 0](l1*sin(theta1)sin(phi1)) + [0; 0; 1](-l1*cos(theta1));
%rm2 = rm1 + [cos(gamma); sin(gamma); 0](l2*sin(theta2)*cos(phi2)) + [-sin(gamma); cos(gamma); 0](l2*sin(theta2)sin(phi2)) + [0; 0; 1](-l2*cos(theta2));

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

R = 0.5*(simplify(76.8*((vrf2(1)^2))) + simplify(76.8*((vrf2(2)^2))) + simplify(76.8*((vrf2(3)^2))));
%R = 0;

%vrel1 = [ vwind - vrf1(1) ; -vrf1(2) ; -vrf1(3) ];
vrel2 = [ vwind - vrf2(1) ; -vrf2(2) ; -vrf2(3) ];


Fd = simplify(0.5 * rho * Cd * A * vrel2 * norm(vrel));

%vrel = [ vrf2(1) - vwind ; vrf2(2) ; vrf2(3) ];

%vrel_sq = simplify(vrel.' * vrel);

L = simplify(T - V);
s = [theta1; theta2; phi1; phi2];
ds= [dtheta1; dtheta2; dphi1; dphi2];

%d_v2d_dtheta1 = diff(vrf2, dtheta1);
%d_v2d_dtheta2 = diff(vrf2, dtheta2);
%d_v2d_dphi1 = diff(vrf2, dphi1);
%d_v2d_dphi2 = diff(vrf2, dphi2);

Q = sym(zeros(length(s),1));

%for i = 1:length(s)
%    Q(i) = (-1)*(diff(R, ds(i)));
%end
%Q = simplify(Q.');

Q(1) = simplify(Fd.' * d_v2d_dtheta1);
Q(2) = simplify(Fd.' * d_v2d_dtheta2);
Q(3) = simplify(Fd.' * d_v2d_dphi1);
Q(4) = simplify(Fd.' * d_v2d_dphi2);

%Q(1) = Fd(1) * d_v2d_dtheta1(1) + Fd(2) * d_v2d_dtheta1(2) + Fd(3) * d_v2d_dtheta1(3);
%Q(2) = Fd(1) * d_v2d_dtheta2(1) + Fd(2) * d_v2d_dtheta2(2) + Fd(3) * d_v2d_dtheta2(3);
%Q(3) = Fd(1) * d_v2d_dphi1(1) + Fd(2) * d_v2d_dphi1(2) + Fd(3) * d_v2d_dphi1(3);
%Q(4) = Fd(1) * d_v2d_dphi2(1) + Fd(2) * d_v2d_dphi2(2) + Fd(3) * d_v2d_dphi2(3);


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

ddtheta1_fun = matlabFunction(sol.ddtheta1, 'Vars', [x, dx, gamma, dgamma, theta1, dtheta1, theta2, dtheta2, phi1, dphi1, phi2, dphi2, l1, m1, m2, l2, g]);
ddtheta2_fun = matlabFunction(sol.ddtheta2, 'Vars', [x, dx, gamma, dgamma, theta1, dtheta1, theta2, dtheta2, phi1, dphi1, phi2, dphi2, l1, m1, m2, l2, g]);
ddphi1_fun   = matlabFunction(sol.ddphi1,   'Vars', [x, dx, gamma, dgamma, theta1, dtheta1, theta2, dtheta2, phi1, dphi1, phi2, dphi2, l1, m1, m2, l2, g]);
ddphi2_fun   = matlabFunction(sol.ddphi2,   'Vars', [x, dx, gamma, dgamma, theta1, dtheta1, theta2, dtheta2, phi1, dphi1, phi2, dphi2, l1, m1, m2, l2, g]);

function dydt = d(~, y, ddtheta1_fun, ddtheta2_fun, ddphi1_fun, ddphi2_fun)
    l1 = 50;
    l2 = 5;
    m1 = 100;
    m2 = 24000;
    g  = 9.8;
    
    x = y(1);
    dx = y(2);
    gamma = y(3);
    dgamma = y(4);
    theta1 = y(5);
    dtheta1 = y(6);
    theta2 = y(7);
    dtheta2 = y(8);
    phi1 = y(9);
    dphi1 = y(10);
    phi2 = y(11);
    dphi2 = y(12);
    
    ddtheta1 = ddtheta1_fun(x, dx, gamma, dgamma, theta1, dtheta1, theta2, dtheta2, phi1, dphi1, phi2, dphi2, l1, m1, m2, l2, g);
    ddtheta2 = ddtheta2_fun(x, dx, gamma, dgamma, theta1, dtheta1, theta2, dtheta2, phi1, dphi1, phi2, dphi2, l1, m1, m2, l2, g);
    ddphi1   =   ddphi1_fun(x, dx, gamma, dgamma, theta1, dtheta1, theta2, dtheta2, phi1, dphi1, phi2, dphi2, l1, m1, m2, l2, g);
    ddphi2   =   ddphi2_fun(x, dx, gamma, dgamma, theta1, dtheta1, theta2, dtheta2, phi1, dphi1, phi2, dphi2, l1, m1, m2, l2, g);
    
    dydt = [dx; 0; dgamma; 0; dtheta1; ddtheta1; dtheta2; ddtheta2; dphi1; ddphi1; dphi2; ddphi2];
end

y0 = [40;%x
      0;%dx
      0;%gamma
      0.1;%dgamma
      pi/180;%theta1
      0;%dtheta1
      pi/180;%theta2
      0;%dtheta2
      pi/180;%phi1
      0;%dphi1
      pi/180;%phi2
      0];%dphi2

tf      = 60;                       % Final time            [s]
fR      = 30;                       % Frame rate            [fps]
dt      = 1/fR;                     % Time resolution       [s]
tspan    = linspace(0,tf,tf*fR);     % Time                  [s]

%tspan = [0 100];

[t,y] = ode113(@(t,y) d(t, y, ddtheta1_fun, ddtheta2_fun, ddphi1_fun, ddphi2_fun), tspan, y0);
%y(isnan(y)) = 0;
max_mod_col = max(abs(y));

L1 = 50;
L2 = 5;

XT = y(:,1).*cos(y(:,3));
YT = y(:,1).*sin(y(:,3));
ZT = y(:,1)*0;

%rT = [x*cos(gamma); x*sin(gamma); 0];
%rm1 = rT  + [l1*sin(theta1)*cos(phi1); l1*sin(theta1)*sin(phi1); -l1*cos(theta1)];
%rm2 = rm1 + [l2*sin(theta2)*cos(phi2); l2*sin(theta2)*sin(phi2); -l2*cos(theta2)];

XP1 = XT + L1*cos(y(:,9)).*sin(y(:,5));
YP1 = YT + L1*sin(y(:,9)).*sin(y(:,5));
ZP1 = -L1*cos(y(:,5));

XP2 = XP1 + L2*cos(y(:,11)).*sin(y(:,7));
YP2 = YP1 + L2*sin(y(:,11)).*sin(y(:,7));
ZP2 = ZP1 - L2*cos(y(:,7));

%XP1 = XT + L1*cos(y(:,3)).*cos(y(:,9)).*sin(y(:,5)) - L1*sin(y(:,3)).*sin(y(:,9)).*sin(y(:,5));
%YP1 = YT + L1*cos(y(:,3)).*sin(y(:,9)).*sin(y(:,5)) + L1*cos(y(:,9)).*sin(y(:,3)).*sin(y(:,5));
%ZP1 = -L1*cos(y(:,5));

%XP2 = XP1 + L2*cos(y(:,3)).*cos(y(:,11)).*sin(y(:,7)) - L2*sin(y(:,3)).*sin(y(:,11)).*sin(y(:,7));
%YP2 = YP1 + L2*cos(y(:,3)).*sin(y(:,11)).*sin(y(:,7)) + L2*cos(y(:,11)).*sin(y(:,3)).*sin(y(:,7));
%ZP2 = ZP1 - L2*cos(y(:,7));

%y0 = [1        %x
%      2        %dx
%      3        %gamma
%      4        %dgamma
%      5;       %theta1
%      6;       %dtheta1
%      7;       %theta2
%      8;       %dtheta2
%      9;       %phi1
%      10;      %dphi1
%      11;      %phi2
%      12];     %dphi2

% Tamanho dos Eixos
min_x = min(XP1)-50;
max_x = max(XP1)+50;
min_y = min(YP1)-50;
max_y = max(YP1)+50;
min_z = min(ZP1)-50;
max_z = 100;
figure
%set(gcf,'Position',[50 50 1920 1080]) % 1080p
%set(gcf,'Position',[50 50 1280 720]) % 720p
%set(gcf,'Position',[50 50 854 480]) % 480p
set(gcf,'Position',[50 50 640 640]) % Low

v = VideoWriter('Guindaste.mp4','MPEG-4');
v.Quality = 100;
open(v);

hold on ; grid on ; axis equal
set(gca,'CameraPosition',[42.0101   30.8293   16.2256])
set(gca,'XLim',[min_x max_x])
set(gca,'YLim',[min_y max_y])
set(gca,'ZLim',[min_z max_z])
set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[])
for i = 1:length(XP1)
    
    cla
    % Linha Vertical em x0 y0
    plot3([0 0],[0 0],[0 -100],'r')
    % Trajetória do Troley 
    plot3(XT(1:i),YT(1:i),ZT(1:i),'b');
    % Trajetória do Pêndulo 1 
    plot3(XP1(1:i),YP1(1:i),ZP1(1:i),'b')
    % Trajetória do Pêndulo 2
    plot3(XP2(1:i),YP2(1:i),ZP2(1:i),'b')
    % Reta ao Troley
    plot3([0 XT(i)],[0 YT(i)],[0 ZT(i)],'r')
    % Reta do Troley ao Pêndulo 1
    plot3([XT(i) XP1(i)],[YT(i) YP1(i)],[ZT(i) ZP1(i)],'r')
    % Reta do Pêndulo 1 ao Pêndulo 2
    plot3([XP1(i) XP2(i)],[YP1(i) YP2(i)],[ZP1(i) ZP2(i)],'r')
    % Esfera do Pêndulo 1
    plot3(XP1(i),YP1(i),ZP1(i),'Marker','o','Color','k','MarkerFaceColor','r','MarkerSize',5);
    % Esfera do Pêndulo 2
    plot3(XP2(i),YP2(i),ZP2(i),'Marker','o','Color','k','MarkerFaceColor','b','MarkerSize',10);
    % Esfera do Troley
    plot3(XT(i),YT(i),ZT(i),'Marker','o','Color','k','MarkerFaceColor','r','MarkerSize',10);
    % Projeções do Pêndulo 2
    plot3(min_x*ones(1,i),YP2(1:i),ZP2(1:i),'g')
    plot3(XP2(1:i),min_y*ones(1,i),ZP2(1:i),'g')
    plot3(XP2(1:i),YP2(1:i),min_z*ones(1,i),'g')
    frame = getframe(gcf);
    writeVideo(v,frame);

end
close(v);
