function F = equilibrio(X)

    theta1 = X(1);
    theta2 = X(2);
    phi1 = X(3);
    phi2 = X(4);
    

    F(1) = -(1960*sin(theta1) + 39200*cos(phi1)^2*sin(theta1) + 39200*cos(phi2)^2*sin(theta1) + 39200*cos(theta2)^2*sin(theta1) - (11466*cos(phi1)*cos(phi2)^2*cos(theta1))/5 - 78400*cos(phi1)^2*cos(phi2)^2*sin(theta1) - 39200*cos(phi1)^2*cos(theta2)^2*sin(theta1) - 39200*cos(phi2)^2*cos(theta2)^2*sin(theta1) + (11466*cos(phi1)*cos(phi2)^2*cos(theta1)*cos(theta2)^2)/5 - (11466*cos(phi2)*cos(theta1)*sin(phi1)*sin(phi2))/5 + (11466*cos(phi2)*cos(theta2)*sin(theta1)*sin(theta2))/5 + 78400*cos(phi1)^2*cos(phi2)^2*cos(theta2)^2*sin(theta1) + (11466*cos(phi2)*cos(theta1)*cos(theta2)^2*sin(phi1)*sin(phi2))/5 - 39200*cos(phi1)*cos(phi2)*cos(theta1)*cos(theta2)*sin(theta2) - 78400*cos(phi1)*cos(phi2)*sin(phi1)*sin(phi2)*sin(theta1) - 39200*cos(theta1)*cos(theta2)*sin(phi1)*sin(phi2)*sin(theta2) + 78400*cos(phi1)*cos(phi2)*cos(theta2)^2*sin(phi1)*sin(phi2)*sin(theta1))/(100*(2000*cos(phi1)^2 + 2000*cos(phi2)^2 + 2000*cos(theta1)^2 + 2000*cos(theta2)^2 - 4000*cos(phi1)^2*cos(phi2)^2 - 2000*cos(phi1)^2*cos(theta1)^2 - 2000*cos(phi1)^2*cos(theta2)^2 - 2000*cos(phi2)^2*cos(theta1)^2 - 2000*cos(phi2)^2*cos(theta2)^2 - 4000*cos(theta1)^2*cos(theta2)^2 + 2000*cos(phi1)^2*cos(theta1)^2*cos(theta2)^2 + 2000*cos(phi2)^2*cos(theta1)^2*cos(theta2)^2 + 4000*cos(phi1)^2*cos(phi2)^2*cos(theta1)^2 + 4000*cos(phi1)^2*cos(phi2)^2*cos(theta2)^2 - 4000*cos(phi1)^2*cos(phi2)^2*cos(theta1)^2*cos(theta2)^2 - 4000*cos(phi1)*cos(phi2)*sin(phi1)*sin(phi2) + 4000*cos(phi1)*cos(phi2)*cos(theta1)^2*sin(phi1)*sin(phi2) + 4000*cos(phi1)*cos(phi2)*cos(theta2)^2*sin(phi1)*sin(phi2) - 4000*cos(phi1)*cos(phi2)*cos(theta1)^2*cos(theta2)^2*sin(phi1)*sin(phi2) - 4000*cos(phi1)*cos(phi2)*cos(theta1)*cos(theta2)*sin(theta1)*sin(theta2) - 4000*cos(theta1)*cos(theta2)*sin(phi1)*sin(phi2)*sin(theta1)*sin(theta2) + 100));

    F(2) = -(82320000*cos(theta1)^2*sin(theta2) - 4815720*cos(phi2)*cos(theta2) + 4586400*cos(phi2)*cos(theta1)^2*cos(theta2) + 4586400*cos(phi1)*cos(phi2)^2*cos(theta1)*sin(theta1)*sin(theta2) - 82320000*cos(phi1)*cos(phi2)*cos(theta1)*cos(theta2)*sin(theta1) - 82320000*cos(theta1)*cos(theta2)*sin(phi1)*sin(phi2)*sin(theta1) + 4586400*cos(phi2)*cos(theta1)*sin(phi1)*sin(phi2)*sin(theta1)*sin(theta2))/(100000*(2000*cos(phi1)^2 + 2000*cos(phi2)^2 + 2000*cos(theta1)^2 + 2000*cos(theta2)^2 - 4000*cos(phi1)^2*cos(phi2)^2 - 2000*cos(phi1)^2*cos(theta1)^2 - 2000*cos(phi1)^2*cos(theta2)^2 - 2000*cos(phi2)^2*cos(theta1)^2 - 2000*cos(phi2)^2*cos(theta2)^2 - 4000*cos(theta1)^2*cos(theta2)^2 + 2000*cos(phi1)^2*cos(theta1)^2*cos(theta2)^2 + 2000*cos(phi2)^2*cos(theta1)^2*cos(theta2)^2 + 4000*cos(phi1)^2*cos(phi2)^2*cos(theta1)^2 + 4000*cos(phi1)^2*cos(phi2)^2*cos(theta2)^2 - 4000*cos(phi1)^2*cos(phi2)^2*cos(theta1)^2*cos(theta2)^2 - 4000*cos(phi1)*cos(phi2)*sin(phi1)*sin(phi2) + 4000*cos(phi1)*cos(phi2)*cos(theta1)^2*sin(phi1)*sin(phi2) + 4000*cos(phi1)*cos(phi2)*cos(theta2)^2*sin(phi1)*sin(phi2) - 4000*cos(phi1)*cos(phi2)*cos(theta1)^2*cos(theta2)^2*sin(phi1)*sin(phi2) - 4000*cos(phi1)*cos(phi2)*cos(theta1)*cos(theta2)*sin(theta1)*sin(theta2) - 4000*cos(theta1)*cos(theta2)*sin(phi1)*sin(phi2)*sin(theta1)*sin(theta2) + 100));

    F(3) = (39200*cos(phi1)*cos(theta1)^3*sin(phi1) - 39200*cos(phi2)*cos(theta1)^3*sin(phi2) + (11466*cos(phi2)^2*sin(phi1)*sin(theta1))/5 - 39200*cos(phi1)*cos(theta1)*sin(phi1) + 39200*cos(phi2)*cos(theta1)*sin(phi2) + 78400*cos(phi1)*cos(phi2)^2*cos(theta1)*sin(phi1) - 78400*cos(phi1)^2*cos(phi2)*cos(theta1)*sin(phi2) + 39200*cos(phi1)*cos(theta1)*cos(theta2)^2*sin(phi1) - 39200*cos(phi2)*cos(theta1)*cos(theta2)^2*sin(phi2) - 78400*cos(phi1)*cos(phi2)^2*cos(theta1)^3*sin(phi1) + 78400*cos(phi1)^2*cos(phi2)*cos(theta1)^3*sin(phi2) - 39200*cos(phi1)*cos(theta1)^3*cos(theta2)^2*sin(phi1) + 39200*cos(phi2)*cos(theta1)^3*cos(theta2)^2*sin(phi2) - (11466*cos(phi2)^2*cos(theta2)^2*sin(phi1)*sin(theta1))/5 - (11466*cos(phi1)*cos(phi2)*sin(phi2)*sin(theta1))/5 + (11466*cos(phi1)*cos(phi2)*cos(theta2)^2*sin(phi2)*sin(theta1))/5 - 78400*cos(phi1)*cos(phi2)^2*cos(theta1)*cos(theta2)^2*sin(phi1) + 78400*cos(phi1)^2*cos(phi2)*cos(theta1)*cos(theta2)^2*sin(phi2) + 78400*cos(phi1)*cos(phi2)^2*cos(theta1)^3*cos(theta2)^2*sin(phi1) - 78400*cos(phi1)^2*cos(phi2)*cos(theta1)^3*cos(theta2)^2*sin(phi2) - 39200*cos(phi1)*cos(theta1)^2*cos(theta2)*sin(phi2)*sin(theta1)*sin(theta2) + 39200*cos(phi2)*cos(theta1)^2*cos(theta2)*sin(phi1)*sin(theta1)*sin(theta2))/(100*(cos(theta1)^2 - 1)*(2000*cos(phi1)^2 + 2000*cos(phi2)^2 + 2000*cos(theta1)^2 + 2000*cos(theta2)^2 - 4000*cos(phi1)^2*cos(phi2)^2 - 2000*cos(phi1)^2*cos(theta1)^2 - 2000*cos(phi1)^2*cos(theta2)^2 - 2000*cos(phi2)^2*cos(theta1)^2 - 2000*cos(phi2)^2*cos(theta2)^2 - 4000*cos(theta1)^2*cos(theta2)^2 + 2000*cos(phi1)^2*cos(theta1)^2*cos(theta2)^2 + 2000*cos(phi2)^2*cos(theta1)^2*cos(theta2)^2 + 4000*cos(phi1)^2*cos(phi2)^2*cos(theta1)^2 + 4000*cos(phi1)^2*cos(phi2)^2*cos(theta2)^2 - 4000*cos(phi1)^2*cos(phi2)^2*cos(theta1)^2*cos(theta2)^2 - 4000*cos(phi1)*cos(phi2)*sin(phi1)*sin(phi2) + 4000*cos(phi1)*cos(phi2)*cos(theta1)^2*sin(phi1)*sin(phi2) + 4000*cos(phi1)*cos(phi2)*cos(theta2)^2*sin(phi1)*sin(phi2) - 4000*cos(phi1)*cos(phi2)*cos(theta1)^2*cos(theta2)^2*sin(phi1)*sin(phi2) - 4000*cos(phi1)*cos(phi2)*cos(theta1)*cos(theta2)*sin(theta1)*sin(theta2) - 4000*cos(theta1)*cos(theta2)*sin(phi1)*sin(phi2)*sin(theta1)*sin(theta2) + 100));

    F(4) = (229320*sin(phi2)*sin(theta2) + 4586400*cos(phi1)^2*sin(phi2)*sin(theta2) + 4586400*cos(theta1)^2*sin(phi2)*sin(theta2) + 4586400*cos(theta2)^2*sin(phi2)*sin(theta2) + 9172800*cos(theta1)*cos(theta2)^3*sin(phi1)*sin(theta1) - 4586400*cos(phi1)^2*cos(theta1)^2*sin(phi2)*sin(theta2) - 4586400*cos(phi1)^2*cos(theta2)^2*sin(phi2)*sin(theta2) - 9172800*cos(theta1)^2*cos(theta2)^2*sin(phi2)*sin(theta2) - 4586400*cos(phi1)*cos(phi2)*sin(phi1)*sin(theta2) - 9172800*cos(theta1)*cos(theta2)*sin(phi1)*sin(theta1) + 4586400*cos(phi1)*cos(phi2)*cos(theta1)^2*sin(phi1)*sin(theta2) + 4586400*cos(phi1)*cos(phi2)*cos(theta2)^2*sin(phi1)*sin(theta2) + 4586400*cos(phi2)^2*cos(theta1)*cos(theta2)*sin(phi1)*sin(theta1) - 4586400*cos(phi2)^2*cos(theta1)*cos(theta2)^3*sin(phi1)*sin(theta1) + 82320000*cos(phi1)*cos(theta1)*sin(phi2)*sin(theta1)*sin(theta2) - 82320000*cos(phi2)*cos(theta1)*sin(phi1)*sin(theta1)*sin(theta2) + 4586400*cos(phi1)^2*cos(theta1)^2*cos(theta2)^2*sin(phi2)*sin(theta2) - 4586400*cos(phi1)*cos(phi2)*cos(theta1)^2*cos(theta2)^2*sin(phi1)*sin(theta2) - 4586400*cos(phi1)*cos(phi2)*cos(theta1)*cos(theta2)*sin(phi2)*sin(theta1) + 4586400*cos(phi1)*cos(phi2)*cos(theta1)*cos(theta2)^3*sin(phi2)*sin(theta1))/(100000*(cos(theta2)^2 - 1)*(2000*cos(phi1)^2 + 2000*cos(phi2)^2 + 2000*cos(theta1)^2 + 2000*cos(theta2)^2 - 4000*cos(phi1)^2*cos(phi2)^2 - 2000*cos(phi1)^2*cos(theta1)^2 - 2000*cos(phi1)^2*cos(theta2)^2 - 2000*cos(phi2)^2*cos(theta1)^2 - 2000*cos(phi2)^2*cos(theta2)^2 - 4000*cos(theta1)^2*cos(theta2)^2 + 2000*cos(phi1)^2*cos(theta1)^2*cos(theta2)^2 + 2000*cos(phi2)^2*cos(theta1)^2*cos(theta2)^2 + 4000*cos(phi1)^2*cos(phi2)^2*cos(theta1)^2 + 4000*cos(phi1)^2*cos(phi2)^2*cos(theta2)^2 - 4000*cos(phi1)^2*cos(phi2)^2*cos(theta1)^2*cos(theta2)^2 - 4000*cos(phi1)*cos(phi2)*sin(phi1)*sin(phi2) + 4000*cos(phi1)*cos(phi2)*cos(theta1)^2*sin(phi1)*sin(phi2) + 4000*cos(phi1)*cos(phi2)*cos(theta2)^2*sin(phi1)*sin(phi2) - 4000*cos(phi1)*cos(phi2)*cos(theta1)^2*cos(theta2)^2*sin(phi1)*sin(phi2) - 4000*cos(phi1)*cos(phi2)*cos(theta1)*cos(theta2)*sin(theta1)*sin(theta2) - 4000*cos(theta1)*cos(theta2)*sin(phi1)*sin(phi2)*sin(theta1)*sin(theta2) + 100));

    %F(1) = ddtheta1_fun(theta1, theta2, phi1, phi2);
    %F(2) = ddtheta2_fun(theta1, theta2, phi1, phi2);
    %F(3) = ddphi1_fun(theta1, theta2, phi1, phi2);
    %F(4) = ddphi2_fun(theta1, theta2, phi1, phi2);

end

options = optimoptions('fsolve', ...
    'FunctionTolerance', 1e-14, ...
    'StepTolerance', 1e-14, ...
    'OptimalityTolerance', 1e-14, ...
    'Display', 'iter');

angulos0 = [0.013376202173787;0.014044926426593;0;0];
fun = @equilibrio;

[X, fval, exitflag, output] = fsolve(fun, angulos0, options);

format long
disp('Solução:')
disp(X)
disp('Valor da função na solução:')
disp(fval)



%Solução:
%   0.054876785277665
%   0.057614700529080
%  -0.007121303499816
%  -0.007121303499816
