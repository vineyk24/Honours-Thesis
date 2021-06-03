% The first model that we use to implement equations 1-5 in Chapter 2.
% Non-dimensionalised versions are implemented in subsequent chapters after
% the technique is known. 
function sol = numerics(tf)


% tf = final time

% Parameters
rho=10;
b=2;
delta=2;

% Initial conditions
S_0   =40;
G_0  = 25;
M_0  = 15;
F_0 = 80;
Fm_0  = 0;

initialvalues = [S_0; G_0; M_0; F_0; Fm_0];

%--------------------------------------------------------------------------

% Simulator
sol = ode45(@odemodel,[0 tf], initialvalues);

%--------------------------------------------------------------------------
% Plot figure

fontsize = 12;

% Plot
figure
set(gca, 'FontSize', fontsize)
hold on

plot(sol.x, sol.y(1,:), 'g--')
plot(sol.x, sol.y(2,:), 'b')
plot(sol.x, sol.y(3,:), 'r', 'LineWidth', 2)
plot(sol.x, sol.y(4,:), 'k')
plot(sol.x, sol.y(5,:), 'k--')

xlabel('Time (yr)')
ylabel('Population')
title('Population Dynamics over Time')
legend('S', 'G', 'M', 'F', 'Fm')
set(gca, 'FontSize', 12)

%--------------------------------------------------------------------------
% Model
    
function dydt = odemodel(t,y)
    % Variables
    S   = y(1);
    G  = y(2);
    M  = y(3);
    F = y(4);
    Fm  = y(5);

    mu = ((b*(Fm+G))/ (S+G+M)) - delta;
    
   
    
    % 5 ODEs, as we have not non-dimensionalised the system yet.
    
    dS= -rho*F*S - delta*S + delta*G + b*G -mu*S;
    
    dG= rho*F*S- 2*delta*G -mu *G;
    
    dM= -delta *M + b* Fm - mu *M;
    
    dF = -rho*F*S-rho*F*M+ 2*b*Fm+ b*G + delta*G-delta*F-mu*F;
    
    dFm= rho*F*M-b*Fm-delta*Fm-mu*Fm;
    
    dydt = [dS; dG; dM; dF; dFm];

end

%--------------------------------------------------------------------------
end