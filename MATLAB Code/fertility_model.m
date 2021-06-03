% Numerical solution to the ODE model in Chapter 4.
function sol = fertility_model(tf)



% tf = final time

% Parameters
rho=2;
b=2;
k=10;
omega=0.2;

% Initial conditions
P_0  = 0;
M_0  = 3;
Fm_0  = 0;
X_0 = 0;

initialvalues = [P_0; M_0; Fm_0; X_0];

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

plot(sol.x, sol.y(1,:), 'b')
plot(sol.x, sol.y(2,:), 'r')
plot(sol.x, sol.y(3,:), 'g', 'LineWidth', 2)
plot(sol.x, sol.y(4,:), 'k')
%plot(sol.x, sol.y(5,:), 'k--')

xlabel('Non-dimensionalised Time (\tau)')
ylabel('Non-dimensionalised Populations')
title('Population Dynamics over Time')
legend('P', 'M', 'Fm','X')
set(gca, 'FontSize', 12)

%--------------------------------------------------------------------------
% Model
    
function dydt = odemodel(t,y)
    % Variables
    P   = y(1);
    M  = y(2);
    Fm  = y(3);
    X= y(4);
    
    % ODEs
        
    dP= ((rho*(k-P-Fm-X)*(k-P-M)) / (k-P)) -(b/k)*P*(P+Fm)-omega*P;
    
    dM= (1/2)*b*Fm -(1/2)*(b/k)*(P+Fm)*M;
    
    dFm= ( ( rho*(k-P-Fm-X)*M)/(k-P)) - b*Fm - omega*Fm -(1/2)* (b/k)*(P+Fm)*Fm;
    
    dX = omega*(k-X) -(1/2)*(b/k)*(P+Fm)*X;
    
    dydt = [dP; dM; dFm;dX];

end

%--------------------------------------------------------------------------
end