% The non-dimensionalised model of Chapter 2.
function sol = model_non_dimensionalised(tf)


% tf = final time

% Parameters (kappa)
kappa=1.8;

% Initial conditions
P_0  = 0;
M_0  = 0.5;
Fm_0  = 0;

initialvalues = [P_0; M_0; Fm_0];

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
plot(sol.x, sol.y(3,:), 'color',[0.9290, 0.6940, 0.1250], 'LineWidth', 2)
%plot(sol.x, sol.y(4,:), 'k')
%plot(sol.x, sol.y(5,:), 'k--')

xlabel('Non-dimensionalised Time (\tau)')
ylabel('Non-dimensionalised Populations')
title('Population Dynamics over Time')
legend('P', 'M', 'Fm')
set(gca, 'FontSize', 12)

%--------------------------------------------------------------------------
% Model
    
function dydt = odemodel(t,y)
    % Variables
    P   = y(1);
    M  = y(2);
    Fm  = y(3);
    
    
    % Simulate the system of ODEs
        
    dP= ((kappa*(1-P-Fm)*(1-P-M)) / (1-P)) -P*(P+Fm);
    
    dM= (1/2)*Fm -(1/2)*(P+Fm)*M;
    
    dFm= ( ( kappa*(1-P-Fm)*M)/(1-P)) - Fm -(1/2)*(P+Fm)*Fm;
    
    dydt = [dP; dM; dFm];

end

%--------------------------------------------------------------------------
end