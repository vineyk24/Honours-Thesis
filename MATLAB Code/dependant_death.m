
% The numerical solution to the model with dependant death in Chapter 5. 

function sol = dependant_death(tf)



% tf = final time

% Parameters

% Hunter-gatherer like parameter estimates
kappa=12;
gamma=2/15;
delta=0.432;

% Initial conditions
P_0  = 0;
M_0  = 0.4;
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
        
    dP= ((kappa*(1-P-Fm-X)*(1-P-M)) / (1-P)) -P*(P+Fm)-gamma*P;
    
    dM= (1/2)*Fm -(1/2)*(P+Fm)*M;
    
    dFm= ( ( kappa*(1-P-Fm-X)*M)/(1-P)) - Fm - gamma*Fm -(1/2)*(P+Fm)*Fm-delta*Fm;
    
    dX = gamma*(1-X) -(1/2)*(P+Fm)*X;
    
    dydt = [dP; dM; dFm;dX];

end

%--------------------------------------------------------------------------
end