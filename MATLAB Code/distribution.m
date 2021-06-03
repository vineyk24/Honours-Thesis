
% Plot the distribution of OSR along the surface alpha(kappa, delta) in
% Chapter 7. 
function distribution()

    all_OSR=zeros(1,10);
    % Get 1000 values of OSR
    for i=1: 1: 1000
        % Sample points
        x=get_x(1000);
        
        % Compute OSR
        OSR = 1/(1-x);
        
        all_OSR(i)=OSR;
    end
    % Make a histogram of it
    all_OSR
    histogram(all_OSR)
    title("Distribution of the OSR on the surface \alpha(\kappa,\delta)")
    xlabel("OSR")
    ylabel("Frequency")
        
    
    




function x = get_x(tf)


% tf = final time

% Parameters sampled at random from particular ranges.
kappa= (20)*rand;
delta = (0.6)*rand;
gamma = ((-3-3*delta) + (3+2*delta)*kappa+(2*kappa*kappa)+(1-2*kappa)*sqrt((delta+1)^2+4*kappa+ 2*delta*kappa+kappa^2)) / (2*(2+delta-kappa+ sqrt((delta+1)^2+4*kappa+ 2*delta*kappa+kappa^2)));

% Initial conditions
P_0  = 0;
M_0  = 0.4;
Fm_0  = 0;
X_0 = 0;

initialvalues = [P_0; M_0; Fm_0; X_0];

%--------------------------------------------------------------------------

% Simulator
sol = ode45(@odemodel,[0 tf], initialvalues);

x = sol.y(4,end);

%--------------------------------------------------------------------------
% Plot figure

% fontsize = 12;

% % Plot
% figure
% set(gca, 'FontSize', fontsize)
% hold on
% 
% plot(sol.x, sol.y(1,:), 'b')
% plot(sol.x, sol.y(2,:), 'r')
% plot(sol.x, sol.y(3,:), 'g', 'LineWidth', 2)
% plot(sol.x, sol.y(4,:), 'k')
% %plot(sol.x, sol.y(5,:), 'k--')
% 
% xlabel('Non-dimensionalised Time (\tau)')
% ylabel('Non-dimensionalised Populations')
% title('Population Dynamics over Time')
% legend('P', 'M', 'Fm','X')
% set(gca, 'FontSize', 12)

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

end