

function sol = plot_projection(tf)


   
    % Sample a large number of points in (kappa,gamma_f, delta) space. 


    for kappa= 1: 1: 20
        
        for gamma_f = 0.01: 0.01: 0.45
            
            for delta = 0.05: 0.05 : 0.6
                disp(kappa+" "+gamma_f+" "+delta)
                gamma_m=0.15;
                
                % Compute the ultimate value of the populations after a
                % long period of time.
                P_0  = 0;
                M_0  = 0.4;
                Fg_0= 0;
                Fm_0  = 0;
                X_0 = 0;
                Y_0= 0;
            
                initialvalues = [P_0; M_0;Fg_0; Fm_0; X_0;Y_0];
                sol = ode45(@odemodel,[0 tf], initialvalues);
                [rownum,colnum]=size(sol.y);
                
                % Compute the 'equilibrium' values (values after a long
                % time)
                P_final = sol.y(1,colnum);
                M_final = sol.y(2,colnum);
                
                % Check tyoe of equilibrium. Colour the point differently
                % based on its type.
                if P_final < 0.01 && M_final < 0.01
                    plot3(kappa,gamma_f,delta,'k*')
                    hold on;
                    
                    %disp("BLACK")
                elseif P_final < 0.01
                        plot3(kappa,gamma_f,delta,'r*')
                        hold on;
                        
                        %disp("BLUE")
                elseif M_final < 0.01
                        plot3(kappa,gamma_f,delta,'b*')
                        hold on;
                        
                        %disp("RED")
                else
                    plot3(kappa,gamma_f,delta,'y*');
                    hold on;
                    
                    %disp("YELLOW")
                end
 
            end
        end
    end



% tf = final time

% Parameters
% kappa=12;
% gamma_f=2/15;
% delta=0.108*4;

% Initial conditions
% P_0  = 0;
% M_0  = 0.4;
% Fg_0= 0;
% Fm_0  = 0;
% X_0 = 0;
% Y_0= 0;
% 
% initialvalues = [P_0; M_0;Fg_0; Fm_0; X_0;Y_0];

%--------------------------------------------------------------------------

% Simulator
% sol = ode45(@odemodel,[0 tf], initialvalues);

%--------------------------------------------------------------------------
% % Plot figure
% 
% fontsize = 12;
% 
% % Plot
% figure
% set(gca, 'FontSize', fontsize)
% hold on


% [rownum,colnum]=size(sol.y);
% P_final = sol.y(1,colnum);
% M_final = sol.y(2,colnum);
% 
% if P_final < 0.01 && M_final < 0.01
%     plot3(kappa,gamma_f,delta,'k*')
%     hold on
%     disp("BLACK")
% elseif P_final < 0.01
%         plot3(kappa,gamma_f,delta,'b*')
%         hold on
%         disp("BLUE")
% elseif M_final < 0.01
%         plot3(kappa,gamma_f,delta,'r*')
%         hold on
%         disp("RED")
% else
%     plot3(kappa,gamma_f,delta,'y*');
%     hold on
%     disp("YELLOW")
% end

% disp(sol.y)
% plot(sol.x, sol.y(1,:), 'b')
% plot(sol.x, sol.y(2,:), 'r')
% plot(sol.x, sol.y(3,:), 'g', 'LineWidth', 2)
% plot(sol.x, sol.y(4,:), 'k')
% plot(sol.x, sol.y(5,:), 'k--')
% plot(sol.x, sol.y(6,:), 'b--')
% 
% xlabel('Non-dimensionalised Time (\tau)')
% ylabel('Non-dimensionalised Populations')
% title('Population Dynamics over Time')
% legend('P', 'M', 'Fg','Fm','X','Y')
% set(gca, 'FontSize', 12)

%--------------------------------------------------------------------------
% Model
    
function dydt = odemodel(t,y)
    % Variables
    P   = y(1);
    M  = y(2);
    Fg  = y(3);
    Fm = y(4);
    X= y(5);
    Y= y(6);
    
    % ODEs
        
    dP= ((kappa*(1-P-Fm-X-Fg)*(1-P-M-Y)) / (1-P-Y)) -P*(P+Fm+Fg)-(gamma_f+gamma_m)*P;
    
    dM= (1/2)*Fm -(1/2)*(P+Fm+Fg)*M-gamma_m*M;
    
    dFg= (1/2)*(P+Fm+Fg)*(P-Fg)- (1+gamma_f+delta)*Fg+gamma_m*P;
    
    dFm= ( ( kappa*(1-P-Fm-Fg-X)*M)/(1-P-Y)) - Fm - gamma_f*Fm -(1/2)*(P+Fm+Fg)*Fm-delta*Fm;
    
    dX = gamma_f*(1-X) -(1/2)*(P+Fm+Fg)*X;
    
    dY = gamma_m*(1-Y) - (1/2)*(P+Fm+Fg)*Y;
    
    dydt = [dP; dM; dFg;dFm;dX;dY];

end

%--------------------------------------------------------------------------
end