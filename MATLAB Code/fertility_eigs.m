
% Mark the points where the coexistence
% equilibrium is stable (coloured blue) and the points where it doesn't
% exist or is unstable ( coloured red). This is from Chapter 4- to generate
% Figure 4.3.

hold on


% First generate the curves from the analysis
x= 1.65:0.01:15.00;

y1= ((7-7.*x-8.*x.^2)./(2*(8.*x-3))) + (1/2)*(( (25+20.*x-231.*x.^2+176.*x.^3+64.*x.^4)./((8.*x-3).^2)).^(1/2));

plot(x,y1,"g");

x2=0.75: 1:15.00;

disp(length(x2));

y2= x2;

for i=1 : length(y2)
   kappa=x2(i);
   r=roots([6,20-2*kappa,18+15*kappa-8*kappa^2,4-7*kappa+2*kappa^2]);
   if imag(r(1))==0 && real(r(1)) > 0.02 && real(r(1)) < 1
       y2(i)=r(1);
   elseif  imag(r(2))==0 && real(r(2)) > 0.02 && real(r(2)) < 1
           y2(i)=r(2);
   else
       y2(i)=r(3);
   end
end

disp(x2)
disp(y2)

plot(x2,y2,"c")

% Plot the upper and lower bounds on coexistence
x3 = 0:0.01:15.00;
y3 = (1/2).* (sqrt(x3.*(x3+2)) - x3);

plot(x3,y3,"b");
           
   
% Check if each point is stable, and corroborate the evidence using
% the plot. We should see that the curves separate the space into clear
% boundaries.


% all_k=0.1: 0.5: 15.00;
% 
% all_gamma = 0.01: 0.01: 0.3;
% 
% for a=1: length(all_k)
%     
%     for b=1: length(all_gamma)
%         
%         k= all_k(a);
%         gamma=all_gamma(b);
%         
%         if gamma > (1/2)*((k*(k+2))^(1/2) -k)
%             
%             plot(k,gamma,"r*");
%             disp("RED");
%             continue
%         end
%         
%       
%         
%         P = (-16*k*gamma^2 + k*(8*gamma+1)*(4*gamma^2+12*gamma+17)^(1/2)-14*k*gamma-9*k+12*gamma^2+28*gamma+8)/(2*(6*gamma^2+14*gamma+4));
%         
%         F= (1+k-P^2+2*gamma- 2*P*gamma- (k^2-4*k*(P-1)+(P-1)^2)^(1/2))/(P-1);
%         
%         M= (F/(P+F));
%         
%         X= (2*gamma)/ (2*gamma+P+F);
%         
%         if P >=0 && P < 1 && F >=0 && F < 1 && M > 0 && M <=1 && X >=0 && X < 1
%             
%              disp(k);
%              disp(gamma);
%         
%         
%             jacobian = zeros(3,3);
% 
%             jacobian(1,1)= (F*(k*M-(P-1)^2)-(P-1)^2*(k+gamma+2*P)+k*M*X)/(P-1)^2;
% 
%             jacobian(1,2)= (-k*(1-P-F-X))/(1-P);
% 
%             jacobian(1,3)= ((-k*(1-P-M))/(1-P))-P;
% 
%             jacobian(1,4)=((-k*(1-P-M))/(1-P));
% 
%             jacobian(2,1)= -M/2;
% 
%             jacobian(2,2)= (1/2)*(-F-P);
% 
%             jacobian(2,3)= (1-M)/2;
% 
%             jacobian(2,4)= 0;
% 
%             jacobian(3,1)= -(F/2)- ((k*M)/(1-P)) + ( (k*M*(1-P-F-X))/ (1-P)^2);
% 
%             jacobian(3,2)= k*(1-P-F-X)/ (1-P); 
% 
%             jacobian(3,3)= -1- F/2 - gamma - ((k*M)/(1-P)) + (1/2)*(-F-P);
% 
%             jacobian(3,4)= (-k*M)/(1-P);
% 
%             jacobian(4,1)= -X/2;
% 
%             jacobian(4,2)= 0;
% 
%             jacobian(4,3)= -X/2;
% 
%             jacobian(4,4)= -gamma + (1/2)*(-F-P);
% 
%             disp(jacobian)
% 
%             try
%                 e= eig(jacobian);
%                 disp(e);
% 
%                 if real(e(1,1)) < 0 && real(e(2,1)) < 0 && real(e(3,1)) < 0 && real(e(4,1)) < 0
% 
%                     disp("BLUE")
%                     plot(k,gamma,"b*")
% 
% 
%                 else
% 
%                     disp("RED")
%                     plot(k,gamma,"r*")
% 
%                 end
% 
%             catch exception
%                 continue
%             end
%       
%         else
%             disp(k)
%             disp(gamma)
%             
%             disp("RED")
%             plot(k,gamma,"r*")
%         end
%         
%         
%     end
%     
% end
% 

xlabel('\kappa = mating rate/ maturation rate')
ylabel('\gamma = rate of female fertility loss/ maturation rate')
title('Region of Stability for each strategy')
legend('$\alpha(\kappa)$', '$\beta(\kappa)$', '$(1/2) ( \sqrt{\kappa(\kappa+2)} - \kappa)$','Interpreter','latex');
        
        
        


