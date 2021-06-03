% Plots the eigenvalues in the region of coexistence (Chapter 3) to show that these
% coexistence equilibria are stable (i.e every eigenvalue has a negative
% real part). We do this by sampling every kappa in the biologically
% relevant domain and plotting the eigenvalues as functions of kappa 
x = 0.719:0.001:1.64;

e1=zeros([1 length(x)]);

e2=zeros([1 length(x)]);

e3= zeros([1 length(x)]);


% Go through each point
for n=1 : length(x)
    
 k = x(n);
 % Compute populations
 P = 1+ k*((sqrt(17)-9)/8);

 Fm = (-k*P +k -(P)^2 + P - sqrt(k-k*P))/(k+P-1);
 
 
  M = Fm / (P + Fm);
 
 
 % Compute Jacobian matrix
 jacobian = zeros(3,3);
 
 jacobian(1,1)= k*( (M*Fm / (P-1)^2) -1) -2*P- Fm;
 
 jacobian(1,2) = -k*(Fm+P-1)/(P-1);
 
 jacobian(1,3) = (-k*(1-P-M)/(1-P));
 
 jacobian(2,1) = -M/2;
 
 jacobian(2,2) = (1/2)*(-Fm-P);
 
 jacobian(2,3) = (1/2)*(1-M);
 
 jacobian(3,1) = Fm*((k*M/(1-P)^2)-(1/2));
 
 jacobian(3,2) = k*(Fm+P-1)/ (P-1);
 
 jacobian(3,3) = -Fm-(P/2) - 1 + (k*M/(P-1));
 
 % Compute real parts of eigenvalues
 
 e=eig(jacobian);
 
 e1(1,n) = real(e(1,1));
 
 e2(1,n) = real(e(2,1));
 
 e3(1,n) = real(e(3,1));
 
 % When are all three eigenvalues negative?
 
 if e1(1,n) < 0 && e2(1,n) < 0 && e3(1,n) < 0
     
    disp(k)
 end
 
 
end

figure

hold on;

% Sort eigenvalues to make sure we're connecting the right curves.

A = [e1' e2' e3'];
sortedA = sort(A,2)
plot(x,sortedA(:,1), 'r')
plot(x,sortedA(:,2), 'g')
plot(x,sortedA(:,3), 'b')

title('Eigenvalue Plot for the Coexistence Equilibrium Points')

xlabel('\kappa')
ylabel('Real parts of eigenvalues Re(\lambda (\kappa))')



% Plot each of the eigenvalues (lambda_1, lambda_2, lambda_3) as functions
% of kappa

%plot(x,e1)

%plot(x,e2)

%plot(x,e3)
