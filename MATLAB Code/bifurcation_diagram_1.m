% Draw the bifurcation diagram in Chapter 2. When P* < F_m *, pair bonding
% becomes the optimal allocation of male reproductive effort. Otherwise,
% multiple mating dominates
x = [0 : 0.01: (9+sqrt(17))/8];
G= (x.^(1/2))./ (1+x.^(1/2));
F= (x.^(2) + 4*x + 1).^(1/2) - (x + 1);

% The bifurcation point
y=[(9+sqrt(17))/8 : 0.01: 5];

A = (y.^(1/2))./ (1+y.^(1/2));
B = (y.^(2) + 4*y + 1).^(1/2) - (y + 1);
% Plot the trajectory, with dotted lines if its an unstable branch
plot(x,G,'b',x,F,'r--',y,A,'b--',y,B,'r')
title('Bifurcation Diagram')
legend('P*','F_m *')
xlabel('\kappa')
ylabel('Equilibrium Population in the Non-dimensionalised Model')


