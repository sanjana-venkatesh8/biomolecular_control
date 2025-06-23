% Title: R-regulator 1: 1 excitatory connection.
% Author: Sanjana Venkatesh
% Date: June 2025
% Description: This code accompanies the LaTeX document "R-regulator 1: 1 Inhibitory connection."
%% Rescaling (S2)
global a1 a2 b1 b2 d1 d2 k1 k2 k3 k4 eta;
reset_params() % set physical parameters to the values from Alexis et al. 2023
k2 = 0.5*(b1*d2+a1*b2)/(b2*d1+a2*b1);

function [A, B, C, D, E, F] = rescale_1exc()
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 k4 eta;
    A = (a1*b2) / (b1*d2);
    B = (a2*b1) / (b2*d1); 
    C = (k3*d2) / (b2*eta);
    D = (k1*b1*eta) / d1^3;
    E = d2 / d1;
    F = (k2*b2*eta) / d2^3;
end

function dydt=sol_1exc(t,y)
    [A, B, C, D, E, F] = rescale_1exc();
    dydt(1,1)= 1 - y(1) + A*y(2);               %y1                                  
    dydt(2,1)= 1 - y(2) + B*y(1) + C*y(4); %y2
    dydt(3,1)= D*y(1) - E*y(3)*y(4);            %z1
    dydt(4,1)= F*y(2) - 1/E*y(3)*y(4);          %z2
end
%% Equilibrium points (S3)
[A, B, C, D, E, F] = rescale_1exc();
y1_ss = E^2*F / (E^2*F - A*D);
y2_ss = D / (E^2*F - A*D);
z2_ss = 1/(C*(E^2*F - A*D))*(-E^2*F + A*D + D - B*E^2*F);
z1_ss = (C*D*E*F) / (-E^2*F + A*D + D - B*E^2*F);
y1_ss, y2_ss, z1_ss, z2_ss
%% Linearization (S4)
J = [-1 A     0           0;
     B  -1    0           C;
     D  0     -E*z2_ss    -E*z1_ss;
     0  F     -1/E*z2_ss  -1/E*z1_ss];
fprintf("Eigenvalues of J:\n")
display(eig(J))
fprintf("Determinant of J: %1.4f\n", det(J))
%% Plot solution
tspace=linspace(0,100,10000); %time
[t,y]=ode23s(@sol_1exc,tspace,[0; 0; 0; 0]); % solve ODE

figure

hold on
subplot(5,1,1); plot(t,y(:, 1),'color','[0 1 0]', LineWidth=2, DisplayName="Y1");
subplot(5,1,2); plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2, DisplayName="Y2");
subplot(5,1,3); plot(t,y(:, 1)./y(:, 2),'--', LineWidth=2, DisplayName="Y1/Y2");
subplot(5,1,4); plot(t,y(:, 3),'color','[0 0 1]', LineWidth=2, DisplayName="Z1");
subplot(5,1,5); plot(t,y(:, 4),'color','[0 0 1]', LineWidth=2, DisplayName="Z2");