% Title: R-regulator 1: 1 Inhibitory connection from Z2-->Y1.
% Author: Sanjana Venkatesh
% Date: June 2025
% Description: This code accompanies the LaTeX document "R-regulator 1: 1 Inhibitory connection."
%% Rescaling (S2)
global a1 a2 b1 b2 d1 d2 k1 k2 k3 k4 eta;
reset_params() % set physical parameters to the values from Alexis et al. 2023
k2 = k1 * d2/a2;

function [A, B, D, E, F, G] = rescale_1inhCross()
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 k4 eta;
    A = (a1*b2) / (b1*d2);
    B = (a2*b1) / (b2*d1);
    D = (k1*b1*eta) / d1^3;
    E = d2 / d1;
    F = (k2*b2*eta) / d2^3;
    G = k3 / eta;
end

function dydt=sol_1inhCross(t,y)
    [A, B, D, E, F, G] = rescale_1inhCross();
    dydt(1,1)= 1 - y(1) + A*y(2) - E*G*y(1)*y(4);%y1                                  
    dydt(2,1)= 1 - y(2) + B*y(1); %y2
    dydt(3,1)= D*y(1) - E*y(3)*y(4);            %z1
    dydt(4,1)= F*y(2) - 1/E*y(3)*y(4);          %z2
end

%% Equilibrium points (S3)
[A, B, D, E, F, G] = rescale_1inhCross();
y1_ss = E^2*F / (D - B*E^2*F);
y2_ss = D / (D - B*E^2*F);
y_ratio_ss = E^2*F / D;
z1_ss = D*E^4*F^2*G / ((D - B*E^2*F)*(D*(1+A) - E^2*F*(1+B)));
z2_ss = (D*(1+A) - E^2*F*(1+B)) / (E^3*F*G);
%% Linearization (S4)
J = [(-1-E*G*z2_ss) A           0           -E*G*y1_ss;
     B              -1          0           0;
     D              0           -E*z2_ss    -E*z1_ss;
     0              F           -1/E*z2_ss  -1/E*z1_ss];
fprintf("Eigenvalues of J:\n")
display(eig(J))
fprintf("Determinant of J: %1.4f\n", det(J))
%% Plot solution
tspace=linspace(0,100,10000); %time
[t,y]=ode23s(@sol_1inhCross,tspace,[0; 0; 0; 0]); % solve ODE
figure

hold on; legend
plot(t,y(:, 1),'color','[0 1 0]', LineWidth=2, DisplayName="y_1");
plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2, DisplayName="y_2");
plot(t,y(:, 1)./y(:, 2),'--', LineWidth=2, DisplayName="y_1/y_2");
plot(t,y(:, 3),'color','[0 0 1]', LineWidth=2, DisplayName="z_1");
plot(t,y(:, 4),'color','[0 0 1]', LineWidth=2, DisplayName="z_2");

yline(y1_ss, '--g', DisplayName="y_1 (ss)"); yline(y2_ss, '--g', DisplayName="y_2 (ss)");
yline(y_ratio_ss, '--', DisplayName="y_1/y_2 (ss)");
yline(z1_ss, '--b', DisplayName="z_1 (ss)"); yline(z2_ss, '--b', DisplayName="z_2 (ss)");
hold off
%% Controller dynamics (S5)
A_ctrl = [-E*z2_ss   -E*z1_ss;
         -1/E*z2_ss   -1/E*z1_ss];
B_ctrl = [D  0;
          0  F];
C_ctrl = [0  1];
D_ctrl = [0  0];

ctrlr = ss(A_ctrl, B_ctrl, C_ctrl, D_ctrl);
ctrlr_tf = tf(ctrlr);
bode(parallel(ctrlr_tf(1), ctrlr_tf(2)))

xline(abs((E*F-D)*z2_ss/F), '--r', DisplayName="Zero") % plot the zero
xline(abs(-(1/E*z1_ss + E*z2_ss)), '--b', DisplayName="Pole") % plot the pole
legend
%% Plant dynamics (S6)
A_plant = [(-1-E*G*z2_ss)    A;
           B                -1];
B_plant = [-E*G*y1_ss; 0];
C_plant = [1 0;
           0 1];
D_plant = [0; 0];

plant = ss(A_plant, B_plant, C_plant, D_plant);
plant_tf = tf(plant);
bode(plant_tf) 

% xline(1, '--r') % plot the zero
xline(abs(1/2 * (-(2+E*G*z2_ss) + sqrt(4*A*B+(E*G*z2_ss)^2))), '--b') % plot pole 1
xline(abs(1/2 * (-(2+E*G*z2_ss) - sqrt(4*A*B+(E*G*z2_ss)^2))), '--b') % plot pole 2
%% Closed-loop system
L = ctrlr_tf*plant_tf;
% s = tf('s');
% L_manual = 1/((s^2+s*(2+C*z2_ss)+(1-A*B+C*z2_ss))*(s^2+s*(1/E*z1_ss+E*z2_ss))) * (A*D*C*y2_ss*z2_ss - C*F*y2_ss*(s+1)*(s+E*z2_ss));
bode(1/(1+L))
% figure; bode(1/(1+L_manual))
% X = AnalysisPoint('X');
% CL = feedback(plant * X, ctrlr);
% CL_tf = tf(CL);
% bode(CL_tf)
% 
% L = getLoopTransfer(CL,'X');
figure; nyquist(L)
% %%
% [S, ~, L, ~] = loopsens(plant_tf, ctrlr_tf);
% figure; bode(S)
% figure; nyquist(L)