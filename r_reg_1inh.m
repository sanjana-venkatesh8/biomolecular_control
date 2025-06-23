% Title: R-regulator 1: 1 Inhibitory connection.
% Author: Sanjana Venkatesh
% Date: June 2025
% Description: This code accompanies the LaTeX document "R-regulator 1: 1 Inhibitory connection."
%% Rescaling (S2)
global a1 a2 b1 b2 d1 d2 k1 k2 k3 k4 eta;
reset_params() % set physical parameters to the values from Alexis et al. 2023

function [A, B, C, D, E, F] = rescale_1inh()
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 k4 eta;
    A = (a1*b2) / (b1*d2);
    B = (a2*b1) / (b2*d1);
    C = k3 / eta;
    D = (k1*b1*eta) / d1^3;
    E = d2 / d1;
    F = (k2*b2*eta) / d2^3;
end

function dydt=sol_1inh(t,y)
    [A, B, C, D, E, F] = rescale_1inh();
    dydt(1,1)= 1 - y(1) + A*y(2);               %y1                                  
    dydt(2,1)= 1 - y(2) + B*y(1) - C*y(2)*y(4); %y2
    dydt(3,1)= D*y(1) - E*y(3)*y(4);            %z1
    dydt(4,1)= F*y(2) - 1/E*y(3)*y(4);          %z2
end

%% Equilibrium points (S3)
[A, B, C, D, E, F] = rescale_1inh();
y1_ss = E^2*F / (E^2*F - A*D);
y2_ss = D / (E^2*F - A*D);
y_ratio_ss = E^2*F / D;
z2_ss = 1/(C*D) * (E^2*F - A*D - D + B*E^2*F);
z1_ss = C*D^2*E*F / ((E^2*F - A*D)*(E^2*F - A*D - D + B*E^2*F));
% %% Linearization (S4)
% J = [-1 A           0           0;
%      C  -B-z2_ss    0           -y2_ss;
%      D  0           -E*z2_ss    -E*z1_ss;
%      0  D*R         -E*z2_ss    -E*z1_ss];
% fprintf("Eigenvalues of J:\n")
% display(eig(J))
% fprintf("Determinant of J: %1.4f\n", det(J))
%% Plot solution
tspace=linspace(0,100,10000); %time
[t,y]=ode23s(@sol_1inh,tspace,[0; 0; 0; 0]); % solve ODE
figure

hold on
plot(t,y(:, 1),'color','[0 1 0]', LineWidth=2, DisplayName="y_1");
plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2, DisplayName="y_2");
plot(t,y(:, 1)./y(:, 2),'--', LineWidth=2, DisplayName="y_1/y_2");
plot(t,y(:, 3),'color','[0 0 1]', LineWidth=2, DisplayName="z_1");
plot(t,y(:, 4),'color','[0 0 1]', LineWidth=2, DisplayName="z_2");

yline(y1_ss, '--g'); yline(y2_ss, '--g');
yline(y_ratio_ss, '--g');
yline(z1_ss, '--b'); yline(z2_ss, '--b');
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

xline(abs((E*F-D)*z2_ss/F), '--r') % plot the zero
xline(abs(-(1/E*z1_ss + E*z2_ss)), '--b') % plot the pole
%% Plant dynamics (S6)
A_plant = [-1    A;
           B     (-1-C*z2_ss)];
B_plant = [0; -C*y2_ss];
C_plant = [1 0;
           0 1];
D_plant = [0; 0];

plant = ss(A_plant, B_plant, C_plant, D_plant);
plant_tf = tf(plant);
bode(plant_tf) 

xline(1, '--r') % plot the zero
xline(abs(1/2 * (-(2+C*z2_ss) + sqrt((2+C*z2_ss)^2-4*(1-A*B+C*z2_ss)))), '--b') % plot pole 1
xline(abs(1/2 * (-(2+C*z2_ss) - sqrt((2+C*z2_ss)^2-4*(1-A*B+C*z2_ss)))), '--b') % plot pole 2
%% Closed-loop system
L = -ctrlr_tf*plant_tf;
s = tf('s');
L_manual = 1/((s^2+s*(2+C*z2_ss)+(1-A*B+C*z2_ss))*(s^2+s*(1/E*z1_ss+E*z2_ss))) * (A*D*C*y2_ss*z2_ss - C*F*y2_ss*(s+1)*(s+E*z2_ss));
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