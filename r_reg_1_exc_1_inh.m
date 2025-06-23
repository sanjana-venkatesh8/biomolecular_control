% Title: R-regulator with 1 inhibitory and 1 excitatory feedback connection.
% Author: Sanjana Venkatesh
% Date: June 2025
% Description:
%% Nondimensionalization (S2)
global a1 a2 b1 b2 d1 d2 k1 k2 k3 k4 eta;
reset_params() % set physical parameters to the values from Alexis et al. 2023
k2 = 0.25; %--> stable for k2 (0, 0.343) but not for other values

function [A, B, C, D, E, F, G] = set_nondim_1exc1inh()
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 k4 eta;
    A = (a1*b2) / (b1*d2);
    B = k4 / eta;
    C = (a2*b1) / (b2*d1);
    D = (k3*d2) / (b2*eta);
    E = k1*b1*eta / (d1^3);
    F = d2/d1;
    G = (k2*b2*eta) / (d2^3);
end

function dydt=sol_1exc1inh(t,y)
    [A, B, C, D, E, F, G] = set_nondim_1exc1inh();
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 k4 eta;
    p = 2;
    % dydt(1,1) = 1 - y(1) + A*y(2) - B*y(1)*y(3);
    % dydt(2,1) = 1 - y(2) + C*y(1) + D*y(4);
    % dydt(3,1) = E*y(1) - F*y(3)*y(4);
    % dydt(4,1) = G*y(2) - 1/F * y(3)*y(4);
    dydt(1,1)=b1 - d1*y(1) + a1*y(2) - k4*y(1)*y(3);%Y1                                  
    dydt(2,1)=b2 - d2*y(2) + a2*y(1) + k3*y(4); %Y2
    dydt(3,1)=k1*y(1) - eta*y(3)*y(4); %Z1
    dydt(4,1)=k2*y(2) - eta*y(3)*y(4); %Z2
end

tspace=linspace(0,100,10000); %time
[t,y]=ode23s(@sol_1exc1inh,tspace,[0; 0; 0; 0]); % solve ODE

% calculate SS values
% lambda1 = d1*k2 - a1*k1;
% lambda2 = b2*(d1*k2 - a1*k1) - b1*(d2*k1 - a2*k2);
% Y1_ss = b1*k2 / lambda1;
% Y2_ss = b1*k1 / lambda1;
% Y_ratio_ss = k2/k1;
% Z1_ss = -(b1*k1*k2*k3)/(eta*lambda2);
% Z2_ss = -lambda2 / (k3*lambda1);

% plot solution
figure

hold on
subplot(5,1,1); plot(t,y(:, 1),'color','[0 1 0]', LineWidth=2, DisplayName="Y1");
subplot(5,1,2); plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2, DisplayName="Y2");
subplot(5,1,3); plot(t,y(:, 1)./y(:, 2),'--', LineWidth=2, DisplayName="Y1/Y2");
subplot(5,1,4); plot(t,y(:, 3),'color','[0 0 1]', LineWidth=2, DisplayName="Z1");
subplot(5,1,5); plot(t,y(:, 4),'color','[0 0 1]', LineWidth=2, DisplayName="Z2");

xline(50, '--r')
% yline(Y1_ss, '--g'); yline(Y2_ss, '--g');
% yline(Y_ratio_ss, '--g');
% yline(Z1_ss, '--b'); yline(Z2_ss, '--b');
hold off


xlabel('Time (min)','FontName', 'Times New Roman','FontSize',12) 
ylabel('Concentration (nM)','FontName', 'Times New Roman','FontSize',12) 
title('Closed loop system 5: 1 excitatory and 1 inhibitory feedback connection', 'FontName', 'Times New Roman','FontSize',12)
legend()
% legend('Y1', 'Y2', 'Y1/Y2', 'dist. time')%, 'Y1*', 'Y2*', 'Y1*/Y2*', 'Z1*', 'Z2*')
%% Equilibrium points (S3)
% have not calculated
%% Linearization (S4)
%% System dynamics (S5)
%%Controller transfer function
[A, B, C, D, E, F, G] = set_nondim_1exc1inh();

A_ctrl = [-F*z2_ss   -F*z1_ss;
         -1/F*z2_ss   -1/F*z1_ss];
B_ctrl = [E  0;
          0  G];
C_ctrl = [1  0;
          0  1];
D_ctrl = [0  0;
          0  0];

ctrlr = ss(A_ctrl, B_ctrl, C_ctrl, D_ctrl);
ctrlr_tf = tf(ctrlr);
% bode(parallel(ctrlr_tf(1), ctrlr_tf(2)))
%% Plant dynamics (S6)
A_plant = [-1    A;
          C     (-B-z2_ss)];
B_plant = [0; -y2_ss];
C_plant = [1 0;
          0 1];
D_plant = [0; 0];

plant = ss(A_plant, B_plant, C_plant, D_plant);
plant_tf = tf(plant);
bode(plant_tf) 
%% Closed-loop system
X = AnalysisPoint('X');
CL = feedback(plant * X, ctrlr);
CL_tf = tf(CL);
bode(CL_tf)

L = getLoopTransfer(CL,'X');
nyquist(L)
%%
[a,b,c,d,e,f] = set_nondim_1exc1inh()

%% Phase plane portrait

y1 = linspace(-2,8,20);
y2 = linspace(-2,2,20);
z1 = linspace(-2,2,20);

% creates two matrices one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrices of the same
% size and shape, in this case 20 rows and 20 columns
[x,y, z] = meshgrid(y1,y2, z1);
size(x)
size(y)

u = zeros(size(x));
v = zeros(size(y));
w = zeros(size(z));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
for i = 1:numel(x)
    Yprime = sol_1exc1inh(t,[x(i); z(i); 0; y(i);]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
    w(i) = Yprime(3);
end

quiver3(x,y,z,u,v,w,'r'); figure(gcf)
xlabel('y_1')
ylabel('y_2')
axis tight equal;