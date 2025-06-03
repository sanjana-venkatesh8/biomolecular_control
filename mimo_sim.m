%% Simulating alternative feedback topologies for the R-regulator system

%% CONSTANTS
global a1 a2 b1 b2 d1 d2 k1 k2 k3 eta;
a1 = 0.1; % min^-1
a2 = 0.4; % min^-1
b1 = 2;%e-9; % M * min^-1
b2 = 1;%e-9; % M * min^-1
d1 = 1; % min^-1
d2 = 1; % min^-1
k1 = 0.5; % min^-1
k2 = 1; % min^-1
k3 = 2;%e-9; % M * min^-1
eta = 10;%e-9; % M * min^-1
%% 0. Open-loop system
% 

%% OL system (Figure 3B)
%time
tspace=linspace(0,100,10000);

[t,y]=ode23s(@sol_OL,tspace,[0;0]);
    
function dydt = sol_OL(t, y)
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 eta;
    p = 2;
    dydt(1, 1) = b1 - d1*y(1) + a1*y(2) + p*heaviside(t-50);%Y1
    dydt(2, 1)= b2 - d2*y(2) + a2*y(1); %Y2
end

% calculate steady-state values
Y1_ss = (a2*b1 + b2*d1) / (d1*d2 - a1*a2);
Y2_ss = (a1*b2 + b1*d2) / (d1*d2 - a1*a2);

figure

hold on
plot(t,y(:, 1),'color','[0 0 0]', LineWidth=2);
plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2);
xline(50, '--r') % show time of disturbance
yline(Y1_ss, '--b'); yline(Y2_ss, '--b') % plot SS values
hold off

xlabel('Time (min)','FontName', 'Times New Roman','FontSize',12) 
ylabel('Concentration (nM)','FontName', 'Times New Roman','FontSize',12) 
title('Open loop response (Figure 3B)', 'FontName', 'Times New Roman','FontSize',12)
legend('Y1', 'Y2', 'dist. time', 'Y1*', 'Y2*')
%% 1. 1 inhibitory feedback connection
% 
% 
% 

%% CL system 1: 1 inhibitory fb connection (R-regulator) (Figure 4A)
% ODE function
function dydt=sol_1inh(t,y)
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 eta;
    p = 2;
    dydt(1,1)=b1 - d1*y(1) + a1*y(2)+ p*heaviside(t-50);%Y1                                  
    dydt(2,1)=b2 - d2*y(2) + a2*y(1) - k3*y(2)*y(4); %Y2
    dydt(3,1)=k1*y(1) - eta*y(3)*y(4); %Z1
    dydt(4,1)=k2*y(2) - eta*y(3)*y(4); %Z2
end

tspace=linspace(0,100,10000); %time
[t,y]=ode23s(@sol_1inh,tspace,[0; 0; 0; 0]); % solve ODE

% calculate SS values
lambda1 = d1*k2 - a1*k1;
lambda2 = b2*(d1*k2 - a1*k1) - b1*(d2*k1 - a2*k2);
Y1_ss = b1*k2 / lambda1;
Y2_ss = b1*k1 / lambda1;
Y_ratio_ss = k2/k1;
Z1_ss = ((b1*k1)^2*k2*k3) / (eta*lambda1*lambda2);
Z2_ss = lambda2 / (k3*b1*k1);

figure

hold on
plot(t,y(:, 1),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 1)./y(:, 2),'--', LineWidth=2);
plot(t,y(:, 3),'color','[0 0 1]', LineWidth=2);
plot(t,y(:, 4),'color','[0 0 1]', LineWidth=2);

xline(50, '--r')
yline(Y1_ss, '--g'); yline(Y2_ss, '--g');
yline(Y_ratio_ss, '--g');
yline(Z1_ss, '--b'); yline(Z2_ss, '--b');
hold off

xlabel('Time (min)','FontName', 'Times New Roman','FontSize',12) 
ylabel('Concentration (nM)','FontName', 'Times New Roman','FontSize',12) 
title('Closed loop system 1: 1 inhibitory feedback connection (R-regulator) (Figure 4A)', 'FontName', 'Times New Roman','FontSize',12)
legend('Y1', 'Y2', 'Y1/Y2', 'dist. time', 'Y1*', 'Y2*', 'Y1*/Y2*', 'Z1*', 'Z2*')
%% 2. 1 excitatory feedback connection
% 
% Question: Why do we apply the disturbance specifically to the Y1_dot term?

%% CL system 2: 1 excitatory fb connection
% ODE function
function dydt=sol_1exc(t,y)
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 eta;
    p = 2;
    dydt(1,1)=b1 - d1*y(1) + a1*y(2) + p*heaviside(t-50);%Y1                                  
    dydt(2,1)=b2 - d2*y(2) + a2*y(1) + k3*y(4); %Y2
    dydt(3,1)=k1*y(1) - eta*y(3)*y(4); %Z1
    dydt(4,1)=k2*y(2) - eta*y(3)*y(4); %Z2
end

tspace=linspace(0,100,10000); %time
[t,y]=ode23s(@sol_1exc,tspace,[0; 0; 0; 0]); % solve ODE

% calculate SS values
lambda1 = d1*k2 - a1*k1;
lambda2 = b2*(d1*k2 - a1*k1) - b1*(d2*k1 - a2*k2);
Y1_ss = b1*k2 / lambda1;
Y2_ss = b1*k1 / lambda1;
Y_ratio_ss = k2/k1;
Z1_ss = -(b1*k1*k2*k3)/(eta*lambda2);
Z2_ss = -lambda2 / (k3*lambda1);

% plot solution
figure

hold on
plot(t,y(:, 1),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 1)./y(:, 2),'--', LineWidth=2);
plot(t,y(:, 3),'color','[0 0 1]', LineWidth=2);
plot(t,y(:, 4),'color','[0 0 1]', LineWidth=2);

xline(50, '--r')
yline(Y1_ss, '--g'); yline(Y2_ss, '--g');
yline(Y_ratio_ss, '--g');
yline(Z1_ss, '--b'); yline(Z2_ss, '--b');
hold off


xlabel('Time (min)','FontName', 'Times New Roman','FontSize',12) 
ylabel('Concentration (nM)','FontName', 'Times New Roman','FontSize',12) 
title('Closed loop system 2: 1 excitatory feedback connection', 'FontName', 'Times New Roman','FontSize',12)
legend('Y1', 'Y2', 'Y1/Y2', 'dist. time', 'Y1*', 'Y2*', 'Y1*/Y2*', 'Z1*', 'Z2*')
%%
J = [-d1  a1          0          0;
      a2 -d2          0         k3;
      k1   0 -eta*Z2_ss -eta*Z1_ss;
       0  k2 -eta*Z2_ss -eta*Z1_ss];
J, eig(J)
%% 
% 
% 
% This system is unstable.
%% 3. 2 inhibitory feedback connections
% 

%% CL system 3: 2 inhibitory fb connections
% ODE function
function dydt=sol_2inh(t,y)
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 eta;
    p = 2; k4 = k3;
    dydt(1,1)=b1 - d1*y(1) + a1*y(2) - k4*y(3)*y(1) + p*heaviside(t-50);%Y1                                  
    dydt(2,1)=b2 - d2*y(2) + a2*y(1) - k3*y(4)*y(2); %Y2
    dydt(3,1)=k1*y(1) - eta*y(3)*y(4); %Z1
    dydt(4,1)=k2*y(2) - eta*y(3)*y(4); %Z2
end

tspace=linspace(0,100,10000); %time
[t,y]=ode23s(@sol_2inh,tspace,[0; 0; 0; 0]); % solve ODE

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
plot(t,y(:, 1),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 1)./y(:, 2),'--', LineWidth=2);
plot(t,y(:, 3),'color','[0 0 1]', LineWidth=2);
plot(t,y(:, 4),'color','[0 0 1]', LineWidth=2);

xline(50, '--r')
% yline(Y1_ss, '--g'); yline(Y2_ss, '--g');
% yline(Y_ratio_ss, '--g');
% yline(Z1_ss, '--b'); yline(Z2_ss, '--b');
hold off


xlabel('Time (min)','FontName', 'Times New Roman','FontSize',12) 
ylabel('Concentration (nM)','FontName', 'Times New Roman','FontSize',12) 
title('Closed loop system 3: 2 inhibitory feedback connections', 'FontName', 'Times New Roman','FontSize',12)
legend('Y1', 'Y2', 'Y1/Y2', 'dist. time');%, 'Y1*', 'Y2*', 'Y1*/Y2*', 'Z1*', 'Z2*')
%% 4. 2 excitatory feedback connections
% 

%% CL system 4: 2 excitatory fb connections
% ODE function
function dydt=sol_2exc(t,y)
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 eta;
    p = 2; k4 = k3;
    dydt(1,1)=b1 - d1*y(1) + a1*y(2) + k4*y(3) + p*heaviside(t-50);%Y1                                  
    dydt(2,1)=b2 - d2*y(2) + a2*y(1) + k3*y(4); %Y2
    dydt(3,1)=k1*y(1) - eta*y(3)*y(4); %Z1
    dydt(4,1)=k2*y(2) - eta*y(3)*y(4); %Z2
end

tspace=linspace(0,100,10000); %time
[t,y]=ode23s(@sol_2exc,tspace,[0; 0; 0; 0]); % solve ODE

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
plot(t,y(:, 1),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 1)./y(:, 2),'--', LineWidth=2);
plot(t,y(:, 3),'color','[0 0 1]', LineWidth=2);
plot(t,y(:, 4),'color','[0 0 1]', LineWidth=2);

xline(50, '--r')
% yline(Y1_ss, '--g'); yline(Y2_ss, '--g');
% yline(Y_ratio_ss, '--g');
% yline(Z1_ss, '--b'); yline(Z2_ss, '--b');
hold off


xlabel('Time (min)','FontName', 'Times New Roman','FontSize',12) 
ylabel('Concentration (nM)','FontName', 'Times New Roman','FontSize',12) 
title('Closed loop system 4: 2 excitatory feedback connections', 'FontName', 'Times New Roman','FontSize',12)
legend('Y1', 'Y2', 'Y1/Y2', 'dist. time');%, 'Y1*', 'Y2*', 'Y1*/Y2*', 'Z1*', 'Z2*')
%% 5. 1 excitatory, 1 inhibitory feedback connection
% 

%% CL system 5: 1 inhibitory, 1 excitatory fb connection
% ODE function
function dydt=sol_1exc1inh(t,y)
    global a1 a2 b1 b2 d1 d2 k1 k2 k3 eta;
    p = 2; k4 = k3;
    dydt(1,1)=b1 - d1*y(1) + a1*y(2) - k4*y(3)*y(1) + p*heaviside(t-50);%Y1                                  
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
plot(t,y(:, 1),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 2),'color','[0 1 0]', LineWidth=2);
plot(t,y(:, 1)./y(:, 2),'--', LineWidth=2);
plot(t,y(:, 3),'color','[0 0 1]', LineWidth=2);
plot(t,y(:, 4),'color','[0 0 1]', LineWidth=2);

xline(50, '--r')
% yline(Y1_ss, '--g'); yline(Y2_ss, '--g');
% yline(Y_ratio_ss, '--g');
% yline(Z1_ss, '--b'); yline(Z2_ss, '--b');
hold off


xlabel('Time (min)','FontName', 'Times New Roman','FontSize',12) 
ylabel('Concentration (nM)','FontName', 'Times New Roman','FontSize',12) 
title('Closed loop system 5: 1 excitatory and 1 inhibitory feedback connection', 'FontName', 'Times New Roman','FontSize',12)
legend('Y1', 'Y2', 'Y1/Y2', 'dist. time')%, 'Y1*', 'Y2*', 'Y1*/Y2*', 'Z1*', 'Z2*')