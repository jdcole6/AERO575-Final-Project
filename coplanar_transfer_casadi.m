function main
% coplanar transfer problem solved using CASADI
clear all; clc; close all

import casadi.*

Re      = 6378;      % Earth Radius (km)
mu      = 398600.4;  % Gravitational constant for Earth (km^3/s^2) 

% Declare model variables
nx = 4; nu=2;
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
x = [x1; x2; x3; x4];
u1 = SX.sym('u1');
u2 = SX.sym('u2');
u = [u1; u2];

% Model equations
xdot = optimalTransfer(x, u);

% Continuous time dynamics
f = Function('f', {x, u}, {xdot});

% Initial conditions
x10 = Re+600;
x20 = 0;
x30 = 0;
x40 = sqrt(mu/(Re+600)^3);

% Terminal conditions
x1f = Re+2000;
x2f = 0;
x4f = sqrt(mu/(Re+2000)^3);

% Constraints
tf = 7200; % final time
dt = 2;
% Start with an empty NLP
w   =  {};
w0  = [];
lbw = [];
ubw = [];
J   = 0;
g   = {};
lbg = [];
ubg = [];

% Set initial conditions
Xk  = MX.sym('X0', nx);
w   = {w{:}, Xk};
lbw = [lbw; x10; x20; x30; x40];
ubw = [ubw; x10; x20; x30; x40];
w0  = [w0;  x10; x20; x30; x40];

N  = tf/dt; % number of steps discretizing the problem

% Formulate the NLP
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],nu);
    w = {w{:}, Uk};
    lbw = [lbw; -Inf; -Inf;];
    ubw = [ubw;  Inf; Inf];
    w0  = [w0;  0; 0];
    
    % New NLP variable for state at end of interval
    Xkp1 = MX.sym(['X_' num2str(k+1)], nx);
    w = [w, {Xkp1}];
    lbw = [lbw; -Inf;    -Inf;  -Inf; -Inf];
    ubw = [ubw;  Inf;     Inf;   Inf;  Inf];
    w0 = [w0; x10; x20; x30; x40];

    % Add equality constraint
    g = [g, {dt*f(Xk, Uk)+Xk-Xkp1}];

    lbg = [lbg; 0; 0; 0; 0];
    ubg = [ubg; 0; 0; 0; 0];
    
    % update X_{k+1}
    Xk = Xkp1;
    
    % Cost function 
     J = J+dt*(Uk'*Uk);
    
end

% Set terminal conditions
lbw(end-3) = x1f;
ubw(end-3) = x1f;
w0(end-3)  = x1f;

lbw(end-2) = x2f;
ubw(end-2) = x2f;
w0(end-2) = x2f;

%lbw(end-1) = x3f;
%ubw(end-1) = x3f;
%w0(end-1) = x3f;

lbw(end) = x4f;
ubw(end) = x4f;
w0(end)  = x4f;

% Supply gradient (not used in this example)
dJdw = gradient(J,vertcat(w{:}));
dGdw = jacobian(vertcat(g{:}),vertcat(w{:}));


% Create an NLP solver
prob = struct('f',  J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, ...
            'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);


% Plot the solution
x1_opt = w_opt(1:6:end);
x2_opt = w_opt(2:6:end);
x3_opt = w_opt(3:6:end);
x4_opt = w_opt(4:6:end);
u1_opt  = w_opt(5:6:end);
u2_opt = w_opt(6:6:end);
tgrid = linspace(0, tf, N+1);

ocpSol.t = tgrid';
ocpSol.X = [x1_opt x2_opt x3_opt x4_opt];
ocpSol.U = [u1_opt,u2_opt; u1_opt(end),u2_opt(end)];



figure(2);

subplot(411); hold on
plot(tgrid, ocpSol.X(:,1),'linewidth',3);set(gca,'fontsize',16)
xlabel('Time')
ylabel('x_1')

subplot(412); hold on
plot(tgrid, ocpSol.X(:,2),'linewidth',3);set(gca,'fontsize',16)
xlabel('Time')
ylabel('x_2')

subplot(413); hold on
plot(tgrid, ocpSol.X(:,3),'linewidth',3);set(gca,'fontsize',16)
xlabel('Time')
ylabel('x_3')

subplot(414); hold on
stairs(tgrid, ocpSol.X(:,4),'linewidth',3);set(gca,'fontsize',16)
xlabel('Time')
ylabel('x_4')
set(gcf,'Position',[0 0 600 700]);

figure(3);
subplot(211)
plot(tgrid, ocpSol.U(:,1),'linewidth',3);set(gca,'fontsize',16)
xlabel('Time')
ylabel('u1')
subplot(212)
plot(tgrid, ocpSol.U(:,2),'linewidth',3);set(gca,'fontsize',16)
xlabel('Time')
ylabel('u2')

% Plot the results
Xp = ocpSol.X;
xT      = [x1f;x2f;NaN;x4f];   % target state
x0      = [x10;x20;x30;x40];     % initial state
figure;plot(Xp(:,1).*cos(Xp(:,3)),Xp(:,1).*sin(Xp(:,3)),'linewidth',2);
t=0:0.1:2*pi;hold on; plot(xT(1)*cos(t),xT(1)*sin(t),'r--');
hold on; plot(x0(1)*cos(t),x0(1)*sin(t),'g--');
xlabel('x,km');ylabel('y,km')
axis('square');set(gca,'fontsize',16)

return;

function [dxdt] = optimalTransfer(x, u)
Re      = 6378;      % Earth Radius (km)
mu      = 398600.4;  % Gravitational constant for Earth (km^3/s^2) 

ar   = u(1); % radial accel
ath  = u(2); % tangential accel

x1=x(1);x2=x(2);x3=x(3);x4=x(4);

% EOMs
xdot1 = x2;
xdot2 = x1*x4^2 - mu/(x1^2) + ar;
xdot3 = x4;
xdot4 = 1/x1*(-2*x2*x4 + ath);
 

dxdt      = [xdot1;xdot2;xdot3;xdot4];

return;