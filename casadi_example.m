addpath('C:\Users\Owner\Documents\MATLAB\AERO 575\Final Project\Casadi Files')
import casadi.*

clear;clc;close all

%%
% 0 for free, 1 for fixed
fixed_tf = 0;

%%

% Discretization variables
T = MX.sym('T');
N = 100;
dt = T/N;

% Model variables 
nx = 3;
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x = [x1; x2; x3];
u = SX.sym('u');

% Model equations
xdot = rocket_model(x, u);

% Continuous time 
f = Function('f', {x, u}, {xdot});

% Initial conditions
v0 = 0;
h0 = 0;
m0 = 214.839;

% Terminal conditions
mf = 67.9833;

% Constraints
Fm = 9.525515;
tf = 100;

% Empty NLP
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g = {};
lbg = [];
ubg = [];

% Final time condition
w = {w{:}, T};

if fixed_tf == 0
    lbw = [lbw; 0];
    ubw = [ubw; inf];
else
    lbw = [lbw; tf];
    ubw = [ubw; tf];
end
w0 = [w0; 0];

% Set initial conditions
Xk = MX.sym('X0', nx);
w = {w{:}, Xk};
lbw = [lbw; v0; h0; m0];
ubw = [ubw; v0; h0; m0];
w0 = [w0; 0; 0; m0];

% Fomrulate the NLP
for k = 0:N-1
    % NLP variable for control
    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw; 0];
    ubw = [ubw; Fm];
    w0 = [w0; 0];
    
    % New NLP variable for state at end of interval
    Xkpl = MX.sym(['X_' num2str(k+1)], nx);
    w = [w, {Xkpl}];
    lbw = [lbw; 0; 0; mf];
    ubw = [ubw; inf; inf; m0];
    w0 = [w0; 0; 0; 0];
    
    % Add equality constraints
    g = [g, {dt*f(Xk, Uk) + Xk - Xkpl}];
    lbg = [lbg; 0; 0; 0];
    ubg = [ubg; 0; 0; 0];
    
    % Update X_{k+1}
    Xk = Xkpl;
end

% Set terminal conditions
ubw(end) = mf;
w0(end) = mf;

% Cost function
J = J + (-Xk(2));

dJdw = gradient(J, vertcat(w{:}));
dGdw = jacobian(vertcat(g{:}), vertcat(w{:}));

% Create NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, ...
            'lbx', lbw, 'ubx', ubw, ...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

%% Plot the solution

T_opt = w_opt(1);
x1_opt = w_opt(2:4:end);
x2_opt = w_opt(3:4:end);
x3_opt = w_opt(4:4:end);
u_opt = w_opt(5:4:end);
tgrid  = linspace(0, T_opt, N+1);

ocpSol.t = tgrid';
ocpSol.X = [x1_opt x2_opt x3_opt];
ocpSol.U = [u_opt; u_opt(end)];

% if fixed_tf == 0    
%     animateRocket(ocpSol, 'rocketOut_free', 1)
% else
%     animateRocket(ocpSol, 'rocketOut_100', 1)
% end

disp(['Hight: ' num2str(x2_opt(end))])
disp(['Final time: ' num2str(T_opt)])

figure(2);

subplot(411); 
hold on
plot(tgrid, ocpSol.X(:,1))
xlabel('Time')
ylabel('Velocity')

subplot(412);
hold on
plot(tgrid, ocpSol.X(:,2))
xlabel('Time')
ylabel('Height')

subplot(413); 
hold on
plot(tgrid, ocpSol.X(:,3))
xlabel('Time')
ylabel('Mass')

subplot(414); 
hold on
stairs(tgrid, ocpSol.U)
xlabel('Time')
ylabel('F (Control)')

% set(gcf,'Position',[0 0 300 800]);

%% Rocket Model
function [dxdt] = rocket_model(x, u)

% States and control
v = x(1);
h = x(2);
m = x(3);
F = u;

% Parameters
D0 = 0.01227;
beta = 0.145e-3;
c = 2060;
g0 = 9.81;
r0 = 6.371e6;

% Drag and gravity
D = D0*exp(-beta*h);
F_D = sign(v)*D*v^2;
g = g0*(r0/(r0+h))^2;

% Dynamics
dv = (F*c - F_D)/m - g;
dh = v;
dm = -F;
dxdt = [dv; dh; dm];

end

