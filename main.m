% Kelsen Case and Jack Cole
% AERO 575
% Final Project

%%
clc; clear;

%% Initialize Variables

global AU m_0 P_0 eta c e_mars i_mars omega_mars w_mars v_mars a_mars

mu = 132712*10^6; % solar gravitational parameter [Km^3/s^2]
g0 = 9.81; % Earth sea-level gravity [m/s^2]
AU = 1.495978707*10^8; % AU conversion [km]

alpha = 0;
beta = 0;

% Direct Earth-Mars parameters
m_0 = 1200; % initial spacecraft mass [kg]
P_0 = 6.5e3; % initial spacecraft power [W] drops at (1/r^2)
Isp = 3100; % specific impulse of engine [s]
eta = 0.65; % efficiency of engine
c = g0*Isp; % exhaust velocity of engine [m/s]

% Launch Date = Jan 1 2009
% Transfer Time = 1.5 years
transfer_time_s = 1.5*365*24*60*60;

% Parameters in ephemeris data
%      JDTDB    Julian Day Number, Barycentric Dynamical Time
% 1      EC     Eccentricity, e                                                   
% 2      QR     Periapsis distance, q (km)                                        
% 3      IN     Inclination w.r.t XY-plane, i (degrees)                           
% 4      OM     Longitude of Ascending Node, OMEGA, (degrees)                     
% 5      W      Argument of Perifocus, w (degrees)                                
% 6      Tp     Time of periapsis (Julian Day Number)                             
% 7      N      Mean motion, n (degrees/sec)                                      
% 8      MA     Mean anomaly, M (degrees)                                         
% 9      TA     True anomaly, nu (degrees)                                        
% 10     A      Semi-major axis, a (km)                                           
% 11     AD     Apoapsis distance (km)                                            
% 12     PR     Sidereal orbit period (sec)   

% Open Earth ephermis data sheet and read to matrix
earth = fopen('earth_ephemeris.txt');
earth_mat = textscan(earth, '%f %s %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'delimiter', ',', 'CollectOutput', true);

% Open Mars ephermis data sheet and read to matrix
mars = fopen('mars_ephemeris.txt');
mars_mat = textscan(mars, '%f %s %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'delimiter', ',', 'CollectOutput', true);

% Decompose Earth elements into separate vectors
e_earth = earth_mat{1,3}(:,1);
i_earth = earth_mat{1,3}(:,3);
omega_earth = earth_mat{1,3}(:,4);
w_earth = earth_mat{1,3}(:,5);
v_earth = earth_mat{1,3}(:,9);
a_earth = earth_mat{1,3}(:,10);

% Decompose Mars elements into separate vectors
e_mars = mars_mat{1,3}(:,1);
i_mars = mars_mat{1,3}(:,3);
omega_mars = mars_mat{1,3}(:,4);
w_mars = mars_mat{1,3}(:,5);
v_mars = mars_mat{1,3}(:,9);
a_mars = mars_mat{1,3}(:,10);

%% Main

global x_0 tspan

% Define initial conditions for spacecraft based on Earth position
a_0 = a_earth(1);
e_0 = e_earth(1);
w_0 = w_earth(1);
omega_0 = omega_earth(1);
i_0 = i_earth(1);
v_0 = v_earth(1);

% Convert orbital elements into equinoctial elements
p_0 = a_0*(1 - e_0^2);
f_0 = e_0*cos(w_0 + omega_0);
g_0 = e_0*sin(w_0 + omega_0);
h_0 = tan(i_0/2)*cos(omega_0);
k_0 = tan(i_0/2)*sin(omega_0);
L_0 = omega_0 + w_0 + v_0;

% Define final time and time span
tf = ceil(1.5*365);
tf = 2;
tspan = 1:ceil(tf);

% Compose state vector and co-state guesses
x_0 = [p_0;f_0;g_0;h_0;k_0;L_0;m_0];
lam_guess = [-1;0;0;0;0;0;0];

% Fsolve function call
% lam_0 = fsolve(@contraints, lam_guess);

% Define elements for fmincon function call
options = optimoptions(@fmincon,'Algorithm','sqp'); % use SQP algorithm
Obj = @(x) -x(7); % Define cost function to maximize mass
transfer = 1.5; % transfer time guess [years]
date = 09206; % 25-Jul-2009, assuming we are using Julian dates
G = [lam_guess;transfer;date]; % these are what I'm guessing are the SQP variables they are talking about
[lambda] = fmincon(Obj,G,[],[],[],[],[],[],@traj,options);

%% Plotting

% Deterime Earth position in x,y coordinates
x_earth = a_earth(tspan).*cosd(v_earth(tspan))/AU;
y_earth = a_earth(tspan).*sind(v_earth(tspan))/AU;

% Deterime Mars position in x,y coordinates
x_mars = a_mars(tspan).*cosd(v_mars(tspan))/AU;
y_mars = a_mars(tspan).*sind(v_mars(tspan))/AU;

% Sun position is at origin
sun_pos = [0,0];

axis_lim = 2; % axis limit of 2 AU

% Define spacecraft position in terms of equinoctial elements
r_sc = zeros(2,1);

alph = sqrt(h^2 - k^2);
s = sqrt(1 + h^2 + k^2);
w = 1 + f*cosd(L) + g*sind(L);
r = p/w;

% Spacecraft x,y coordinates
r_sc(1) = (r/s^2)*(cosd(L) + alph^2*cosd(L) + 2*h*k*sind(L));
r_sc(2) = (r/s^2)*(sind(L) + alph^2*sind(L) + 2*h*k*cosd(L));

figure(1)
plot(sun_pos(1),sun_pos(2), 'ro'); % plot Sun as red circle
text(sun_pos(1),sun_pos(2) - axis_lim*0.1, 'Sun'); % label Sun
hold on
plot(x_earth,y_earth, '--k'); % plot Earth trajectory
plot(x_earth(end),y_earth(end), 'go'); % plot Earth as green circle
text(x_earth(end),y_earth(end) - axis_lim*0.1, 'Earth'); % label Earth
plot(x_mars,y_mars, '--k'); % plot Mars trajectory
plot(x_mars(end),y_mars(end), 'o', 'Color', [0.8500 0.3250 0.0980]); % plot Mars as orange circle
text(x_mars(end),y_mars(end) - axis_lim*0.1, 'Mars'); % label Mars

plot(r_sc,r_sc, '-k'); % plot spacecraft trajectory
plot(r_sc(end),r_sc(end), 'bo'); % plot spacecraft as blue circle
text(r_sc(end),r_sc(end) - axis_lim*0.1, 'Spacecraft'); % label spacecraft

% Set axis limits and label axes
axis([-axis_lim axis_lim -axis_lim axis_lim]);
xlabel('x, AU');
ylabel('y, AU');
axis square

%% Functions

function [c,ceq] = traj(X)

global x_0 tspan e_mars i_mars omega_mars w_mars v_mars a_mars

lam = X(1:7); % co-state vector
tt = X(8)*365; % transfer time in [days]
start = X(9); % start date
XX = [x_0;lam]; % initial state vector for integration

[X] = RungeKutta(@sys, XX); % calculate trajectory

xx = X(:, end); % pull out states
tf = ceil(tspan(end)); % pull out final time

% Convert equinoctial elements to orbital elements
a_sc = xx(1)/(1-xx(2)^2 - xx(3)^2);
e_sc = sqrt(xx(2)^2 + xx(3)^2);
i_sc = (180/pi())*(2*inv(tan(sqrt(xx(4)^2 + xx(5)^2))));
w_sc = (180/pi())*(inv(tan(xx(3)/xx(2))) - inv(tan(xx(5)/xx(4))));
omega_sc = (180/pi())*inv(tan(xx(5)));
v_sc = (180/pi())*xx(6) - (omega_sc + w_sc);

c = [];

% Define equality constraints on equinoctial elements and terminal conditions
ceq(1) = xx(1) - a_sc*(1 - e_sc^2);
ceq(2) = xx(2) - e_sc*cos(w_sc + omega_sc);
ceq(3) = xx(3) - e_sc*sin(w_sc + omega_sc);
ceq(4) = xx(4) - tan(i_sc/2)*cos(omega_sc);
ceq(5) = xx(5) - tan(i_sc/2)*sin(omega_sc);
ceq(6) = xx(6) - omega_sc + w_sc + v_sc;
ceq(7) = a_sc - a_mars(tf);
ceq(8) = e_sc - e_mars(tf);
ceq(9) = i_sc - i_mars(tf);
ceq(10) = w_sc - w_mars(tf);
ceq(11) = omega_sc - omega_mars(tf);
ceq(12) = v_sc - v_mars(tf);

end

function f = contraints(p_0)

global x_0

[xp] = RungeKutta(@sys, [x_0;p_0]); % calculate trajectory
x = xp(end,:); % pull out states

% Decompose states into equinoctial elements
p = x(1);
f = x(2);
g = x(3);
h = x(4);
k = x(5);
L = x(6);

% Define equinoctial element constraints
f(1) = p - a*(1 - e^2);
f(2) = f - e*cos(w + omega);
f(3) = g - e*sin(w + omega);
f(4) = h - tan(i/2)*cos(omega);
f(5) = k - tan(i/2)*sin(omega);
f(6) = L - omega + w + v;

end

function [state_dot] = sys(xp)

global P_0 eta c

x = xp(1:6); % pull out states
m = xp(7); % pull out mass state
lam = xp(8:13); % pull out co-states
lam_m = xp(14); % pull out mass co-state

% Function call to get matrices
[M,D,DMDQ,DDDQ] = MatSet(x);

alpha_vec = -(lam'*M)'/norm(lam'*M); % define optimal control vector

% P = P_0/r^2;

T_max = (2*eta*P_0)/c; % calculate thrust

H = lam'*M*(T_max/m)*alpha_vec + lam'*D' - lam_m*(T_max/c); % define Hamiltonian

% Define switching funciton
if H < 0
    T = 0;
else
    T = T_max;
end

% Equations of motion for equinoctial elements
xdot = M*(T/m)*alpha_vec + D';
mdot = -T/c;

% Adjoint equations for co-states
lam_xdot = zeros(6,1);

for i = 1:6
    lam_xdot(i) = -(lam'*DMDQ(i:6:end,:)*(T/m)*alpha_vec + lam'*DDDQ(i:6:end)');
end

lam_mdot = -norm(lam'*M)*(T/m^2);

state_dot = [xdot;mdot;lam_xdot;lam_mdot]; % formulate differentials into vector

end

function [y] = RungeKutta(func, y0)

h = 1; % step size
x = 0:h:10; % number of interations                             
y = zeros(length(y0),length(x));
y(:,1) = y0;  % initial condition

for i = 1:(length(x)-1) % integation loop
    k_1 = func(y(:,i));
    k_2 = func(y(:,i)+0.5*h*k_1);
    k_3 = func(y(:,i)+0.5*h*k_2);
    k_4 = func(y(:,i)+k_3*h);

    y(:,i+1) = y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end

end
