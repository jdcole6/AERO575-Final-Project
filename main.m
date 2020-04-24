% Kelsen Case and Jack Cole
% AERO 575
% Final Project

%%
clc; clear;

%% Initialize Variables

global a e omega w v i mu eta c alpha beta

a = 1000;
e = 0.8;
omega = 15;
w = 30;
v = 30;
i = 5;
g0 = 9.81;
mu = 3.98e16;

alpha = 0;
beta = 0;

% Direct Earth-Mars
m0 = 1200; % kg
P0 = 6.5e3; % W drops at (1/r^2)
Isp = 3100; % s
eta = 0.65;
c = g0*Isp;

% Launch Date = Jan 1 2009
% Transfer Time = 1.5 years
transfer_time = 1.5*365*24*60*60;

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

earth = fopen('earth_ephemeris.txt');
earth_mat = textscan(earth, '%f %s %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'delimiter', ',', 'CollectOutput', true);

mars = fopen('mars_ephemeris.txt');
mars_mat = textscan(mars, '%f %s %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'delimiter', ',', 'CollectOutput', true);

e_earth = earth_mat{1,3}(:,1);
i_earth = earth_mat{1,3}(:,3);
omega_earth = earth_mat{1,3}(:,4);
w_earth = earth_mat{1,3}(:,5);
v_earth = earth_mat{1,3}(:,9);
a_earth = earth_mat{1,3}(:,10);

e_mars = mars_mat{1,3}(:,1);
i_mars = mars_mat{1,3}(:,3);
omega_mars = mars_mat{1,3}(:,4);
w_mars = mars_mat{1,3}(:,5);
v_mars = mars_mat{1,3}(:,9);
a_mars = mars_mat{1,3}(:,10);

%% Main

beta = [alpha beta];
phi = [lambda(t0) beta];

J = phi;

%% Plotting

x = a_earth;

%% Functions

function f = contraints(p_0)

global a e omega w v i eta P c mu

[xp] = RungeKutta(sys,[0 t1], [x0;p_0]);
xx = xp(end,:);

f(1) = p - a*(1 - e^2);
f(2) = f - e*cos(w + omega);
f(3) = g - e*sin(w + omega);
f(4) = h - tan(i/2)*cos(omega);
f(5) = k - tan(i/2)*sin(omega);
f(6) = L - omega + w + v;
f(7) = T - (2*eta*P)/c;
f(8) = ww - 1 + f*cos(L) + g*sin(L);
f(9) = s - sqrt(1 + h^2 + k^2);

end

function [xdot, mdot] = sys(t,x)

global mu

alpha_vec = [sin(alpha)*cos(beta) cos(alpha)*cos(beta) sin(beta)]';

% M = [ 0                 ((2*p)/ww)*sqrt(p/mu)                ;
%       sqrt(p/mu)*sin(L) sqrt(p/mu)*((ww + 1)*cos(L) + f)/w  -sqrt(p/mu)*(h*sin(L) - k*cos(L))*g/ww;
%      -sqrt(p/mu)*cos(L) sqrt(p/mu)*((ww + 1)*sin(L) + g)/w  sqrt(p/mu)*(h*sin(L) - k*cos(L))*f/ww; 
%       0                 0                                   sqrt(p/mu)*s^2/(2*ww)*cos(L);
%       0                 0                                   sqrt(p/mu)*s^2/(2*ww)*sin(L);
%       0                 0                                   sqrt(p/mu)*(h*sin(L) - k*cos(L))/w];
%   
% D = [0 0 0 0 0 sqrt(mu*p)*(w/p)^2]';

x = [1 1 1 1 1 1];
[M,D,DM,DD] = MatSet(x)

xdot = M*(T/m)*alpha_vec + D;
mdot = -T/c;
lam_
lam_xdot = -(lam'*dMdx*(T/m)*alpha_vec + lam'*dDdx);
lam_mdot = -norm(lam'*M)*(T/m^2);

end

function [y] = RungeKutta(func, tspan, y0)

h=1;                                               % step size
x = 0:h:5;                                         % Calculates upto y(3)
y = zeros(1,length(x)); 
y(1) = y0;                                         % initial condition

for i = tspan(1):tspan(2)                          % calculation loop
    k_1 = func(x(i),y(i));
    k_2 = func(x(i)+0.5*h,y(i)+0.5*h*k_1);
    k_3 = func((x(i)+0.5*h),(y(i)+0.5*h*k_2));
    k_4 = func((x(i)+h),(y(i)+k_3*h));

    y(i+1) = y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end

end