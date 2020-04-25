% Kelsen Case and Jack Cole
% AERO 575
% Final Project

%%
clc; clear;

%% Initialize Variable

global AU m_0 P_0 eta c

mu = 132712*10^6; % solar gravitational parameter [Km^3/s^2]
g0 = 9.81;
AU = 1.495978707*10^8; % km

alpha = 0;
beta = 0;

% Direct Earth-Mars
m_0 = 1200; % kg
P_0 = 6.5e3; % W drops at (1/r^2)
Isp = 3100; % s
eta = 0.65;
c = g0*Isp;

% Launch Date = Jan 1 2009
% Transfer Time = 1.5 years
transfer_time_s = 1.5*365*24*60*60;

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

global x_0 tspan

a_0 = a_earth(1);
e_0 = e_earth(1);
w_0 = w_earth(1);
omega_0 = omega_earth(1);
i_0 = i_earth(1);
v_0 = v_earth(1);

p_0 = a_0*(1 - e_0^2);
f_0 = e_0*cos(w_0 + omega_0);
g_0 = e_0*sin(w_0 + omega_0);
h_0 = tan(i_0/2)*cos(omega_0);
k_0 = tan(i_0/2)*sin(omega_0);
L_0 = omega_0 + w_0 + v_0;

tf = ceil(1.5*365);
tf = 2;
tspan = 1:ceil(tf);

x_0 = [p_0;f_0;g_0;h_0;k_0;L_0;m_0];
p_guess = [-1;0;0;0;0;0;0];

p_0 = fsolve(@contraints, p_guess);

%% Plotting

x_earth = a_earth(tspan).*cosd(v_earth(tspan))/AU;
y_earth = a_earth(tspan).*sind(v_earth(tspan))/AU;

x_earth_final = a_earth(tspan(end)).*cosd(v_earth(tspan(end)))/AU;
y_earth_final = a_earth(tspan(end)).*sind(v_earth(tspan(end)))/AU;

x_mars = a_mars(tspan).*cosd(v_mars(tspan))/AU;
y_mars = a_mars(tspan).*sind(v_mars(tspan))/AU;

x_mars_final = a_mars(tspan(end)).*cosd(v_mars(tspan(end)))/AU;
y_mars_final = a_mars(tspan(end)).*sind(v_mars(tspan(end)))/AU;

sun_pos = [0,0];

axis_num = 2;

figure(1)
plot(sun_pos(1), sun_pos(2), 'ro');
text(sun_pos(1), sun_pos(2) - axis_num*0.1, 'Sun');
hold on
plot(x_earth,y_earth, '--k');
plot(x_earth_final,y_earth_final, 'go');
text(x_earth_final,y_earth_final - axis_num*0.1, 'Earth');
plot(x_mars,y_mars, '--k');
plot(x_mars_final,y_mars_final, 'o', 'Color', [0.8500 0.3250 0.0980]);
text(x_mars_final,y_mars_final - axis_num*0.1, 'Mars');
axis([-axis_num axis_num -axis_num axis_num]);
xlabel('x, AU');
ylabel('y, AU');
axis square

%% Functions

function f = contraints(p_0)

global x_0

[xp] = RungeKutta(@sys, [x_0;p_0]);
x = xp(end,:);

p = x(1);
f = x(2);
g = x(3);
h = x(4);
k = x(5);
L = x(6);

f(1) = p - a*(1 - e^2);
f(2) = f - e*cos(w + omega);
f(3) = g - e*sin(w + omega);
f(4) = h - tan(i/2)*cos(omega);
f(5) = k - tan(i/2)*sin(omega);
f(6) = L - omega + w + v;

end

function [state_dot] = sys(t, xp)

global P_0 eta c

x = xp(1:6);
m = xp(7);
lam = xp(8:13);
lam_m = xp(14);

% M = [ 0                 ((2*p)/ww)*sqrt(p/mu)                0;
%       sqrt(p/mu)*sin(L) sqrt(p/mu)*((ww + 1)*cos(L) + f)/ww -sqrt(p/mu)*(h*sin(L) - k*cos(L))*g/ww;
%      -sqrt(p/mu)*cos(L) sqrt(p/mu)*((ww + 1)*sin(L) + g)/ww  sqrt(p/mu)*(h*sin(L) - k*cos(L))*f/ww; 
%       0                 0                                    sqrt(p/mu)*s^2/(2*ww)*cos(L);
%       0                 0                                    sqrt(p/mu)*s^2/(2*ww)*sin(L);
%       0                 0                                    sqrt(p/mu)*(h*sin(L) - k*cos(L))/ww];
%   
% D = [0 0 0 0 0 sqrt(mu*p)*(ww/p)^2]';

[M,D,DM,DD] = MatSet(x);

alpha_vec = -(lam'*M)'/norm(lam'*M);

dMdp = DM{1,1};
dMdf = DM{2,1};
dMdg = DM{3,1};
dMdh = DM{4,1};
dMdk = DM{5,1};
dMdL = DM{6,1};

dDdp = DD{1,1}';
dDdf = DD{2,1}';
dDdg = DD{3,1}';
dDdh = DD{4,1}';
dDdk = DD{5,1}';
dDdL = DD{6,1}';

% P = P_0/r^2;

T = (2*eta*P_0)/c;

xdot = M*(T/m)*alpha_vec + D';
mdot = -T/c;
lam_pdot = -(lam'*dMdp*(T/m)*alpha_vec + lam'*dDdp);
lam_fdot = -(lam'*dMdf*(T/m)*alpha_vec + lam'*dDdf);
lam_gdot = -(lam'*dMdg*(T/m)*alpha_vec + lam'*dDdg);
lam_hdot = -(lam'*dMdh*(T/m)*alpha_vec + lam'*dDdh);
lam_kdot = -(lam'*dMdk*(T/m)*alpha_vec + lam'*dDdk);
lam_ldot = -(lam'*dMdL*(T/m)*alpha_vec + lam'*dDdL);
lam_mdot = -norm(lam'*M)*(T/m^2);

lam_xdot = [lam_pdot; lam_fdot; lam_gdot; lam_hdot; lam_kdot; lam_ldot];

state_dot = [xdot;mdot;lam_xdot;lam_mdot];

end

function [y] = RungeKutta(func, y0)

h = 1;                                            % step size
x = 0:h:10;                                         
y = zeros(14,length(x)); 
y(:,1) = y0;                                      % initial condition

for i = 1:(length(x)-1)                           % calculation loop
    k_1 = func(y(:,i));
    k_2 = func(y(:,i)+0.5*h*k_1);
    k_3 = func(y(:,i)+0.5*h*k_2);
    k_4 = func(y(:,i)+k_3*h);

    y(:,i+1) = y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end

end