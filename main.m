% Kelsen Case and Jack Cole
% AERO 575
% Final Project

%%
clc; clear;

%% Initialize Variables

global mu eta c alpha beta

mu = 132712*10^6; % solar gravitational parameter [Km^3/s^2]

g0 = 9.81;

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

a = 100;
h = 8000;
t1 = 100;
x0 = [0;0;0;0];
p_guess = [0;0;-1;0];

p_0 = fsolve(@fun, p_guess);

%% Plotting

tspan = 1:240;

x_earth = a_earth(tspan).*cosd(v_earth(tspan));
y_earth = a_earth(tspan).*sind(v_earth(tspan));

x_earth_final = a_earth(tspan(end)).*cosd(v_earth(tspan(end)));
y_earth_final = a_earth(tspan(end)).*sind(v_earth(tspan(end)));

x_mars = a_mars(tspan).*cosd(v_mars(tspan));
y_mars = a_mars(tspan).*sind(v_mars(tspan));

x_mars_final = a_mars(tspan(end)).*cosd(v_mars(tspan(end)));
y_mars_final = a_mars(tspan(end)).*sind(v_mars(tspan(end)));

sun_pos = [0,0];

axis_num = 2.5e8;

figure(1)
plot(sun_pos, 'ro');
text(sun_pos(1), sun_pos(2) - axis_num*0.1, 'Sun');
hold on
plot(x_earth,y_earth, '--k');
plot(x_earth_final,y_earth_final, 'go');
text(x_earth_final,y_earth_final - axis_num*0.1, 'Earth');
plot(x_mars,y_mars, '--k');
plot(x_mars_final,y_mars_final, 'o', 'Color', [0.8500 0.3250 0.0980]);
text(x_mars_final,y_mars_final - axis_num*0.1, 'Mars');
axis([-axis_num axis_num -axis_num axis_num]);

%% Functions

function f = contraints(x_0,p_0)
[xp] = RungeKutta(sys,[0 t1],[x0;p_0]);
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
% f(7) = T - (2*eta*P)/c;
% f(8) = ww - 1 + f*cos(L) + g*sin(L);
% f(9) = s - sqrt(1 + h^2 + k^2);

end

function [xdot, mdot, lam_xdot, lam_mdot] = sys(t,x)

x = x(1:6);
m = x(7);
lam = x(8:13);
lam_m = x(14);

% M = [ 0                 ((2*p)/ww)*sqrt(p/mu)               0;
%       sqrt(p/mu)*sin(L) sqrt(p/mu)*((ww + 1)*cos(L) + f)/w -sqrt(p/mu)*(h*sin(L) - k*cos(L))*g/ww;
%      -sqrt(p/mu)*cos(L) sqrt(p/mu)*((ww + 1)*sin(L) + g)/w  sqrt(p/mu)*(h*sin(L) - k*cos(L))*f/ww; 
%       0                 0                                   sqrt(p/mu)*s^2/(2*ww)*cos(L);
%       0                 0                                   sqrt(p/mu)*s^2/(2*ww)*sin(L);
%       0                 0                                   sqrt(p/mu)*(h*sin(L) - k*cos(L))/w];
%   
% D = [0 0 0 0 0 sqrt(mu*p)*(w/p)^2]';

[M,D,DM,DD] = MatSet(x);

alpha_vec = -(lam'*M)'/norm(lam'*M);

dMdp = DM{1,1};
dMdf = DM{2,1};
dMdg = DM{3,1};
dMdh = DM{4,1};
dMdk = DM{5,1};
dMdL = DM{6,1};

dDdp = DD{1,1};
dDdf = DD{2,1};
dDdg = DD{3,1};
dDdh = DD{4,1};
dDdk = DD{5,1};
dDdL = DD{6,1};

P = P0/r^2;

T = (2*eta*P)/c;

xdot = M*(T/m)*alpha_vec + D;
mdot = -T/c;
lam_pdot = -(lam'*dMdp*(T/m)*alpha_vec + lam'*dDdp);
lam_fdot = -(lam'*dMdf*(T/m)*alpha_vec + lam'*dDdf);
lam_gdot = -(lam'*dMdg*(T/m)*alpha_vec + lam'*dDdg);
lam_hdot = -(lam'*dMdh*(T/m)*alpha_vec + lam'*dDdh);
lam_kdot = -(lam'*dMdk*(T/m)*alpha_vec + lam'*dDdk);
lam_ldot = -(lam'*dMdL*(T/m)*alpha_vec + lam'*dDdL);
lam_mdot = -norm(lam'*M)*(T/m^2);

lam_xdot = [lam_pdot; lam_fdot; lam_gdot; lam_hdot; lam_kdot; lam_ldot];

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