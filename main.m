% Kelsen Case and Jack Cole
% AERO 575
% Final Project

%%
clc; clear;

%% Initialize Variables

global a e omega w v i mu eta c 

a = 1000;
e = 0.8;
omega = 15;
w = 30;
v = 30;
i = 5;
g0 = 9.81;
mu = 3.98e16;

alpha0 = 0;
beta0 = 0;

% Direct Earth-Mars
m0 = 1200; % kg
P0 = 6.5e3; % W drops at (1/r^2)
Isp = 3100; % s
eta = 0.65;
c = g0*Isp;

% Launch Date = Jan 1 2009
% Transfer Time = 1.5 years
transfer_time = 1.5*365*24*60*60;

%% Main

beta = [alpha beta];
phi = [lambda(t0) beta];

J = phi;

%% Functions

function [xdot, mdot] = clp(t,x)

global a e omega w v i eta P c mu

p = a*(1 - e^2);
f = e*cos(w + omega);
g = e*sin(w + omega);
h = tan(i/2)*cos(omega);
k = tan(i/2)*sin(omega);
L = omega + w + v;

T = (2*eta*P)/c;
alpha = [sin(alpha)*cos(beta) cos(alpha)*cos(beta) sin(beta)]';

w = 1 + f*cos(L) + g*sin(L);
s = sqrt(1 + h^2 + k^2);

M = [ 0                 ((2*p)/w)*sqrt(p/mu)               0;
      sqrt(p/mu)*sin(L) sqrt(p/mu)*((w + 1)*cos(L) + f)/w -sqrt(p/mu)*(h*sin(L) - k*cos(L))*g/w;
     -sqrt(p/mu)*cos(L) sqrt(p/mu)*((w + 1)*sin(L) + g)/w  sqrt(p/mu)*(h*sin(L) - k*cos(L))*f/w; 
      0                 0                                  sqrt(p/mu)*s^2/(2*w)*cos(L);
      0                 0                                  sqrt(p/mu)*s^2/(2*w)*sin(L);
      0                 0                                  sqrt(p/mu)*(h*sin(L) - k*cos(L))/w];
  
D = [0 0 0 0 0 sqrt(mu*p)*(w/p)^2]';

xdot = M*(T/m)*alpha + D;
mdot = -T/c;

end
