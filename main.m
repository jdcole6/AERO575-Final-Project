% Kelsen Case and Jack Cole
% AERO 575
% Final Project

%%
clc; clear;

w = 1 + f*cos(L) + g*sin(L);
s = sqrt(1 + h^2 + h^2);

p = a*(1 - e^2);
f = e*cos(w + omega);
g = e*sin(w + omega);
h = tan(i/2)*cos(omega);
k = tan(i/2)*sin(omega);
L = omega + w + v;

T = (2*eta*P)/(g*Isp);
alpha = [sin(alpha)*cos(beta) cos(alpha)*cos(beta) sin(beta)]';

M = [ 0                 ((2*p)/w)*sqrt(p/mu)               0;
      sqrt(p/mu)*sin(L) sqrt(p/mu)*((w + 1)*cos(L) + f)/w -sqrt(p/mu)*(h*sin(L) - k*cos(L))*g/w;
     -sqrt(p/mu)*cos(L) sqrt(p/mu)*((w + 1)*sin(L) + g)/w  sqrt(p/mu)*(h*sin(L) - k*cos(L))*f/w; 
      0                 0                                  sqrt(p/mu)*s^2/(2*w)*cos(L);
      0                 0                                  sqrt(p/mu)*s^2/(2*w)*sin(L);
      0                 0                                  sqrt(p/mu)*(h*sin(L) - k*cos(L))/w];
  
D = [0 0 0 0 0 sqrt(mu*p)*(w/p)^2]';

xdot = M*(T/m)*alpha + D;
mdot = -T/(g*Isp);
