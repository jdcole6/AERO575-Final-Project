function [M,D,DM,DD] = MatSet(x)
%% Declare Function Variables
mu = 132712*10^6; % solar gravitational parameter [Km^3/s^2]
syms p(t) f(t) g(t) h(t) k(t) L(t) M D
q = [p f g h k L];
ww = 1 + f*cos(L) + g*sin(L);
s2 = 1 + h^2 + k^2;

%% Setup Matrices
% M = [ 0,                 ((2*p)/ww)*sqrt(p/mu),               0;
%       sqrt(p/mu)*sin(L), sqrt(p/mu)*((ww + 1)*cos(L) + f)/ww, -sqrt(p/mu)*(h*sin(L) - k*cos(L))*g/ww;
%      -sqrt(p/mu)*cos(L), sqrt(p/mu)*((ww + 1)*sin(L) + g)/ww,  sqrt(p/mu)*(h*sin(L) - k*cos(L))*f/ww; 
%       0,                 0,                                  sqrt(p/mu)*s2/(2*ww)*cos(L);
%       0,                 0,                                  sqrt(p/mu)*s2/(2*ww)*sin(L);
%       0,                 0,                                  sqrt(p/mu)*(h*sin(L) - k*cos(L))/ww];
%D = [0 0 0 0 0 sqrt(mu*p)*(w/p)^2]'
D(1) = 0;
D(2) = 0;
D(3) = 0;
D(4) = 0;
D(5) = 0;
D(6) = (mu*p)^.5 * (ww/p)^2;
M(1,1)= 0;
M(1,2)= (2*p/ww)*(p/mu)^.5;
M(1,3)= 0;
M(2,1)= (p/mu)^(.5)*sin(L);
M(2,2)= (p/mu)^(.5)*((ww + 1)*cosd(L) + f)/ww;
M(2,3)= -(p/mu)^(.5)*(h*sind(L) - k*cosd(L))*g/ww;
M(3,1)= -(p/mu)^(.5)*cosd(L);
M(3,2)= (p/mu)^(.5)*((ww + 1)*sind(L) + g)/ww;
M(3,3)= (p/mu)^(.5)*(h*sind(L) - k*cosd(L))*f/ww;
M(4,1)= 0;
M(4,2)= 0;
M(4,3)= (p/mu)^(.5)*s2/(2*ww)*cosd(L);
M(5,1)= 0;
M(5,2)= 0;
M(5,3)= (p/mu)^(.5)*s2/(2*ww)*sind(L);
M(6,1)= 0;
M(6,2)= 0;
M(6,3)= (p/mu)^(.5)*(h*sind(L) - k*cosd(L))/ww;

for k1 = 1:size(M,1)
    for k2 = 1:size(M,2)
        dM_dq{k1,k2} = functionalDerivative(M(k1,k2),[p,f,g,h,k,L]);
    end
end
for kk = 1:6
        dD_dq{kk} = functionalDerivative(D(kk),[p,f,g,h,k,L]);
end
%% Establish Variable values
p = x(1);
f = x(2);
g = x(3);
h = x(4);
k = x(5);
L = x(6);

%% Convert symbolics to variable
double(p);
double(f);
double(g);
double(k);
double(h);
double(L);
s2 = subs(s2);
ww = subs(ww);
q = subs(q);

%% Convert Partial differential matrix into numeric and discretize into components
DMDQ = double(subs(cell2sym(dM_dq)));
DDDQ = double(subs(cell2sym(dD_dq)));
dMdp = DMDQ(1:6:end,:);
dMdf = DMDQ(2:6:end,:);
dMdg = DMDQ(3:6:end,:);
dMdh = DMDQ(4:6:end,:);
dMdk = DMDQ(5:6:end,:);
dMdL = DMDQ(6:6:end,:);
dDdp = DDDQ(1:6:end);
dDdf = DDDQ(2:6:end);
dDdg = DDDQ(3:6:end);
dDdh = DDDQ(4:6:end);
dDdk = DDDQ(5:6:end);
dDdL = DDDQ(6:6:end);

%% recompose matrices and derivative
DM = {dMdp;dMdf;dMdg;dMdh;dMdk;dMdL};
DD = {dDdp;dDdf;dDdg;dDdh;dDdk;dDdL};
M = double(subs(M));
D = double(subs(D));
end