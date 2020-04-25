options = optimoptions(@fmincon,'Algorithm','sqp');
Obj = @(x) -x(7);
lambda0 = [-1;0;0;0;0;0];
transfer = 1.5; %[years], transfer time guess
date = 09206; %25-Jul-2009, assuming we are using Julian dates
G = [lambda0;transfer;date]; %these are what I'm guessing are the SQP variables they are talking about
[costates] = fmincon(Obj,G,[],[],[],[],[],[],@Traj,options)




%% functions
function [c,ceq] = Traj(X)
global x0 m0 
lm = 1; %mass costate
l = [G(1:6);lm]; %costate vector
tt = G(7)*365; %transfer time in [days]
Start = G(8); %start date
XX = [x0;m0;l;lm]; %initial state vector for integration

[t, X] = ode45(@INT,[0 tt],XX) % calculate trajectory
xx = X(end,:)
tf = ciel(t(end));
af = xx(1)/(1-xx(2)^2 - xx(3)^2);
ef = sqrt(xx(2)^2 + xx(3)^2);
iif = (180/pi())*(2*inv(tan(sqrt(xx(4)^2 + xx(5)^2))));
wf = (180/pi())*(inv(tan(xx(3)/xx(2)))-inv(tan(xx(5)/xx(4))));
omegaf = (180/pi())*inv(tan(xx(5),xx(4)));
thetaf = (180/pi())*xx(6)-(omegaf+wf);

ceq(1) = af - a_mars(tf);
ceq(2) = ef - e_mars(tf);
ceq(3) = iif - i_mars(tf);
ceq(4) = wf - w_mars(tf);
ceq(5) = omegaf - omega_mars(tf);
ceq(6) = thetaf - theta_mars(tf);
end



function dxdt = INT(t,xml) %Primary ODE solver with unknown costate values
x = xml(1:6); % state vector
ll= xml(8:end-1); % costate vector
Lm = xml(end); % mass costate
m = xml(7); %mass
P = p0; %assume constant power for now
T = 2*eta*P/(g*Isp);
[M,D,DMDQ,DDDQ] = MatSet(x);
% determine derivatives
dx = M*(T/m)*alpha+D;
dm = -T/(g*Isp);
for i = 1:6
dl(i) = -(ll'*DMDQ(i:6:end,:)*(T/m)*alpha+ll'*DDDQ(i:6:end);
end
dlm = -norm(ll'*M)*(T/m^2);
%Recompose vecotr
dxdt = [dx;dm;dl;dlm]; %recombinging states,mass,costates
end
