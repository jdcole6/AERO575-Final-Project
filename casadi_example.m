addpath('C:\Users\Owner\Documents\MATLAB\AERO 575\Final Project\Casadi Files')
import casadi.*
x = MX.sym('x')
disp(jacobian(sin(x),x))