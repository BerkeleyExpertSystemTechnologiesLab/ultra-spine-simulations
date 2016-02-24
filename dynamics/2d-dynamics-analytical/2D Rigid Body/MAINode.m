%% ODE for simulation of 2D Spine
% Author: Zeerek Ahmad
% Date: 07/6/2014
% Modified: 7/6/2014
% BEST Lab Berkeley CA


%% Set up initial conditions

%BUILD SPINE
x1 = [0 0];
l = 1;

A = Spine2(l,109.5);

l1 = [A.bars(1).lx(2), A.bars(1).ly(2)];
l2 = [A.bars(2).lx(2), A.bars(2).ly(2)];
l3 = [A.bars(3).lx(2), A.bars(3).ly(2)];

angles = [A.th];

%parameters of x2
x_input(1) = 0; %x
x_input(2) = 0; %xdot

x_input(3) = 2; %y
x_input(4) = 0; %ydot

x_input(5) = 0; %theta
x_input(6) = 0; %thetadot


%CONSTANTS
m = 3;
I = (0.5^2*0.33)*3; %Approximation
r0 = 0.5;
k =  1500;
ks = 1500;
 
t= 0;

SpineODE(t,x_input,l1,l2,l3,l,angles,k,I,m,r0);

tspan = [0 10];

%% SIMULATE
options = odeset('reltol',1.e-10,'abstol',1.e-10);
[T,Y] = ode45(@(t,xx) SpineODE(t,xx,l1,l2,l3,l,angles,k,I,m,r0),tspan,x_input);

%%
plot(T,Y(:,3))

%%
plot(T,(m./2).*Y(:,4).^2)

