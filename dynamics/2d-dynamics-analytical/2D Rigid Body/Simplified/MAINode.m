%% ODE for simulation of 2D Spine
% Author: Zeerek Ahmad
% Date: 07/6/2014
% Modified: 7/6/2014
% BEST Lab Berkeley CA

%%

clear all
clc 
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
%I = (0.5^2*0.33)*3; %Approximation
I = 3*(0.5^2);
r0 = 0.5;
k =  2500;
ks = 2500;
% c = 200;
c = 0; %damping coefficient
t= 0;

SpineODE(t,x_input,l1,l2,l3,l,angles,k,I,m,r0,c);

tspan = 0:0.01:5;

%% SIMULATE
options = odeset('reltol',1.e-10,'abstol',1.e-10);
[T,Y] = ode45(@(t,xx) SpineODE(t,xx,l1,l2,l3,l,angles,k,I,m,r0,c),tspan,x_input);

% %%
% plot(T,Y(:,3))

%%
kEnergy = m./2 .* Y(:,4).^2;
kEnergy = kEnergy';
% plot(T,(m./2).*Y(:,4).^2)

%% Build spine segments over time
x = Y(:,1);
y = Y(:,3);
theta = Y(:,5);

L(1).x = A.bars(1).lx;
L(1).y = A.bars(1).ly;
L(2).x = A.bars(2).lx;
L(2).y = A.bars(2).ly;
L(3).x = A.bars(3).lx;
L(3).y = A.bars(3).ly;

for n = 1:size(Y,1)
    L(4).x{n} = [x(n); x(n) + l.*cos(angles(1) - theta(n))];
    L(4).y{n} = [y(n); y(n) + l.*sin(angles(1) - theta(n))];
    L(5).x{n} = [x(n); x(n) + l.*cos(angles(2) - theta(n))];
    L(5).y{n} = [y(n); y(n) + l.*sin(angles(2) - theta(n))];
    L(6).x{n} = [x(n); x(n) + l.*cos(angles(3) - theta(n))];
    L(6).y{n} = [y(n); y(n) + l.*sin(angles(3) - theta(n))];
    
    r = zeros(5,2);
    r(1,:) = [L(4).x{n}(2) - L(1).x(2), L(4).y{n}(2) - L(1).y(2)];
    r(2,:) = [L(5).x{n}(2) - L(2).x(2), L(5).y{n}(2) - L(2).y(2)];
    r(3,:) = [L(6).x{n}(2) - L(3).x(2), L(6).y{n}(2) - L(3).y(2)];
    r(4,:) = [L(6).x{n}(2) - L(1).x(2), L(6).y{n}(2) - L(1).y(2)];
    r(5,:) = [L(6).x{n}(2) - L(2).x(2), L(6).y{n}(2) - L(2).y(2)];
%     keyboard;
    
    for q = 1:5
       R(q) = norm(r(q,:)); 
       pE(q) = 0.5*k.*(R(q) - r0)^2;
    end
    
    pEnergy(n) = sum(pE);
    
end

%% Plot energies
figure;
subplot(3,1,1)
plot(T,pEnergy + kEnergy)
title('Total Energy vs Time')

subplot(3,1,2)
plot(T,pEnergy)
title('Potential Energy vs Time')

subplot(3,1,3)
plot(T,kEnergy)
title('Total Kinetic Energy vs Time');


%% Play the Simulation
figure;

numz = 0;
for i = 2:size(Y,1)
    if round(T(i),4) ~= round(T(i-1),4)
        for n = 1:3
            plot(L(n).x,L(n).y,'b');
            hold on
        end
        
        for n = 4:6
            plot(L(n).x{i},L(n).y{i},'r');
        end
        
        axis([-4 4 -4 4]);
        hold off
        legend(num2str(T(i)))
        numz = numz+1;
        drawnow;
        Frames(numz) = getframe;
    end
end
